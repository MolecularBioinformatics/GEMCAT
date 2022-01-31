import cobra
from numpy.lib.index_tricks import fill_diagonal
import pandas as pd
import numpy as np
import warnings
from typing import Tuple, Union, List, Dict
import abc
import re
from pyreporter import utils, verification, regexes

regex_recon3d_old = regexes.REGEX['Recon3D_old']

class Expression(abc.ABC):
    """
    Abstract base class for Expression to define interface. 
    Do NOT use.
    """
    @abc.abstractmethod
    def __init__(
        self,
        data: pd.Series,
        gpr: Dict[str, List[str]],
    ):
        pass
    
    @staticmethod
    def _verify_data(data):
        verification.raise_for_duplicated_index(data)
        verification.raise_for_non_int_float_64_dtype(data)
        data = verification.convert_df_index_to_str(data)
        return data

    @abc.abstractmethod
    def load_gpr():
        pass

    @abc.abstractmethod
    def map():
        pass

    def get_mapped_values(self) -> np.array:
        """
        Returns the mapped values as a NumPy array.
        :return: Gene expression values
        :rtype: 1-D np.array (r)
        """
        return self.mapped_values.values

    def fillna_mean(self):
        """
        Fills all NaN values in self.mapped_values with the mean
        of the other data.
        """
        mean = self.mapped_values.mean()
        self.mapped_values.fillna(mean, inplace=True)


class ExpressionFang2012(Expression):
    """
    Integration of gene expression according to Fang et al. (2012)
    doi.org/10.1371/journal.pcbi.1002688
    """
    

    def __init__(
        self,
        gpr: dict[str, str],
        reaction_gene_mapping: Dict[str, List[str]],
        data: pd.Series,
        gene_fill = 0.,
        ):
        """
        
        :param model: [description]
        :type model: cobra.Model
        :param data: [description]
        :type data: pd.Series
        :param re_gene: [description], defaults to '\d+\.\d'
        :type re_gene: str, optional
        """
        self.quant_gpr: pd.Series = None
        self.mapped_values: pd.Series = None
        self.gene_fill = gene_fill
        self.data = self._verify_data(data)
        self.load_gpr(gpr)
        self.reaction_gene_list = reaction_gene_mapping
        self.rewrite_gpr()
        self.map()
        self.fillna_mean()

    def load_gpr(self, gpr):
        self.gpr = pd.Series(gpr)
    
    def rewrite_gpr(self):
        quant_gpr = self.gpr.index.map(
            lambda rxn: self.rewrite_single_gpr(rxn) 
            )
        self.quant_gpr = pd.Series(list(quant_gpr), self.gpr.index)
        
    def rewrite_single_gpr(self, rxn):
        gpr = self.gpr[rxn]
        genes = self.reaction_gene_list.get(rxn, '')
        if not gpr:
            return ''
        if not isinstance(gpr, str):
            return ''
        
        for gene in genes:
            if isinstance(gene, float):
                gene = str(gene)
            gene_val = float(self.data.get(gene, self.gene_fill))
            gpr = gpr.replace(gene, f'{gene_val}')
        
        gpr = gpr.replace('or', '+')
        re_float = '\d*\.\d*'
        re_and = re.compile(f'{re_float}(?: and {re_float})+')
        hits = re_and.findall(gpr)
        if not hits:
            return gpr
        for hit in hits:
            gids = hit.split('and')
            gids = [x.strip() for x in gids]
            gids = ', '.join(gids)
            new = f'utils.geometric_mean({gids})'
            gpr = gpr.replace(hit, new)
        
        return gpr

    def map(self):
        self.mapped_values = self.quant_gpr.map(
            lambda gpr: self.eval_single_gpr(gpr)
            )

    @staticmethod
    def eval_single_gpr(
        gpr: str, 
        ):
        if len(gpr) == 0:
            return np.nan
        try: 
            result = eval(gpr)
        except SyntaxError:
            print(f'Failed: {gpr}')
            result = np.nan
        return result


class ExpressionMapSingleAverage(Expression):
    """
    Simple integration of expression data by mapping each reaction
    onto the mean expression of the genes involved in it.
    """
    def __init__(
        self, 
        data: pd.Series, 
        gpr: Dict[str, List[str]],
        ):
        """
        Initialize object and map expression data onto reactions.
        :param data: Gene expression data
        :type data: pd.Series
        :param gpr: Gene product rules
        :type gpr: Dict[str, List[str]]
        """
        self.gpr = None
        self.load_gpr(gpr)
        self.mapped_values = None
        data = self._verify_data(data)
        self.map(data)

    def load_gpr(
        self, 
        gpr: Dict[str, List[str]],
        ):
        """
        Load gene product rule into the object.
        :param gpr: [description]
        :type gpr: Dict[str, List[str]]
        """
        self.gpr = gpr

    def map(self, data: pd.Series):
        """
        Map the expression data using the gene product rule.
        Save mapped values in object as a pandas series.
        :param data: Gene expression data
        :type data: pd.Series
        :raises ValueError: Raised if no GPR is set on the object
        """
        if not self.gpr:
            raise ValueError('GPR not yet set')
        vals_dict = {}
        for reaction, genes in self.gpr.items():
            expression = []
            for gene in genes:
                try:
                    expression.append(data.loc[gene])
                except KeyError:
                    expression.append(0)
            if not expression:
                vals_dict[reaction] = np.nan
            else:
                vals_dict[reaction] = sum(expression) / len(expression)
        self.mapped_values = pd.Series(vals_dict)
        self.fillna_mean()


def read_simple_gpr_from_cobra(
    model: cobra.Model,
    ) -> Dict[str, List[str]]:
    """
    Reads a simple version of GPR from a COBRA model:
    All reactions are mapped to all involved genes.
    :return: Dictionary with keys reactions, values associated genes 
    :rtype: Dict[str: List[str]]
    """
    genes = []
    for r in model.reactions:
        genes.append([g.id for g in r.genes])
    
    genes_dict = dict(zip([r.id for r in model.reactions], genes))
    return genes_dict

def read_gpr_strings_from_cobra(
    model : cobra.Model,
) -> Tuple[Dict[str, str], list[str]]:
    """
    Extracts the gene product rules (GPRs) from a model.
    :param model: Model from which to extract GPR
    :type model: cobra.Model
    :return: Tuple of extracted GPR (dict) and list of all genes.
    :rtype: Tuple[ Dict[str, str], list[str] ]
    """
    gpr = {rxn.id: rxn.gene_reaction_rule for rxn in model.reactions}
    genes = {rxn.id: [gene.id for gene in rxn.genes] for rxn in model.reactions}
    return gpr, genes
