import cobra
from numpy.lib.index_tricks import fill_diagonal
import pandas as pd
import numpy as np
import warnings
from typing import Union, List, Dict
import abc
import re
from pyreporter import utils, verification

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
        model: cobra.Model,
        data: pd.Series,
        re_gene = '\d+\.\d',
        ):
        """
        
        :param model: [description]
        :type model: cobra.Model
        :param data: [description]
        :type data: pd.Series
        :param re_gene: [description], defaults to '\d+\.\d'
        :type re_gene: str, optional
        """
        self.re_gene = re_gene
        self.data = self._verify_data(data)
        self.load_gpr(model)
        self.rewrite_gpr()
        self.mapped_values = None
        self.map()
        
    def load_gpr(self, model: cobra.Model):
        rxns = [r.id for r in model.reactions]
        gpr = [r.gene_reaction_rule for r in model.reactions]
        gpr = dict(zip(rxns, gpr))
        self.gpr = pd.Series(gpr)
    
    def rewrite_gpr(self):
        self.gpr = self.gpr.map(
            lambda x: self.rewrite_single_gpr(x, self.re_gene) 
            )
        
    @staticmethod
    def rewrite_single_gpr(gpr, re_gene):
        if not gpr:
            return ''
        if not isinstance(gpr, str):
            return ''
        re_and = re.compile(f'{re_gene}(?: and {re_gene})+')
        gpr = gpr.replace('or', '+')
        res = re_and.findall(gpr)
        if not res:
            return gpr
        for hit in res:
            gids = hit.split('and')
            gids = [x.strip() for x in gids]
            gids = ', '.join(gids)
            new = f'utils.geometric_mean({gids})'
            gpr = gpr.replace(hit, new)
        return gpr

    def map(self):
        fill_val = self.data.mean()
        self.mapped_values = self.gpr.map(
            lambda gpr: self.eval_single_gpr(
                gpr, 
                self.data, 
                fill_val, 
                self.re_gene
                )
            )

    @staticmethod
    def eval_single_gpr(
        gpr: str, 
        gene_vals: pd.Series, 
        fill_val: float, 
        re_gene: str,
        ):
        n_plus = max(gpr.count('+'), 1) 
        # arithmetic mean of alternative enzymes
        # divide by 1 if no plusses
        re_gene = re.compile(re_gene)
        gid_finds = re_gene.findall(gpr)
        for gid in gid_finds:
            gid_val = gene_vals.get(gid, fill_val)
            gpr = gpr.replace(gid, str(gid_val))
        try: 
            result = eval(gpr)
        except SyntaxError:
            print(f'Failed: {gpr}')
            result = fill_val
        return result


class ExpressionModifiedFang2012Single(ExpressionFang2012):
    """
    Integration of gene expression according to Fang et al. (2012)
    doi.org/10.1371/journal.pcbi.1002688
    """
    @staticmethod
    def eval_single_gpr(
        gpr: str, 
        gene_vals: pd.Series, 
        fill_val: float, 
        re_gene: str,
        ):
        # no dividing by number of enzymes as we want absolute vals
        re_gene = re.compile(re_gene)
        gid_finds = re_gene.findall(gpr)
        for gid in gid_finds:
            gid_val = gene_vals.get(gid, fill_val)
            gpr = gpr.replace(gid, str(gid_val))
        try: 
            result = eval(gpr)
        except SyntaxError:
            print(f'Failed: {gpr}')
            result = fill_val
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