# GEMCAT: Gene Expression-based Metabolite Centrality Analyses Tool
A computational toolbox associated with the manuscript entitled "GEMCAT - A new algorithm for gene expression-based prediction of metabolic alterations". 
Cite using: https://www.biorxiv.org/content/10.1101/2024.01.15.575710v1

## Compatibility
The package is tested for compatibility with Python >= 3.10 on Ubuntu and Windows.

## Installation
Install from pip:

```pip install gemcat```

Or clone the repository and install GEMCAT from there using:  

```pip install .```


## Usage

### Standard workflow from the Command-Line Interface (CLI)

Use a single file containing per-gene fold-changes to calculate the resulting differential centralities:
``` gemcat ./expression_file.csv ./model_file.xml -e column_name -a pure -o <result_file.csv>```
Make sure the .csv file is either comma- or tab-delimited.
`column_name` is the name of the column in the file containing the fold-change.

Alternatively, use two files (or one file) with expression values for condition and baseline:
``` gemcat <./condition_file.csv> <./model_file.xml> -e <condition_column_name> -b <./baseline_file> -c <baseline_column_name> -o <result_file.csv>```

Currently only models in XML/SBML format are supported in the CLI.
Further models can used from the Python library.
Support will come to the CLI.

positional arguments:
- expression file path
- model file path

All parameters:
`-i --integration` method to use for expression integration into GPR rules
  - `means`: geometric means for AND, arithmetic means for OR [default, recommended]
  - `average`: simply average across all genes involved in reaction
`-e --expressioncolumn` name of column containing condition expression data
`-r RANKING, --ranking` method to use for ranking 
  - currently only `Pagerank`
`-b BASELINE, --baseline` file containing baseline expression data
`-c BASELINECOLUMN, --baselinecolumn` name of column containing baseline expression data
`-a ADJACENCY, --adjacency` method for adjacency calculation
  - `pure`: interaction between metabolites has unit strength [default, recommended]
  - `half`: interaction between metabolites depends on product stoichiometry
  - `full`: interaction between metabolites depends on product and educt stoichiometry 
`-g GENEFILL, --genefill` value to missing expression entries with
`-v VERBOSE, --verbose` enables verbose output
`-o OUTFILE, --outfile` write output to this file
`-l LOGFILE, --logfile` write logs to this file

### Standard workflow in Python using a CobraPy model
```
import gemcat as gc
results = gc.workflows.workflow_standard(cobra_model: cobra.Model,
                                        mapped_genes_baseline: pd.Series,
                                        mapped_genes_comparison: pd.Series,
                                        adjacency = gc.adjacency_transformation.ATPureAdjacency,
                                        ranking = gc.ranking.PagerankNX,
                                        gene_fill = 1.0
                                        )
```
This will return the changes in centrality relative to the baseline in a Pandas Series.
When using fold-changes as the mapped expression, use a vector of all ones as a comparison.

## Modularity and Configuration
GEMCAT is designed to be modular, and its central components can easily be swapped out or appended by other components 
adhering to the specifications laid out in the module base classes (primarily adjacency transformation, expression integration, and ranking components).
All classes inheriting from the abstract base classes laid out in the modules are swappable.

## Core modules
### Model
The core of the package is the GEMCAT model structure that contains the model data, integrates the workflow, and calculates the results.
### Adjacency
Different approaches can be used to calculate adjacency in the networks.
We offer alternatives and a platform to create custom algorithms for the model.
### Expression
Module covering the mapping of gene values onto reactions in the model via gene product rules.
Providing different algorithms along with a platform to create alternatives.
### Pagerank
Module providing ranking algorithms for the models along with a platform to include custom algorithms.
### workflows
The workflow module contains example workflows.
To customize the workflow to your needs simply copy the provided functions and switch out the desired steps.
### io
Input and output functions that create GEMCAT models from different sources.
## utils
Contains common utility functions used throughout the package.
## verification
Functions to verify data integrity.
