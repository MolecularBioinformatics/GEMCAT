# GEMCAT: Gene Expression-based Metabolite Centrality Analyses Tool
A computational toolbox associated with the manuscript entitled _GEMCAT — A new algorithm for gene expression-based prediction of metabolic alterations_. 
Cite using: https://www.biorxiv.org/content/10.1101/2024.01.15.575710v1

Note: We are still refining the tool. Particularly, GEMCAT does not yet provide guidance for significance of predicted changes or any other measure of prediction quality. We suggest filtering the predictions for consistency. We do not recommend pre-filtering the transcriptomics and proteomics data based on significance as this is affecting the network coverage which might negatively impact the prediction quality as genes/proteins not present in the dataset should be unchanged. 

## Compatibility
We tested the package for compatibility with Python >= 3.10 on Ubuntu and Windows.

## Installation
Install from pip:

```pip install gemcat```

Or clone the repository and install GEMCAT from there using:  

```pip install .```


## Usage

### Standard workflow from the Command-Line Interface (CLI)

Use a single file containing per-gene fold-changes to calculate the resulting differential centralities:
``` gemcat <./expression_file.csv> <./model_file.xml> -e <column_name> -o <result_file.csv>```
Make sure the .csv file is either comma- or tab-delimited.
`column_name` is the name of the column in the file containing the fold-change.

Alternatively, use two files (or one file) with expression values for condition and baseline:
``` gemcat <./condition_file.csv> <./model_file.xml> -e <condition_column_name> -b <./baseline_file> -c <baseline_column_name> -o <result_file.csv>```

If you do not have a model file ready, some models can be automatically accessed using their names:
``` gemcat ./expression_file.csv <model_name> -e column_name -o <result_file.csv>```

Model names currently supported are:
- ```recon3d```: [Recon3D](http://bigg.ucsd.edu/models/Recon3D)
- ```ratgem```: [Rat-GEM](https://github.com/SysBioChalmers/Rat-GEM)


Currently, GEMCAT supports models in SBML, JSON, and MAT formats.

Important points to remember:
Your gene or protein identifiers should be the first column of the expression file.
Make sure the gene or protein identifiers in your expression data file exactly match those in the model.
A results list of all 1.0 is a sure sign of no identifier matching.

Positional arguments:
- expression file path
- model file path

All parameters:
`-e --expressioncolumn` name of column containing condition expression data
`-b BASELINE, --baseline` file containing baseline expression data
`-c BASELINECOLUMN, --baselinecolumn` name of column containing baseline expression data
`-v VERBOSE, --verbose` enables verbose output
`-o OUTFILE, --outfile` write output to this file
`-l LOGFILE, --logfile` write logs to this file


### Standard workflow in Python using a CobraPy model
```
import gemcat as gc
results = gc.workflows.workflow_standard(
  cobra_model: cobra.Model,
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
GEMCAT is modular, and its central components can easily be swapped out or appended by other components 
adhering to the specifications laid out in the module base classes (primarily adjacency transformation, expression integration, and ranking components).
All classes inheriting from the abstract base classes laid out in the modules are exchangeable.

## Core modules
### Model
The core of the package is the GEMCAT model structure that contains the model data, integrates the workflow, and calculates the results.
### adjacency_transformation
Different approaches can be used to calculate adjacency in the networks.
We offer alternatives and a platform to create custom algorithms for the model.
### expression
Module covering the mapping of gene values onto reactions in the model via gene product rules.
Providing different algorithms along with a platform to create alternatives.
### ranking
Module providing ranking algorithms for the models along with a platform to include custom algorithms.
### workflows
The workflow module contains example workflows.
To customize the workflow to your needs simply copy the provided functions and switch out the desired steps.
### cli
Command-line interface for GEMCAT.
### io
Input and output functions that create GEMCAT models from different sources.
### utils
Contains common utility functions used throughout the package.
### verification
Functions to verify data integrity.
### model_manager
Functionality for automatic downloading, storing, and retrieving of common models.


## Development
You can run all local tests with `pytest .`. Default behavior is to also run integration tests, which takes time.
You can exclude slow running tests by using `pytest . -m "not slow"`.
These slow running tests are integration tests with _real world data_ and will take 10-30s each according to your hardware.

To run tests, make sure you have [git lfs](https://git-lfs.com/) installed and all the Tests are running.
Make sure to run `isort` and `black` to have properly formatted code.

The CI pipeline in GitHub will check with isort, black, and pytest.
