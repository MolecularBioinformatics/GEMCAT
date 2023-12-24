# GEMCAT: Gene Expression based Metabolite Centrality Analyses Tool

## Compatibility
The package is tested for compatibility with Python >= 3.10 on Ubuntu and Windows.

## Installation
Install from pip:
```pip install gemcat```
NOTE: This will only work after publishing. 

Install from the locally cloned repository using pip: 
Create a standard installation for usage: ```pip install .```
Or an editable installation for development: ```pip install -e .``` (currently not working with pyproject.toml)


## Usage

### Standard workflow from the Command-Line Interface (CLI)
` gemcat ./expression_file.csv ./model_file.xml -e column_name -a pure -g 1.0 -o result_file.csv`

will run the standard workflow.

### Standard workflow from a CobraPy model
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
This will return the relative changes in centrality in the comparison relative to the baseline in a Pandas Series.

## Modularity and Configuration
GEMCAT is designed to be modular, and its central components can easily be swapped out or appended by other components adhering to the specifications laid out in the module base classes. (primarily adjacency transformation, expression integration, and ranking components)
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
