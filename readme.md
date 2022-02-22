# Pyreporter: a tool for reporter metabolite analysis

## Compatibility
Python 3.10 is not currently supported due to issues with NumPy.
The package has been tested for and will run on Python 3.6, 3.7, 3.8 and 3.9.
It has been tested on Windows 11 and the latest release of Ubuntu.

## Installation
Install using pip. Clone the repository, then simply navigate to this folder. 
Choose either of the following installation routes depending on your needs.

Create a standard installation for usage: ```pip install --user pyreporter ```
Create an editable installation for development: ```pip install --user -e pyreporter```

## Usage

### Standard workflow from a CobraPy model
```
import pyreporter as pr
results = workflow_Fang2012(
    cobra_model: cobra.Model,
    mapped_genes_baseline: pd.Series,
    mapped_genes_comparison: pd.Series,
    adjacency = pr.AdjacencyTransformation.ATPureAdjacency,
    ranking = pr.PageRank.PageRankNX,
    gene_fill = 0.,
    )
```
This will return the relative changes in centrality in the comparison relative to the baseline in a Pandas Series.

## Modularity and Configuration
Pyreporter is designed to be modular and its central components can easily be swapped out or appended by other components adhering to the specifications laid out in the module base classes. (primarily adjacency transformation, expression integration, and ranking components)
All classes inheriting from the abstract base classes laid out in the modules are swappable.

## Core modules
### Model
Core of the package is the PyReporter model structure that contains the model data, integrates the workflow and calculates the results.
### Adjacency
Different approaches can be used to calculate adjacency in the networks.
We offer alternatives and a platform to create custom algorithms for the model.
### Expression
Module covering the mapping of gene values onto reactions in the model via gene product rules.
Providing different algorithms along with a platform to create alternatives.
### PageRank
Module providing ranking algorithms for the models along with a platform to include custom algorithms.
### workflows
The workflow module contains example workflows.
To customize the workflow to your needs simply copy the provided functions and switch out the desired steps.
### io
Input and output functions that create PyReporter models from different sources.
## utils
Common utility functions used throughout the package.
## verification
Functions to verify data integrity.