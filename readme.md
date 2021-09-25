# pyreporter tools for reporter metabolite analysis

## Dependencies
Dependencies can all be installed from pip: 
pip install --user cobra pandas numpy

## Installation
Installation proceeds via pip. Simply navigate to this folder, then call either of the following:

pip install --user pyreporter (to create a standard install)
pip install --user -e pyreporter (to create an editable install)

## Usage 
import pyreporter
pyreporter.calc_reporters_from_expression(model, data, p_col, f_col)

where model is your cobra.Model object, data is your pandas.DataFrame with columns p_col containing p-values and f_col containing log2-fold-changes

The pagerank and local_changes subpackages have the functions separated out in case you need to adapt the code. 
Just follow the public functions (functions not starting in _).