# GEMCAT: Gene Expression-based Metabolite Centrality Analyses Tool

GEMCAT is a computational toolbox designed to predict metabolic alterations based on gene expression data. It's the 
accompanying software for our manuscript, "_GEMCAT — A new algorithm for gene expression-based prediction of metabolic alterations_."

## Quick links:
* **How to Cite:** [https://doi.org/10.1093/nargab/lqaf003](https://doi.org/10.1093/nargab/lqaf003)
* **PyPI:** [https://pypi.org/project/gemcat/](https://pypi.org/project/gemcat/)
* **Source Code (GitHub):** [https://github.com/MolecularBioinformatics/GEMCAT](https://github.com/MolecularBioinformatics/GEMCAT) 

## Important Considerations

* **Prediction Quality:** GEMCAT is still under refinement. It doesn't yet provide guidance on the
  statistical significance of predicted changes or any other measure of prediction quality. We
  recommend **filtering predictions for consistency** based on your domain knowledge.
* **Data Pre-filtering:** We **don't recommend pre-filtering transcriptomics and proteomics data based
  on significance**. This can negatively impact network coverage, as genes/proteins not present in the
  filtered dataset are implicitly considered "unchanged" by GEMCAT.
* **Graphical User Interface (GUI):** We are actively developing a user-friendly GUI for GEMCAT,
  which will be released soon. Stay tuned for updates on our GitHub repository and PyPI page! A **development version**
  of the GUI is currently hosted in a private repository; if you're interested in gaining early access,
  please contact **suraj.sharma@uib.no**.

---

## Compatibility
GEMCAT has been tested and is compatible with **Python >= 3.10** on Ubuntu and Windows operating systems.

## Installation
You can install GEMCAT in two ways:

1.  **Using pip (recommended):**
    ```bash
    pip install gemcat
    ```
2.  **From source (for developers or specific versions):**
    First, clone the repository, then install:
    ```bash
    git clone https://github.com/MolecularBioinformatics/GEMCAT.git
    cd gemcat
    pip install .
    ```
---

## How to Use GEMCAT

GEMCAT offers both a Python API for flexible, programmatic access and a command-line interface (CLI) for straightforward, scriptable use.

### Python Workflow with CobraPy

For more control and integration into existing Python projects, use the `workflow_standard` function:

```python
import gemcat as gc
import cobra # Assuming cobrapy is installed for model handling
import pandas as pd # For pd.Series

# Example usage (replace with your actual data and model)
# Make sure your mapped_genes_baseline and mapped_genes_comparison are pandas Series
# with gene/protein identifiers as the index.

# Example: Load a CobraPy model
# model = cobra.io.read_sbml_model("your_model.xml")

# Example: Create dummy mapped gene series
# mapped_genes_baseline = pd.Series([10, 20, 30], index=['geneA', 'geneB', 'geneC'])
# mapped_genes_comparison = pd.Series([15, 25, 35], index=['geneA', 'geneB', 'geneC'])

results = gc.workflows.workflow_standard(
    cobra_model=your_cobra_model, # Your loaded cobra.Model object
    mapped_genes_baseline=your_baseline_series, # pd.Series of baseline expression
    mapped_genes_comparison=your_comparison_series, # pd.Series of comparison expression
    adjacency=gc.adjacency_transformation.ATPureAdjacency, # Optional: Customize adjacency method
    ranking=gc.ranking.PagerankNX, # Optional: Customize ranking algorithm
    gene_fill=1.0 # Value to fill for genes not present in mapped_genes_comparison
)

print(results)
```
This function returns the changes in centrality relative to the baseline as a Pandas Series. If you're 
using fold-changes as your mapped_genes_comparison, you should provide a vector of all 1.0s for mapped_genes_baseline.

For further examples of using genome-scale metabolic networks from two different organisms refer: 
[An engineered human cell line with a functional deletion of the mitochondrial NAD transporter](https://github.com/MolecularBioinformatics/prm_manuscript/blob/main/jupyter_notebooks/pr_SLC25A51ko.ipynb),
[Patients with inflammatory bowel disease](https://github.com/MolecularBioinformatics/prm_manuscript/blob/main/jupyter_notebooks/pr_UC.ipynb),
[Training-induced metabolic changes in rats](https://github.com/MolecularBioinformatics/prm_manuscript/blob/main/jupyter_notebooks/pr_rats.ipynb),

### Command-Line Interface (CLI)

The CLI allows you to calculate differential centralities using gene expression data.

**Key Requirements for Input Files:**

* Your gene or protein identifiers **must be in the first column** of your expression file.
* These identifiers **must exactly match** those in your metabolic model. If you see a results list of all 1.0, it's
  a strong indicator of an identifier mismatch.
* Expression `.csv` files can be either comma- or tab-delimited.

**Common Workflows:**

1.  **Using a single file with pre-calculated fold-changes:**
    ```bash
    gemcat <model_file.xml> <expression_file.csv> -e <column_name> -o <result_file.csv>
    ```
    * `<model_file.xml>`: Path to your metabolic model file (SBML, JSON, or MAT format).
    * `<expression_file.csv>`: Path to your input file.
    * `<column_name>`: The name of the column in your CSV containing the fold-change values.
    * `<result_file.csv>`: The desired output file path.

2.  **Using two files (or one) with condition and baseline expression values:**
    ```bash
    gemcat <model_file.xml> <condition_file.csv> -e <condition_column_name> -b <baseline_file.csv> -c <baseline_column_name> -o <result_file.csv>
    ```
    * `<model_file.xml>`: Path to your metabolic model file (SBML, JSON, or MAT format).
    * `<condition_file.csv>`: Path to the file with expression values for your experimental condition.
    * `<baseline_file.csv>`: Path to the file with baseline expression values. If this is the same as the condition file, you can omit the `-b` flag and just use `<condition_file.csv>` as the second positional argument.
    * `<condition_column_name>`: Name of the column with condition expression data.
    * `<baseline_column_name>`: Name of the column with baseline expression data.

3.  **Using built-in models:**
    If you don't have a model file, GEMCAT can automatically access some common models by name:
    ```bash
    gemcat <model_name> <expression_file.csv> -e <column_name> -o <result_file.csv>
    ```
    Currently supported model names:
    * `recon3d`: [Recon3D](http://bigg.ucsd.edu/models/Recon3D)
    * `ratgem`: [Rat-GEM](https://github.com/SysBioChalmers/Rat-GEM)

**All CLI Parameters:**

* **Positional Arguments:**
    * `expression_file_path`: Path to your expression data file.
    * `model_file_path`: Path to your metabolic model file (or model name).
* **Optional Arguments:**
    * `-e --expressioncolumn`: Name of the column containing condition expression data (required for expression files).
    * `-b BASELINE, --baseline`: Path to the file containing baseline expression data.
    * `-c BASELINECOLUMN, --baselinecolumn`: Name of the column containing baseline expression data.
    * `-o OUTFILE, --outfile`: Path to write the output results.
    * `-v VERBOSE, --verbose`: Enables verbose output for detailed execution information.
    * `-l LOGFILE, --logfile`: Path to write logs.
 
---

## Modularity and Customization

GEMCAT is designed with a modular architecture, allowing you to easily swap out or append central components 
to customize its behavior. This is achieved by adhering to specifications laid out in the module base classes, particularly for:

* **Adjacency Transformation:** Defines how network adjacencies are calculated.
* **Expression Integration:** Handles mapping gene expression values onto reactions.
* **Ranking Components:** Implements different centrality ranking algorithms.

Any class inheriting from the abstract base classes in these modules can be exchanged.

---

## Core Modules Overview

* **`model`**: The central GEMCAT model structure, responsible for integrating workflows and calculating results.
* **`adjacency_transformation`**: Provides various approaches for calculating network adjacency and a platform for custom algorithms.
* **`expression`**: Manages the mapping of gene values onto reactions in the model via gene product rules, offering different algorithms along with a platform to create alternatives.
* **`ranking`**: Offers various ranking algorithms for the models along with a platform to include custom algorithms.
* **`workflows`**: Contains example workflows. To customize the workflow to your needs simply copy the provided functions and switch out the desired steps.
* **`cli`**: Command-line interface for GEMCAT.
* **`io`**: Input and output functions that create GEMCAT models from different sources.
* **`utils`**: Contains common utility functions used throughout the package.
* **`verification`**: Functions to verify data integrity.
* **`model_manager`**: Functionality for automatic downloading, storing, and retrieving of common models.

---

## Development

If you're contributing to GEMCAT:

* **Running Tests:**
    * Run all local tests with `pytest .`.
    * You can exclude slow-running tests by using `pytest . -m "not slow"`. These slow-running tests are
      integration tests with *real-world data* and will take 10-30 seconds each depending on your hardware.
* **Prerequisites:** Ensure you have [git lfs](https://git-lfs.com/) installed for tests that rely on large files.
* **Code Formatting:** Before committing, make sure your code is properly formatted using `isort` and `black`.
* **CI Pipeline:** The GitHub CI pipeline automatically checks for `isort`, `black`, and `pytest` compliance.

---

## Contact and Support

For questions, bug reports, or support, please open an issue on the 
[GitHub Issues page](https://github.com/MolecularBioinformatics/GEMCAT/issues). We will do our best to respond promptly.

For direct inquiries about the **development version of the GEMCAT GUI** or other specific questions, you can also contact:

* **Suraj Sharma:** suraj.sharma@uib.no
