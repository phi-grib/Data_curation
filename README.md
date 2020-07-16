# SMILES_curation
Python code for handling SMILES of problematic molecules: organometallics, mixtures, peptides.

Installation:
* `git clone https://github.com/phi-grib/SMILES_curation.git`
* `conda activate <working_environment>`
* `cd SMILES_curation`
* `python setup.py install`

Now you can use the code for automatically curating and classifiying compounds from its structure.
Still under construction.

Future features to add:
- Mixture checking
- Sanitization error handler
- Optimizing the filtering function. Now it's a bit nested, which makes it more difficult to update and mantain for other developers

Usage example:
In the Jupyter Notebook that comes in the package there is an example of usage.
