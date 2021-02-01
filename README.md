# Data_curation
Python code for data preparation for modelling:
- handling SMILES of problematic molecules: organometallics, mixtures, peptides.
- data sampling
- data selection

Installation:
* `git clone https://github.com/phi-grib/Data_curation.git`
* `conda activate <working_environment>`
* `cd Data_curation`
* `python setup.py install`

Now you can use the code for automatically curating and classifiying compounds from its structure.
Still under construction.

Future features to add:
- Mixture checking
- Optimizing the filtering function. Now it's a bit nested, which makes it more difficult to update and mantain for other developers

Usage example:<br>
In the Jupyter Notebook that comes in the package there is an example of usage.
