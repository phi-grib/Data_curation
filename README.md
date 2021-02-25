# Data_curation

This Data curation tool provides a user-friendly CLI for treating raw data. It classifies the substances passed as input by filtering the SMILES and also applies a pre-processing of those structures to make them available for QSAR modelling.
It can be used in Jupyter Notebook as well, including data selection funcionalities that still need to be implemented in the CLI.


## Installation:

This tool has been designed for Linux. It hasn't been tested neither for macOS configurations nor for Windows. 
It needs a suitable conda working environment where to be installed. 

Download the repository:

```bash
git clone https://github.com/phi-grib/Data_curation.git
```

Go to the repository directory:

```bash
cd Data_curation/
```

Create the **conda environment**:
```bash
conda env create -f environment.yml
```

Then activate the environment:

```bash
conda activate datacuration
```
And install the pip dependencies. Datacuration is included there:

``` bash
pip install -r requirements.txt
```

## Configuration

After installation is completed, you must run the configuration command to configure the directory where datacur will place the curated files. If it hasn't been configured previously the following command

```bash
datacur -c config
```
will suggest a default directory structure following the XDG specification in GNU/Linux.

To specify a custom path use the `-d` parameter to enter the root folder where the curated files will be placed:

```bash
datacur -c config -d /my/custom/path
```
will set up the curation repository to `/my/custom/path/curation`. 

Once Data curation has been configured, the current setting can be displayed using again the command 

```bash
datacur -c config
```

As a fallback, Data curation can also be configured using the following command

```bash
datacur -c config -a silent
```

This option sets up the curation repository within the Data curation installation directory (`curate\curation`). Unlike other options, this command does not ask permision to the end-user to create the directories or set up the repositories and is used internally by automatic installers and for software development. 

## Quickstarting

Data curation provides a command-line interface (CLI), curate.py, which allows for a fast and easy to implement curation of a file containing, an least, molecules with SMILES and an identifier (CAS, EC, name of the molecule, internal id etc...).

You can run the following commands from any terminal, in a computer where Data curation has been installed and the environment (datacuration) was activated (`source activate datacuration` in Linux).

Let's curate a sample file:

```sh
datacur -i sample_file.xlsx -o output_file -f xlsx -c curate -r
```

This will take the input file sample_file.xlsx and return output_file.xlsx since we specified the format using -f option. With -r we asked the program to remove problematic structures and store them in a separate file for further revision. Since we haven't specified SMILES column nor ID column, the program uses a predifined name for each, being 'structure' for SMILES and 'name' for ID. If we want to specify those columns, which is recommended, we have to type:

```sh
datacur -i sample_file.xlsx -o output_file -f xlsx -c curate -s smiles_colname -id id_colname -r
```

In that case, our input is an Excel file and the code handles this internally using Pandas option read_excel().
If we want to use another accepted format, like csv or tsv and we know we have a specific separator that is not a comma nor a tab, we can also specify the separator using the -sep option:

```sh
datacur -i sample_file.csv -o output_file -sep ':' -f xlsx -c curate -s smiles_colname -id id_colname -r
```

Also, if we want an sdf file to use it directly in a QSAR modelling tool, like Flame, we can put that in the -f option:

```sh
datacur -i sample_file.xlsx -o output_file -f sdf -c curate -s smiles_colname -id id_colname -r
```

The tool is still under construction, but thus far it works with the specified options from above.

## Data curation commands

| Command | Description |
| --- | --- |
| -i/ --infile |  Name of the input file used by the command. |
| -o/ --outfile | Name of the output file. |
| -f/ --format | Output file formats that can be provided. Acceptable values are *xlsx*, *csv*, *tsv* and *sdf*. |
| -a/ --action | Management action to be carried out. Acceptable value is *silent*. |
| -c/ --command | Specific action to be done. Acceptable values are *curate*, *split* and *config*. |
| -d/ --directory | Defines the root directory for the curation repository. |
| -id/ --id_column | Column name containing the molecule identifier. |
| -s/ --smiles_col | Column name containing the SMILES string. |
| -sep/ --separator | If added, uses this argument as the input file separator. |
| -r/ --remove | If added, removes problematic structures after SMILES curation. |
| -h/ --help | Shows a help message on the screen |

## Technical details


### Using Data curation

Data curation was designed to be used in different ways, using diverse interfaces. For example:
- Using the `curate.py` command described above
- As a Python package in a Jupyter Notebook
- As a Python package, importing the classes and the functions.

## Licensing

TODO
