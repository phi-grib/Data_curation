# Data_curation

This Data curation tool provides a user-friendly CLI for treating raw data. It classifies the substances passed as input by filtering the SMILES and also applies a pre-processing of those structures to make them available for QSAR modelling.
It can be used in Jupyter Notebook as well, including data selection funcionalities that still need to be implemented in the CLI.

In the future, we plan to implement it as a part of Flame modelling software (https://github.com/phi-grib/flame).

## Installation:

This tool has been designed as a standalone application for Linux and also works for MacOS and Windows.
It needs a suitable conda working environment where to be installed. 

Download the repository:

```bash
git clone https://github.com/phi-grib/Data_curation.git
```

Go to the repository directory:

```bash
cd Data_curation/
```

Create a **new conda environment** and install the core scientific libraries using `conda`. This is the most crucial step as it will handle the complex binary dependencies correctly:
```bash
conda create -n datacuration python=3.10 numpy pandas rdkit -c conda-forge
```

Then activate the environment:

```bash
conda activate datacuration
```

Install Data curation:

```bash
pip install ".[dev]"
```

If you want to install it in editable mode, just add -e
```bash
pip install -e ".[dev]"
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

Data curation provides a command-line interface (CLI), dataset_curation.py, which allows for a fast and easy to implement curation of a file containing, at least, molecules with SMILES and an identifier (CAS, EC, name of the molecule, internal id etc...).

You can run the following commands from any terminal, in a computer where Data curation has been installed and the environment (datacuration) was activated (`source activate datacuration` in Linux).

First of all, we need to define an endpoint for our curated files:

```sh
datacur -c manage -a new -e myEndpoint
```

This creates a new entry in the curation repository. From now on, all our curation results will be stored there. In order to check the contents of the repository we can use the following command:

```sh
datacur -c manage -a list
```
Now the curation repository is totally configured and ready to store the outputs.
Let's curate a sample file:

```sh
datacur -i sample_file.xlsx -e myEndpoint -c curate -r -a chem
```

This will take the input file sample_file.xlsx and store the curated data as a pickle file in the curation repository (curated_data.pkl). With -r we asked the program to remove problematic structures and store them in a separate file for further revision. With -a chem we perform only a chemical curation (i.e. only the SMILES are treated). Since we haven't specified SMILES column nor ID column, the program uses a predifined name for each, being 'structure' for SMILES and 'name' for ID. If we want to specify those columns, which is recommended, we have to type:

```sh
datacur -i sample_file.xlsx -e myEndpoint -c curate -a chem -s smiles_colname -id id_colname -r
```

In that case, our input is an Excel file and the code handles this internally using Pandas option read_excel().
If we want to use another accepted format, like csv or tsv and we know we have a specific separator that is not a comma nor a tab, we can also specify the separator using the -sep option:

```sh
datacur -i sample_file.csv -e myEndpoint -sep ':' -c curate -a chem -s smiles_colname -id id_colname -r
```

If we have a large file containing lots of columns but we only want to keep some of them, then the --metadata or -m option is available. It will generate
an output only containing the most important columns for the curation plus the ones selected as metadata. Imagine that our file contains the columns meta1 and meta 2:
```sh
datacur -i sample_file.csv -e myEndpoint -c curate -a chem -s smiles_colname -id id_colname -r -m 'meta1,meta2'
```

Our output will be stored containig the columns id_colname, smiles_colname, structure_curated, substance_type_name, meta1 and meta2. If this option is not selected, all the columns are stored by default.

Also, there's an option to list all the output files in the endpoint directory using the following command:

```sh
datacur -c manage -e myEndpoint -a list
```

Finally, if we want to retrieve the curated data, we have the option download, where we can specify one of the accepted formats: tsv, csv, xlsx, sdf or json:

```sh
datacur -a download -c manage -e myEndpoint -f sdf
```

The output file will be stored in the local directory where the command has been executed.

## HTT-like files

This special curation option is specified with -a htt. It basically detects which columns include floats and puts them into a file called x_matrix which is returned as a tsv.
It performs the chemical curation only taking into account the molecule identifier and the SMILES column and it can be executed typing the following:

```sh
datacur -i sample_file_htt.xlsx -e myEndpoint -c curate -a htt -s smiles_colname -id id_colname -r
```

We can recover the output using the download option mentioned above.

## ChEMBL download

To use this option we have to pass to the -i command a valid ChEMBL ID (i.e CHEMBL230) or a tabular file (csv, tsv, xlsx) containing a field with ChEMBL IDs. Then select the -a option chembl. This will connect with the ChEMBL API and download the associated data to the ID in a pandas dataframe and the program will curate the structures present there.
If the input is a file, the dataframes will be concatenated and each register will be tagged with the corresponding ChEMBL ID.

Important note: this option works only for ChEMBL targets or proteins with enough assayed compounds. If you pass a compound or a non assayed target it will return a warning but it will continue working.

- Single ID:
```sh
datacur -i CHEMBL230 -a chembl -e myEndpoint -c curate -r
```

- Input file:
```sh
datacur -i input_file.csv -a chembl -e myEndpoint -c curate -r
```

We can choose to add -r or not to remove the problematic structures. To retrieve the results we just have to use the download option mentioned above.

## Data curation commands

| Command | Description |
| --- | --- |
| -i/ --infile |  Name of the input file used by the command. |
| -e/ --endpoint |  Name of the endpoint of our curation files. |
| -f/ --format | Output file formats that can be provided. Acceptable values are *xlsx*, *csv*, *tsv*, *sdf* and *json*. |
| -a/ --action | Management action to be carried out. Acceptable value are *silent*, *new*, *list*, *remove*, *export* and *download*. The meaning of these actions and examples of use are provided below. |
| -c/ --command | Specific action to be done. Acceptable values are *curate*, *split*, *config* and *manage*. |
| -d/ --directory | Defines the root directory for the curation repository. |
| -id/ --id_column | Column name containing the molecule identifier. |
| -s/ --smiles_col | Column name containing the SMILES string. |
| -m/ --metadata | Column names containing metadata of interest. |
| -sep/ --separator | If added, uses this argument as the input file separator. |
| -r/ --remove | If added, removes problematic structures after SMILES curation. |
| -h/ --help | Shows a help message on the screen |

Management commands deserve further description:


### Management commands

| Command | Example | Description |
| --- | --- | ---|
| silent | *datacur -c config -a silent* | Sets up the curation repository within the Data curation installation directory  |
| new | *datacur -c manage -a new -e MyEndpoint* | Creates a new entry in the curation repository named MyEndpoint  |
| remove | *datacur -c manage -a remove -e MyEndpoint* | Removes the specified endpoint from the curation repository |
| list | *datacur -c manage -a list* | Lists the endpoints present in the curation repository. If the name of an endpoint is provided, lists only the files within that endpoint directory  |
| export | *datacur -c manage -a export* | Exports the full curation directory as a tarball to the current working directory  |
| download | *datacur -c manage -a download -e myEndpoint -f myFormat* | Gets curated data in the specified format and stores it in the current working directory  |

## Technical details


### Using Data curation

Data curation was designed to be used in different ways, using diverse interfaces. For example:
- Using the `curate.py` command described above
- As a Python package in a Jupyter Notebook
- As a Python package, importing the classes and the functions.

## Licensing

Data curation was produced at the PharmacoInformatics lab (http://phi.upf.edu), in the framework of the eTRANSAFE project (http://etransafe.eu). eTRANSAFE has received support from IMI2 Joint Undertaking under Grant Agreement No. 777365. This Joint Undertaking receives support from the European Union’s Horizon 2020 research and innovation programme and the European Federation of Pharmaceutical Industries and Associations (EFPIA). 

Copyright 2021 Eric March (eric.march@upf.edu)

Data curation is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License as published by the Free Software Foundation version 3**.

Data curation is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Data curation. If not, see <http://www.gnu.org/licenses/>.
