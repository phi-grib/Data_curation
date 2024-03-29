"""
    Managing functions for Data curation CLI.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 18/02/2021, 17:49 PM
"""

import os
import pandas as pd
import pathlib
import pickle
import shutil
import tarfile

from typing import Tuple, Union

from curate.parameters import Parameters
from curate.util import utils, get_logger

LOG = get_logger(__name__)

def set_curation_repository(path: str = None):
    """
        Set the path to the curation repository

        :param path: string indicating the path of the curation repository
    """

    utils.set_curation_repository(path)

    LOG.info('Model repository updated to {}\n'.format(path))

    return True, 'curation repository updated\n'

#### API functions

def action_new(curation_path: str) -> Tuple[bool, str]:
    """
        Create a new curation endpoint tree, using the given name.
        
        :param curation_path: curation endpoint in curation repository where output will be saved

        :return bool: True when evertyhing has workded, otherwise False.
        :return str: strings that would be the equivalent to the standard error.
    """

    if not curation_path:
        return False, 'empty endpoint curation label\n'

    # importlib does not allow using 'test' and issues a misterious error when we
    # try to use this name. This is a simple workaround to prevent creating paths 
    # with this name 
    if curation_path == 'test':
        return False, 'the name "test" is disallowed, please use any other name\n'

    # curation endpoint directory
    ndir = pathlib.Path(utils.curation_tree_path(curation_path))

    # check if there is already a tree for this endpoint
    if ndir.exists():
        return False, "Endpoint {} already exists\n".format(curation_path)

    try:
        ndir.mkdir(parents=True)
        LOG.info("{} created\n".format(ndir))
    except:
        return False, "Unable to create path for {} endpoint\n".format(curation_path)

    # Copy classes skeletons to ndir
    wkd = pathlib.Path(os.path.dirname(os.path.abspath(__file__)))
    
    # copy parameter yml file
    params_path = wkd / 'children' / 'curation_parameters.yaml'
    shutil.copy(params_path, ndir)

    LOG.info("New endpoint {} created\n".format(curation_path))
    
    return True, "new endpoint {} created".format(curation_path)

def get_creation_date(endpoint_path: str) -> str:
    """
        Returns the creation date of the selected endpoint dir in the curation repository.

        :param endpoint_path: complete path to endpoint directory.

        :return creation_date: string with the creation date of the endpoint dir.
    """

    import datetime

    creation_date = datetime.datetime.fromtimestamp(os.stat(endpoint_path).st_mtime).strftime("%d-%b-%Y")
    
    return creation_date

def action_list(curation_dir: str) -> Tuple[bool, str]:
    """
        In no argument is provided lists all endpoints present at the repository 
        otherwyse lists all files for the endpoint provided as argument.

        :param curation_dir: path to the endpoint in curation repo
    """

    # if no name is provided, just list the different curation dirs
    if not curation_dir:
        rdir = utils.curation_repository_path()
        if os.path.isdir(rdir) is False:
            return False, 'the curation repository path does not exist. Please run "datacur -c config".\n'

        num_curs = 0
        LOG.info('Curation endpoints found in repository:\n')
        for x in os.listdir(rdir):
            xpath = os.path.join(rdir,x)
            # discard if the item is not a directory
            if not os.path.isdir(xpath):
                continue
            num_curs += 1
            creation_date = get_creation_date(xpath)
            LOG.info("\n{} {}\n".format(x, creation_date))
            
        LOG.info("\nRetrieved list of curation endpoints from {}\n".format(rdir))

        return True, "{} endpoints found".format(num_curs)

    else:
        # if a path name is provided, list files
        base_path = utils.curation_tree_path(curation_dir)
        num_files = 0
        LOG.info('Files found in curation endpoint {}:\n'.format(curation_dir))
        for x in os.listdir(base_path):
            num_files += 1
            xpath = os.path.join(base_path,x)
            creation_date = get_creation_date(xpath)
            LOG.info("\n{} {}\n".format(x, creation_date))

        return True, "Endpoint {} has {} files".format(curation_dir, num_files)

def action_remove(curation_endpoint: str) -> Tuple[bool, str]:
    """
        Remove the curation endpoint directory indicated as 
        argument

        :param curation_endpoint: curation endpoint to be removed
    """

    if not curation_endpoint:
        return False, 'Empty curation endpoint'

    rdir = utils.curation_tree_path(curation_endpoint)
    if not os.path.isdir(rdir):
        return False, '{} not found'.format(curation_endpoint)

    shutil.rmtree(rdir, ignore_errors=True)
    LOG.info("Curation endpoint dir {} has been removed\n".format(curation_endpoint))

    return True, "Curation endpoint dir {} has been removed".format(curation_endpoint)

def action_dir() -> Tuple[bool,Union[str,list]]:
    """
        Returns a list of curation endpoints and files

        :return bool:
        :return str:
        :return results:
    """

    # get curation repo path

    cur_path = pathlib.Path(utils.curation_repository_path())
    if cur_path.is_dir() is False:
        return False,  'Curation repository path does not exist. Please run "datacur -c config".\n'

    # get directories in curation repo path

    dirs = [x for x in cur_path.iterdir() if x.is_dir()]
    results = []

    for directory in dirs:
        dir_dict = {}
        # I convert directory, which is a PosixPath object, into a string
        # this will not work in Windows (MP)
        # directory_string = str(directory).split('/')[-1]
        directory_string = directory.parts[-1]

        # Not showing statistics files in the list of files within the directory
        dir_dict['curation_endpoint'] = directory_string
        dir_dict['creation_date'] = get_creation_date(directory)
        dir_dict['curation_output'] = 'unk'

        for file_ in os.listdir(directory):
            if file_.startswith('curated_data') and 'head' not in file_:
                dir_dict['curation_output'] = file_
        
        results.append(dir_dict)
    
    return True, results

def action_kill(curation_endpoint: str) -> Tuple[bool,str]:
    """
        Removes the endpoint tree described by the argument.

        :param curation_endpoint: path to curation endpoint in the repo.

        :return bool:
        :return str:
    """

    if not curation_endpoint:
        return False, 'Empty endpoint name'

    ndir = utils.curation_tree_path(curation_endpoint)

    if not os.path.isdir(ndir):
        return False, "Model {} not found".format(curation_endpoint)

    try:
        shutil.rmtree(ndir, ignore_errors=True)
    except:
        return False, "Failed to remove model {}".format(curation_endpoint)

    LOG.info("Model {} removed\n".format(curation_endpoint))
    
    return True, "Model {} removed".format(curation_endpoint)

def action_export(curation_endpoint: str) -> Tuple[bool,str]:
    """
        Exports the whole curation endpoint tree indicated in the argument as a single
        tarball file with the same name.

        :param curation_endpoint: path to curation endpoint in the repo.
    """

    if not curation_endpoint:
        return False,  'Empty endpoint name'

    current_path = os.getcwd()
    exportfile = os.path.join(current_path,curation_endpoint+'.tgz')

    base_path = utils.curation_tree_path(curation_endpoint)

    if not os.path.isdir(base_path):
        return False, 'Unable to export, endpoint directory not found'

    # change to curation repository to tar the file from there
    os.chdir(base_path)

    itemend = os.listdir()
    itemend.sort()

    with tarfile.open(exportfile, 'w:gz') as tar:
        for iversion in itemend:
            if not os.path.isdir(iversion):
                continue
            tar.add(iversion)

    # return to current directory
    os.chdir(current_path)
    LOG.info("Endpoint {} exported as {}.tgz\n".format(curation_endpoint,curation_endpoint))

    return True, "Endpoint {} exported as {}.tgz".format(curation_endpoint,curation_endpoint)

def action_info_curation(endpoint: str) -> Tuple[bool, Union[str,dict]]:
    """
        Returns a list of curation statistics

        :param endpoint: curation endpoint to collect statistics from
        
        :return bool:
        :return str:
        :return stats_dict:
    """

    # get curation endpoint path

    endpoint_curation = pathlib.Path(utils.curation_tree_path(endpoint))
    if endpoint_curation.is_dir() is False:
        return False,  'Curation endpoint path does not exist.\n'

    # get statisics files in curation endpint
    stats_file = os.path.join(endpoint_curation,'statistics.pkl')
    
    if not os.path.isfile(stats_file):
        return False,  'Statistics file does not exist.\n'

    stats = []
    with (open(stats_file, "rb")) as openfile:
        stats = pickle.load(openfile)

    return True, stats

def action_header_curation(endpoint: str) -> Tuple[bool, Union[str,dict]]:
    """
        Returns a dataframe of curation output header

        :param endpoint: curation endpoint
        
        :return bool:
        :return str:
        :return head_:
    """

    # get curation endpoint path

    endpoint_curation = pathlib.Path(utils.curation_tree_path(endpoint))
    if endpoint_curation.is_dir() is False:
        return False,  'Curation endpoint path does not exist.\n'

    # get header file in curation endpint
    header_file = os.path.join(endpoint_curation,'curated_data_head.pkl')
    
    if not os.path.isfile(header_file):
        return False,  'Curation header file does not exist.\n'

    head_ = []
    with (open(header_file, "rb")) as openfile:
        head_ = pickle.load(openfile)
    
    return True, head_

def action_curation_results(args: list) -> Tuple[bool, Union[dict,str]]:
    """
        Returns the output file in the specified format and problematic structures file if the option was selected

        :param args:
        
        :return bool:
        :return str:
        :return dict:
    """

    # get current working directory to store curated data
    current_path = os.getcwd()
    
    # get curation endpoint path
    endpoint_curation = utils.curation_tree_path(args.endpoint)
    if not os.path.exists(os.path.dirname(endpoint_curation)):
        return False,  'Curation endpoint path does not exist.\n'
    
    # get curation parameters
    params = Parameters()
    success, curation_parameters = params.get_parameters(endpoint_curation)
    
    if not success:
        return success, curation_parameters

    identifier = curation_parameters['molecule_identifier']
    smiles_column = curation_parameters['structure_column']
    log = 'curations not found for {} directory'.format(args.endpoint)

    output_handling(endpoint=endpoint_curation,
                    filename='curated_data',
                    log=log,
                    currpath=endpoint_curation,
                    outfile_type=args.format,
                    smiles=smiles_column,
                    id=identifier)

    # check if remove problematic is true or false and if curation type is htt.
    # if flag = 1, it downloads problematic structures file and/or x matrix in a tarball.

    flag = 0
    list_of_files = []

    if curation_parameters['curation_type'] == 'htt':
        flag = 1
        log = 'X matrix pickle does not exist. Please use -a htt option when curating a dataset.\n'

        output_handling(endpoint=endpoint_curation,
                        filename='x_matrix',
                        log=log,
                        currpath=endpoint_curation,
                        outfile_type='tsv',
                        smiles=None,
                        id=None)
        
        list_of_files.append('x_matrix.tsv')

    if curation_parameters['remove_problematic'] == 'true':
        flag = 1
        log = 'Problematic structures pickle does not exist. Please use -r option when curating a dataset.\n'

        output_handling(endpoint=endpoint_curation,
                        filename='problematic_structures_removed',
                        log=log,
                        currpath=endpoint_curation,
                        outfile_type='xlsx',
                        smiles=smiles_column,
                        id=identifier)

        list_of_files.append('problematic_structures_removed.xlsx')
    
    if flag == 1:
        exportfile = os.path.join(current_path,'curation.tgz')
        list_of_files.append('.'.join(['curated_data',args.format]))
        
        os.chdir(endpoint_curation)

        with tarfile.open(exportfile, 'w:gz') as tar:
            for file_ in list_of_files:
                if not os.path.isfile(file_):
                    continue
                tar.add(file_)
                os.remove(file_)
        
        outfile_name = 'curation.tgz'
        os.chdir(current_path)
    else:
        outfile_name = '.'.join(['curated_data',args.format])
        source_curation = os.path.join(endpoint_curation,outfile_name)
        output_curation = os.path.join(current_path,outfile_name)
        shutil.copy(source_curation, output_curation)
        os.remove(source_curation)

    return True, "Curated data downloaded successfully as {}".format(outfile_name)

def output_handling(endpoint: str, filename: str, log: str, currpath: str, outfile_type: str, smiles: str, id: str):
    """
        Writes the proper output file given the parameters.

        :param endpoint: endpoint directory path
        :param filename: name of the output file wihtout the format
        :param log: log to show if the file is not found
        :param currpath: current working directory path
        :param outfile_type: format of the output file
        :param smiles: smiles column name
        :param id: modelcule identifier column name
    """

    pickle = os.path.join(endpoint, '.'.join([filename,'pkl']))
    
    if not os.path.isfile(pickle):
        return False, log
    
    data = pd.read_pickle(pickle)
    output_file = os.path.join(currpath,filename)
    
    utils.format_output(data = data, 
                    outfile_type = outfile_type, 
                    outfile_path = output_file, 
                    smiles_column = smiles, 
                    identifier = id)

def action_parameters(curation_path: str, oformat: str = 'text') -> Union[Tuple[bool, str],Tuple[bool, object]]:
    """
        Returns an object with the curation parameters for a given endpoint

        :param curation_path:
        :param oformat:
    """

    if curation_path is None:
        return False, 'Empty curation label'

    from curate.parameters import Parameters

    param = Parameters()
    success, results = param.loadYaml_curation(curation_path)

    if not success:
        LOG.error("Error obtaining parametes for curation endpoint {} : {}\n".format(curation_path, results))
        return False, results

    if oformat != 'text':
        return True, param

    else:
        yaml = param.dumpYAML_curation()
        for line in yaml:
            print(line)

        return True, 'Parameters listed', yaml