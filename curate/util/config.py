"""
    Handles config.yaml and curation repository paths.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 22/02/2021, 13:00 PM
"""

import appdirs
import os
import sys

from pathlib import Path
from typing import Union, Tuple

from curate.util import utils

def ask_user() -> bool:
    """
        Utility function to obtain binary confirmation (yes/no) from the user
    """
    userinput = input()

    if userinput.lower() not in ['yes', 'no', 'y', 'n']:
        sys.stderr.write('Please write "yes", "no", "y" or "n"\n')
        return False
    elif userinput.lower() in ['yes', 'y']:
        return True
    return False

def configure(path: str = None, silent: bool = False) -> Union[Tuple[bool, str], Tuple[bool, dict]]:
    """
        Configures curation repository.
        Loads config.yaml and writes a correct curation repository path
        with the path provided by the user or a default from appdirs
        if the path is not provided.

        :param path: path to use as repo
        :param silent: boolean. If True, path is configured in the background.
    """
    
    success, config = utils.read_config()
    if not success:
        return False, config

    if silent:
        if path is not None:  
            source_dir = os.path.abspath(path)
        else:
            source_dir = os.path.dirname(os.path.dirname(__file__)) 

        curation_path = os.path.join(source_dir,'curation')

        try:
            if not os.path.isdir(curation_path):
                os.mkdir(curation_path)
        except Exception as e:
            return False, e
        
        utils.set_repositories(curation_path)
        
        sys.stderr.write('Curation repository set to {}\n'.format(curation_path))

        return True, config

    if path is None:  # set default

        # If datacur has been already configured, just show values in screen and return values
        if config['config_status'] == True:

            curation_path = config["curation_repository_path"]
            
            try:
                if not os.path.isdir(curation_path):
                    os.mkdir(curation_path)
            except Exception as e:
                return False, e

            sys.stderr.write("Current curation repository is {}\n".format(curation_path))

            return True, config

        # If datacur has not been already configured assign defaults
        curation_path = appdirs.user_data_dir('curation', 'datacur')
        
    else :

        try:
            source_dir = os.path.realpath(path)
        except:
            return False, "Input path {} is not recognized as a valid path".format(path)

        curation_path = os.path.join(source_dir,'curation')

    # at this point, paths must has been assigned
    sys.stderr.write("Curation repository will be set to {}\n".format(curation_path))
    sys.stderr.write("Continue?(y/n)\n")

    if ask_user():

        try:
            if not os.path.isdir(curation_path):
                os.mkdir(curation_path)
        except Exception as e:
            return False, e
        
        utils.set_repositories(curation_path)

        sys.stderr.write("Curation repository set to {}\n".format(curation_path))
        return True, config

    else:
        sys.stderr.write("Aborting...\n")
        return False, "Configuration aborted"


