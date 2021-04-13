"""
    Copied from idata.py from Flame

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 26/03/2021, 11:05 PM
"""

import os
import sys
import pickle
import shutil
import tempfile
from joblib import Parallel, delayed
import pathlib

import numpy as np
from rdkit import Chem

from standardiser import standardise

from curate.util import utils

class Idata:

    def __init__(self, input_source: str):
        """
            Input data class to handle curation parameters

            Parameters
            ----------
            parameters: dict
                dict with model parameters

            conveyor: class
                main class to store the workflow results

            input_source: str
                SDF file with the molecules to use as training or predict        

        """

        # parameters and conveyor should have been initialized by the
        # parent class calling idata
        self.param = parameters
        self.conveyor = conveyor

        # self.format can inform if we are running in ghost mode
        # as part of an ensemble (low ensemble models)
        self.format = self.param.getVal('output_format')

        # path for temp files (fallback default)
        self.dest_path = '.'

        # add metainformation
        self.conveyor.addMeta('endpoint',self.param.getVal('endpoint'))
        self.conveyor.addMeta('version',self.param.getVal('version'))
        self.conveyor.addMeta('quantitative',self.param.getVal('quantitative'))
        
        input_type = self.param.getVal('input_type')
        self.conveyor.addMeta('input_type',input_type)


        # in case of top ensemble models...
        if input_type == 'model_ensemble':
            self.idata = input_source
            self.ifile = None
            randomName = 'flame-'+utils.id_generator()
            self.dest_path = os.path.join(tempfile.gettempdir(), randomName)

            #analyze first result to get the name of the input file
            ifile = 'ensemble input'
            try:
                ifile = input_source[0].getMeta('input_file')
            except:
                pass

            self.conveyor.addMeta('input_file',ifile)

        else:
            self.idata = None
            self.ifile = input_source
            self.dest_path = os.path.dirname(self.ifile)

            self.conveyor.addMeta('input_file',self.ifile)

