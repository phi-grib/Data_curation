"""
    Idata_child copied from Flame

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 26/03/2021, 10:33 PM
"""

from curate.idata import Idata


class IdataChild(Idata):

    def __init__(self, parameters, conveyor, input_source):

        Idata.__init__(self, parameters, conveyor, input_source)
