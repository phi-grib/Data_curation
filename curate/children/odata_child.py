"""
    Odata_child copied from Flame

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 26/03/2021, 10:33 PM
"""

from curate.odata import Odata


class OdataChild(Odata):

    def __init__(self, parameters, conveyor):

        Odata.__init__(self, parameters, conveyor)
