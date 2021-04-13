"""
    Learn_child copied from Flame

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 13/04/2021, 09:53 AM
"""

from flame.learn import Learn


class LearnChild(Learn):

    def __init__(self, parameters, conveyor):

        Learn.__init__(self, parameters, conveyor)
