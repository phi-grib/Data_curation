"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 21/09/2020, 13:45 PM
"""

import pandas as pd
import sys

from sklearn.model_selection import train_test_split

class Selection(object):
    
    """
        This object aims to authomatise the dataset selection from CII database to build the models.
        It also includes the imbalance correction, which is applied if the user needs to.
    """

    def __init__(self, dataframe: pd.DataFrame, train_prop: float, test_prop: float, activity_field: str):
        """
            Initializes class

            :param dataframe:
            :param train_prop:
            :param test_prop:
            :param activity_field: column name with the activity data
        """

        train_test_proportion = train_prop + test_prop
        if train_test_proportion != 1.0:
            sys.stderr.write('Please introduce a valid proportion of train and test set. The sum of both should be equal to 1.0.\n')
            sys.exit()
        else:
            self.train_prop = train_prop
            self.test_prop = test_prop

        self.main_data = dataframe  
        self.activity_field = activity_field

    ### Selection main function

    def split_main_dataset(self) -> pd.DataFrame:
        """
            This is the main function that returns the training set and the test set
            after applying the different proportions.
            

            :return train_set:
            :return test_set:
        """

        train_set, test_set = self.get_sets(self.main_data, self.train_prop, self.test_prop)

        return train_set, test_set

    ### Set selection

    def get_sets(self, df: pd.DataFrame, train_prop: float, test_prop: float) -> pd.DataFrame:
        """
            This function performs the train and test set selection as explained in the link below:
            https://stats.stackexchange.com/questions/394056/splitting-into-train-and-test-sets-keeping-class-proportions

            :param df:
            :param train_prop:
            :param test_prop:

            :return train_set, test_set:
        """

        y = df[self.activity_field]
        x = df.drop(columns=self.activity_field)
        
        X_train, X_test, y_train, y_test = train_test_split(x, y, train_size=train_prop, test_size=test_prop, stratify = y, random_state=42)
        
        train_set = pd.concat([X_train, y_train], axis=1).reindex(X_train.index)
        test_set = pd.concat([X_test, y_test], axis=1).reindex(X_test.index)
        
        return train_set, test_set