"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 21/09/2020, 13:45 PM
"""

import imblearn
import numpy as np
import pandas as pd
import sys

from imblearn.combine import SMOTEENN, SMOTETomek
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
from sklearn.model_selection import train_test_split
from typing import Union, Optional, Tuple

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

    def split_main_dataset(self, imbalance_algorithm: str = None) -> pd.DataFrame:
        """
            This is the main function that returns the training set and the test set
            after applying the different proportions and the imbalance correction to the main set
            
            :param imbalance_algorithm: default is None. Otherwise, it can be Oversampling, Subsampling or SMOTEEN

            :return train_set, test_set:
        """

        train_set, test_set = self.get_sets(self.main_data, self.train_prop, self.test_prop)

        if imbalance_algorithm:
            if imbalance_algorithm.lower() == 'oversampling':
                train_set, test_set = self.random_oversampler_subsampler(train_set, test_set, 'oversampling')
            elif imbalance_algorithm.lower() == 'subsampling':
                train_set, test_set = self.random_oversampler_subsampler(train_set, test_set, 'subsampling')
            elif imbalance_algorithm.lower() == 'smoteen':
                train_set, test_set = self.smoteen_resample_sets(train_set, test_set)

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
    
    ### Imbalance correction

    def random_oversampler_subsampler_single_dataframe(self, df: pd.DataFrame, sampler: imblearn.base.BaseSampler) -> pd.DataFrame:
        """
            This function applies either Random Subsampling or Random Oversampling.
            Only for one dataframe used as training or test set.

            Random Undersampling
            https://imbalanced-learn.readthedocs.io/en/stable/under_sampling.html

            Random Oversampling
            https://imbalanced-learn.readthedocs.io/en/stable/over_sampling.html

            :param train_:
            :param test_:
            :param sampler: selected sampler. (Oversampling, Subsampling)

            :return resampled_train_set:
        """

        y_train = df[self.activity_field]
        x_train = df.drop(columns=self.activity_field)
        
        regular_cols = x_train.columns
        
        x_train_resampled, y_train_resampled = sampler.fit_resample(x_train, y_train)
        
        resampled_train_set = pd.DataFrame(data=x_train_resampled, columns=regular_cols)
        resampled_train_set.loc[:,self.activity_field] = y_train_resampled
        
        return resampled_train_set

    def random_oversampler_subsampler(self, train_: pd.DataFrame, test_: pd.DataFrame, algorithm: str) -> pd.DataFrame:
        """
            This function applies either Random Subsampling or Random Oversampling.

            Random Undersampling
            https://imbalanced-learn.readthedocs.io/en/stable/under_sampling.html

            Random Oversampling
            https://imbalanced-learn.readthedocs.io/en/stable/over_sampling.html

            :param train_:
            :param test_:
            :param algorithm: this can either be 'oversampling' or 'subsampling'

            :return resampled_train_set, resampled_test_set:
        """

        if algorithm == 'oversampling':
            sampler = RandomOverSampler(random_state=42)
        elif algorithm == 'subsampling':
            sampler = RandomUnderSampler(random_state=42)
        
        resampled_train_set = self.random_oversampler_subsampler_single_dataframe(train_, sampler)
        resampled_test_set = self.random_oversampler_subsampler_single_dataframe(test_, sampler)
        
        return resampled_train_set, resampled_test_set

    #### Double positive condition

    def double_positive_condition(self, df: pd.DataFrame) -> pd.DataFrame:
        """
            This function doubles the rows of the positive condition compounds.
            Is a simple way of oversampling

            :param df:

            :return positive_doubled_df:
        """

        pos_cond = df[df[self.activity_field] == 1]
        positive_doubled_df = pd.concat([df,pos_cond], ignore_index=True, axis=0)
        
        return positive_doubled_df

    #### Double negative condition

    def double_negative_condition(self, df: pd.DataFrame) -> pd.DataFrame:
        """
            This functions doubles the rows of the negative condition compounds.
            Simple way of oversampling.

            :param df:

            :return doubled_df:
        """

        neg_cond = df[df[self.activity_field] == 0]
        negative_doubled_df = pd.concat([df,neg_cond], ignore_index=True, axis=0)
        
        return negative_doubled_df

    #### Triple condition

    def triple_condition(self, df: pd.DataFrame, condition: str) -> pd.DataFrame:
        """
            This function triplicates the selected contidion (positive or negative).
            Simple way of oversampling.

            :param df:
            :param condtion: a string that can either be 'positive' or 'negative'

            :return triple_cond_df
        """

        if condition == 'positive':
            cond_df = df[df[self.activity_field] == 1]
        elif condition == 'negative':
            cond_df = df[df[self.activity_field] == 0]

        triple_cond_df = pd.concat([df, cond_df, cond_df], ignore_index=True, axis=0)

        return triple_cond_df

    #### SMOTEEN part

    def process_datasets(self, dataset: pd.DataFrame) -> Tuple[np.ndarray, pd.DataFrame]:
        """
            this function is used to reshape the dataset in order to prepare it for SMOTEEN.

            :param dataset:

            :return x_reshaped, y_:
        """

        y_ = dataset[self.activity_field]
        x_ = dataset.drop(columns=self.activity_field)
        x_reshaped = self.reshape_datasets(x_)
        
        return x_reshaped, y_

    def reshape_datasets(self, dataset: pd.DataFrame) -> np.ndarray:
        """
            This functions applies numpy's reshape to the dataset

            :param dataset:

            :return reshaped_dataset:
        """

        reshaped_dataset = np.reshape(dataset.index, (len(dataset.index), 1))
        
        return reshaped_dataset

    def get_proper_datasets(self, dataset: pd.DataFrame, indexes: np.ndarray) -> pd.DataFrame:
        """
            Removes duplicate registers after SMOTEEN algorithm is applied.
            Some of the registers that were being created had the same id in both datasets, thus
            giving errors.
            We've decided to remove those registers in order to keep datasets cleaned of noise.

            :param dataset:

            :return proper_dataset:
        """

        proper_dataset = pd.DataFrame(columns=dataset.columns,index=indexes).reset_index()
        index_to_drop = set()
        
        for i, row in proper_dataset.iterrows():
            idx = row['index']
            if idx in dataset.index:
                proper_dataset.loc[i] = dataset.loc[idx]
            else:
                index_to_drop.add(idx)
        
        proper_dataset.drop(proper_dataset.loc[proper_dataset['index'].isin(index_to_drop)].index, inplace=True)
        
        return proper_dataset
                
    def smoteen_resample_sets(self, train_: pd.DataFrame, test_: pd.DataFrame) -> pd.DataFrame:
        """
            This function applies SMOTEEN once train and test set are split.

            :param train_:
            :param test_:

            :return resampled_train_set, resampled_test_set:
        """

        smote_enn = SMOTEENN(random_state=42)
        
        x_train_index_reshape, y_train = self.process_datasets(train_)
        x_test_index_reshape, y_test = self.process_datasets(test_)

        x_train_resampled, y_train_resampled = smote_enn.fit_resample(x_train_index_reshape, y_train)
        x_test_resampled, y_test_resampled = smote_enn.fit_resample(x_test_index_reshape, y_test)
        
        x_train_res_flatten = x_train_resampled.flatten()
        x_test_res_flatten = x_test_resampled.flatten()
        
        resampled_train_set = self.get_proper_datasets(train_, x_train_res_flatten)
        resampled_test_set = self.get_proper_datasets(test_, x_test_res_flatten)
        
        return resampled_train_set, resampled_test_set

    def resample_main_set(self, main_set:pd.DataFrame) -> pd.DataFrame:
        """
            This function applies SMOTEEN to the whole dataset before splitting into train and test sets.
            It is deprecated since we found that is better to apply it after splitting, but we keep the function.

            :param main_set:

            :return resampled_set:
        """

        smote_enn = SMOTEENN(random_state=42)
        
        x, y = self.process_datasets(main_set)
        
        x_resampled, y_resampled = smote_enn.fit_resample(x,y)
        
        x_flatten = x_resampled.flatten()
        
        resampled_set = self.get_proper_datasets(main_set, x_flatten)
        
        return resampled_set