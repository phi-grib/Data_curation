"""
    This code aims to handle the class imbalance by
    providing different functions.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 09/02/2021, 17:29 PM
"""

import imblearn
import numpy as np
import pandas as pd

from imblearn.combine import SMOTEENN, SMOTETomek
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler
from typing import Tuple

class ImbalanceData(object):

    """
        Class that includes the different imbalance correction algorithms, mainly
        from imblearn module and some handcrafted at phi-lab.
    """

    def __init__(self, imbalanced_data: pd.DataFrame, activity_field: str, imbalance_algorithm: str):
        """
            Initializes class. Needs the data to be processed and the algorithm/technique to be applied.

            :param imbalanced_data: dataframe with imbalanced classes
            :param activity field: column name that includes the activity
            :param imbalance_algorithm: algorithm to be used in the datasets. {oversampling, subsampling, SMOTEEN, SMOTETomek}

            TODO: SMOTETomek and ensemble modelling will be added as options in imbalance_algorithms.
        """

        self.imbalanced_data = imbalanced_data
        self.activity_field = activity_field
        self.imbalance_algorithm = imbalance_algorithm
        self.sampler = self.select_sampler_algorithm(imbalance_algorithm)

    def select_sampler_algorithm(self, imbalance_algorithm: str) -> imblearn.base.BaseSampler:
        """
            Checks the imbalance algorithm input and returns the proper sampler
            to be applied to the data

            :param imbalance_algorithm: name of the imbalance algorithm

            :return sampler: sampler to be used in the data
        """

        if imbalance_algorithm:
            if imbalance_algorithm.lower() == 'oversampling':
                sampler = RandomOverSampler(random_state=42)
            elif imbalance_algorithm.lower() == 'subsampling':
                sampler = RandomUnderSampler(random_state=42)
            elif imbalance_algorithm.lower() == 'smoteenn':
                sampler = SMOTEENN(random_state=42)
            elif imbalance_algorithm.lower() == 'smotetomek':
                sampler = SMOTETomek(random_state=42)
        
        return sampler
            
    ### Imbalance correction

    def imbalance_correction(self) -> pd.DataFrame:
        """
            Applies the selected sampler into the data and returns a resampled dataset
            
            Random Undersampling
            https://imbalanced-learn.org/stable/under_sampling.html

            Random Oversampling
            https://imbalanced-learn.org/stable/over_sampling.html

            SMOTEEEN/SMOTETomek
            https://imbalanced-learn.org/stable/combine.html
            
            :return resampled_set
        """
        
        if self.imbalance_algorithm.lower() in ('smoteenn','smotetomek'):
            x_, y_ = self.process_datasets(self.imbalanced_data)
        else:
            y_ = self.imbalanced_data[self.activity_field]
            x_ = self.imbalanced_data.drop(columns=self.activity_field)
            regular_cols = x_.columns
        
        x_resampled, y_resampled = self.sampler.fit_resample(x_, y_)
        
        if self.imbalance_algorithm.lower() in ('smoteenn','smotetomek'):
            x_flatten = x_resampled.flatten()
            resampled_set = self.get_proper_datasets(self.imbalanced_data, x_flatten)
        else:
            resampled_set = pd.DataFrame(data=x_resampled, columns=regular_cols)
            resampled_set.loc[:,self.activity_field] = y_resampled
        
        return resampled_set

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

    #### SMOTEEN/SMOTETomek part

    def process_datasets(self, dataset: pd.DataFrame) -> Tuple[np.ndarray, pd.DataFrame]:
        """
            This function is used to reshape the dataset in order to prepare it for SMOTEEN.

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