"""
    This code aims to handle the class imbalance by
    providing different functions.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 09/02/2021, 17:29 PM
"""

import imblearn
import numpy as np
import pandas as pd
import sys

from imblearn.combine import SMOTEENN, SMOTETomek
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import RandomUnderSampler, NearMiss, ClusterCentroids, EditedNearestNeighbours, RepeatedEditedNearestNeighbours, AllKNN, InstanceHardnessThreshold
from typing import Tuple, Optional, Union

from curate.util import get_logger

LOG = get_logger(__name__)
class ImbalanceData(object):

    """
        Class that includes the different imbalance correction algorithms, mainly
        from imblearn module and some written at phi-lab.
    """

    def __init__(self, imbalanced_data: pd.DataFrame, molecule_id: str, activity_field: str, imbalance_algorithm: str, descriptors: Optional[str] = None,
                 sampling_strategy: Union[float, str, dict, callable] = 'auto'):
        """
            Initializes class. Needs the data to be processed and the algorithm/technique to be applied.

            :param imbalanced_data: dataframe with imbalanced classes
            :param molecule_id: column name for molecule identifier, such as CAS, ChEMBL ID, InChi
            :param activity field: column name that includes the activity
            :param imbalance_algorithm: algorithm to be used in the datasets. {oversampling, subsampling, SMOTEEN, SMOTETomek, NearMiss1,
            NearMiss2, NearMiss3, cluster_centroid, edited_knn, rep_edited_knn, all_knn, iht}
            :param descriptors: it can either be empty of filled with a string. Currently only works with rdkit.

            TODO: SMOTETomek and ensemble modelling will be added as options in imbalance_algorithms.
        """

        self.imbalanced_data = imbalanced_data
        self.molecule_id = molecule_id
        self.activity_field = activity_field
        self.imbalance_algorithm = imbalance_algorithm
        self.sampling_strategy = sampling_strategy
        self.sampler = self.select_sampler_algorithm(imbalance_algorithm)

        if descriptors:
            if descriptors.lower() == 'rdkit':
                from rdkit.Chem import Descriptors
                self.descriptors = [x[0] for x in Descriptors._descList]
                self.descriptors.remove('Ipc')

    def select_sampler_algorithm(self, imbalance_algorithm: str) -> imblearn.base.BaseSampler:
        """
            Checks the imbalance algorithm input and returns the proper sampler
            to be applied to the data

            :param imbalance_algorithm: name of the imbalance algorithm

            :return sampler: sampler to be used in the data
        """

        if imbalance_algorithm:
            if imbalance_algorithm.lower() == 'oversampling':
                sampler = RandomOverSampler(random_state=42, sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'subsampling':
                sampler = RandomUnderSampler(random_state=42, sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'smoteenn':
                sampler = SMOTEENN(random_state=42, sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'smotetomek':
                sampler = SMOTETomek(random_state=42, sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'NearMiss1':
                sampler = NearMiss(version=1, sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'NearMiss2':
                sampler = NearMiss(version=2, sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'NearMiss3':
                sampler = NearMiss(version=3, sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'cluster_centroid':
                """ TODO: check how to recover indexes from cluster centroid sampler"""
                LOG.error('Cluster centroid still not available. Please consider another sampling strategy\n')
                # sampler = ClusterCentroids(random_state=42, sampling_strategy=self.sampling_strategy)
                sys.exit(1)
            elif imbalance_algorithm.lower() == 'edited_knn':
                sampler = EditedNearestNeighbours(sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'rep_edited_knn':
                sampler = RepeatedEditedNearestNeighbours(sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'all_knn':
                sampler = AllKNN(sampling_strategy=self.sampling_strategy)
            elif imbalance_algorithm.lower() == 'iht':
                from sklearn.linear_model import LogisticRegression
                sampler = InstanceHardnessThreshold(random_state=42,
                                estimator=LogisticRegression(
                                solver='lbfgs', multi_class='auto'), sampling_strategy=self.sampling_strategy)

        return sampler
            
    ### Imbalance correction

    def imbalance_correction(self) -> pd.DataFrame:
        """
            Applies the selected sampler into the data and returns a resampled dataset
            
            Random Undersampling, NearMiss, Cluster Centroids, Edited KNN, Repited Edited KNN, All KNN, IHT
            https://imbalanced-learn.org/stable/under_sampling.html

            Random Oversampling
            https://imbalanced-learn.org/stable/over_sampling.html

            SMOTEEEN/SMOTETomek
            https://imbalanced-learn.org/stable/combine.html
            
            :return resampled_set
        """
        
        if self.imbalance_algorithm.lower() in ('smoteenn','smotetomek'):
            x_, y_ = self.process_datasets(self.imbalanced_data)
        elif self.imbalance_algorithm.lower() in ('nearmiss1','nearMmss2', 'nearmiss3', 
                                                  'cluster_centroid', 'edited_knn', 'rep_edited_knn', 'all_knn', 'iht'):
            x_, y_ = self.get_numerical_descriptors_only(self.imbalanced_data)
        else:
            y_ = self.imbalanced_data[self.activity_field]
            x_ = self.imbalanced_data.drop(columns=self.activity_field)
            regular_cols = x_.columns

        x_resampled, y_resampled = self.sampler.fit_resample(x_, y_)
        
        if self.imbalance_algorithm.lower() in ('smoteenn','smotetomek'):
            x_flatten = x_resampled.flatten()
            resampled_set = self.get_proper_datasets(self.imbalanced_data, x_flatten)
        elif self.imbalance_algorithm.lower() in ('nearmiss1','nearMmss2', 'nearmiss3', 
                                                  'cluster_centroid', 'edited_knn', 'rep_edited_knn', 'all_knn', 'iht'):
            resampled_set = self.recover_former_indexes_from_sampler(self.imbalanced_data, self.sampler)
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
    

    #### Undersampling algorithms dataset preparation

    def get_numerical_descriptors_only(self, dataset:pd.DataFrame) -> pd.DataFrame:
        """
            Splits the dataset into X and y, both containing only numerical data

            :return X:
            :return y:
        """
        
        X = dataset[self.descriptors]
        y = dataset[self.activity_field]

        return X, y

    def recover_former_indexes_from_sampler(self, dataset: pd.DataFrame, sampler: imblearn.base.BaseSampler) -> pd.DataFrame:
        """
            Uses the built-in function sample_indices_ to recover the indexes from the sampler
            and rebuilt the new dataframe once resampled

            :return resampled_dataset:
        """

        resampled_dataset = dataset.loc[dataset.index.isin(sampler.sample_indices_)]

        return resampled_dataset