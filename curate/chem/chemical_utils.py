"""
    Different functions to handle chemical miscelaneous issues such as plotting

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 14/06/2024, 04:28 PM
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def plot_space(dataframe: pd.DataFrame, x: str, y: str, hue: str, figure_title: str = None, 
               figure_name: str = None, selected_comps: pd.DataFrame = None, selected_comps_id: str = None):
    """
        Function to plot the chemical space from a Dataframe after processing it with a PCA.
        It also includes the option to plot specific compounds within that dataframe if another dataframe is passes with the indicated compounds and 
        the PCA coordinates i.e. a subset of the main dataframe containint compounds of interest.

        :param dataframe:
        :param x:
        :param y:
        :param hue:
        :param figure_title:
        :param figure_name:
        :param selected_comps:
        :param selected_comps_id:

    """

    plt.figure(figsize=(8,8))
    sns.scatterplot(data=dataframe,
                    x=x,
                    y=y,
                    hue=hue,
                    alpha=0.5,
                   legend=True)
    if figure_title:
        plt.title(figure_title)

    if isinstance(selected_comps, pd.DataFrame):
        for cas, x,y in zip(selected_comps[selected_comps_id],selected_comps[x], selected_comps[y]):
            plt.text(x+0.01, y, cas, horizontalalignment='left', size='medium', color='black')
            plt.scatter(x, y, color='red')

    if figure_name:
        plt.savefig("{}.png".format(figure_name), bbox_inches='tight')
        
    plt.show()