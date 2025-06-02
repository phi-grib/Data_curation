"""
    Different functions to handle chemical miscelaneous issues such as plotting

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 14/06/2024, 04:28 PM
"""

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns

def plot_space(dataframe: pd.DataFrame, x: str, y: str, hue: str, style_var: str = None, figure_title: str = None, 
               figure_name: str = None, selected_comps: pd.DataFrame = None, selected_comps_id: str = None):
    """
        Function to plot the chemical space from a Dataframe after processing it with a PCA.
        It also includes the option to plot specific compounds within that dataframe if another dataframe is passes with the indicated compounds and 
        the PCA coordinates i.e. a subset of the main dataframe containint compounds of interest.

        :param dataframe:
        :param x:
        :param y:
        :param hue:
        :param style_var:
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
                    style=style_var,
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

def plot_space_interactive(dataframe: pd.DataFrame, x: str, y: str, hue_var: str, style_var: str,
               figure_title: str = None, figure_name: str = None,
               selected_comps: pd.DataFrame = None, selected_comps_id: str = None):
    """
    Function to plot the chemical space from a Dataframe after processing it with a PCA.
    This version uses Plotly for interactive visualization, allowing for zooming, panning,
    and hovering over data points. It visualizes two categorical variables: one using color
    (hue_var) and another using marker symbol (style_var).
    It also includes the option to highlight specific compounds within that dataframe.

    :param dataframe: The main DataFrame containing PCA coordinates and categorical variables.
    :param x: The name of the column in 'dataframe' to use for the x-axis (e.g., 'PC1').
    :param y: The name of the column in 'dataframe' to use for the y-axis (e.g., 'PC2').
    :param hue_var: The name of the column in 'dataframe' to use for coloring the points
                    (e.g., 'activity').
    :param style_var: The name of the column in 'dataframe' to use for changing marker symbols
                      (e.g., 'original_chemspace').
    :param figure_title: Optional title for the plot.
    :param figure_name: Optional name to save the figure as an HTML file (e.g., 'my_pca_plot').
    :param selected_comps: Optional DataFrame containing a subset of compounds to highlight.
                           This DataFrame must also contain 'x', 'y' coordinates and 'selected_comps_id'.
    :param selected_comps_id: The column name in 'selected_comps' that contains the identifiers
                              for the selected compounds (e.g., 'CAS_number').
    """

    # Create a copy to avoid modifying the original DataFrame directly
    df_plot = dataframe.copy()

    # Add a column to mark selected compounds for highlighting in Plotly
    # Initialize a column to False, then mark True for selected compounds
    df_plot['is_selected'] = False
    if isinstance(selected_comps, pd.DataFrame) and not selected_comps.empty:
        if not all(col in selected_comps.columns for col in [selected_comps_id, x, y]):
            print(f"Warning: 'selected_comps' DataFrame is missing required columns: {selected_comps_id}, {x}, or {y}. Skipping highlighting.")
        else:
            # Get the IDs of selected compounds
            selected_ids = selected_comps[selected_comps_id].tolist()
            # Mark rows in df_plot that match selected_ids
            df_plot['is_selected'] = df_plot[selected_comps_id].isin(selected_ids)

    # Create the scatter plot using plotly.express
    # color: for hue_var
    # symbol: for style_var
    # hover_data: to show additional information on hover
    fig = px.scatter(df_plot,
                     x=x,
                     y=y,
                     color=hue_var,
                     symbol=style_var, # Use symbol for the second categorical variable
                     hover_data=[selected_comps_id, hue_var, style_var], # Show these on hover
                     title=figure_title,
                     opacity=0.7,
                     height=700, # Set a fixed height for the plot
                     labels={x: f'{x} (PCA Component)', y: f'{y} (PCA Component)'} # More descriptive labels
                    )

    # Highlight selected compounds using a different trace or by adjusting properties
    # A more robust way to highlight is to add a separate trace for selected compounds
    if df_plot['is_selected'].any():
        selected_df_for_plot = df_plot[df_plot['is_selected']].copy()

        # Add a custom hover_name for selected compounds to show their ID prominently
        selected_df_for_plot['hover_name'] = selected_df_for_plot[selected_comps_id]

        fig.add_trace(
            go.Scatter(
                x=selected_df_for_plot[x],
                y=selected_df_for_plot[y],
                mode='markers+text', # Show markers and text
                marker=dict(
                    color='red',
                    size=15,
                    line=dict(width=2, color='DarkSlateGrey')
                ),
                text=selected_df_for_plot[selected_comps_id], # Text labels for selected compounds
                textposition="top center",
                name='Selected Compounds', # Name for the legend entry
                hoverinfo='text', # Only show text on hover for this trace
                hovertext=[f'<b>{row[selected_comps_id]}</b><br>PC1: {row[x]:.2f}<br>PC2: {row[y]:.2f}' for index, row in selected_df_for_plot.iterrows()]
            )
        )

    # Update layout for better appearance
    fig.update_layout(
        title_font_size=20,
        legend_title_text=f'{hue_var} / {style_var}', # Combined legend title
        hovermode="closest", # Improves hover experience
        margin=dict(l=40, r=40, t=80, b=40), # Adjust margins
        paper_bgcolor='LightSteelBlue', # Light background color for the plot area
        plot_bgcolor='white', # White background for the plotting region
        xaxis_title_font=dict(size=14),
        yaxis_title_font=dict(size=14)
    )

    # Show the plot
    fig.show()

    # Save the plot as an HTML file
    if figure_name:
        fig.write_html(f"{figure_name}.html")
        print(f"Interactive plot saved as {figure_name}.html")