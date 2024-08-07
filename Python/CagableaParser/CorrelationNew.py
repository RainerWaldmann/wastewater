import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def plotCorrelations(dfp: pd.DataFrame):
    nCols = len(dfp.columns)
    corrPlot_fig = make_subplots(rows=nCols, cols=nCols, start_cell="top-left",shared_xaxes=True,shared_yaxes=True,vertical_spacing=0.005, horizontal_spacing=0.005)
    for row_index, row_name in enumerate(dfp.columns):
        row_data = dfp[row_name]
        for column_index, column_name in enumerate(dfp.columns):
            column_data = dfp[column_name]

            trace = go.Scatter(y=column_data,
                               x=row_data,
                               mode="markers")
            #if column_index != row_index:
            if column_index <= row_index:
                corrPlot_fig.add_trace(trace,
                                         row=row_index+1, col=column_index+1)
                if column_index == 0:
                    corrPlot_fig.update_yaxes(title_text=row_name, row=row_index+1, col=column_index+1)
                if row_index == nCols - 1:
                    corrPlot_fig.update_xaxes(title_text=column_name, row=row_index+1, col=column_index+1)
                corrPlot_fig.update_traces(textposition='top center', row=row_index+1, col=column_index+1)

    corrPlot_fig.update_layout(width=400 * nCols,
                                 height=600 * nCols, title_text="Sequencing Depth (nReads per position)",
                                 yaxis=dict(
                                     showgrid=True
                                 ), title=""
                                 )
    #corrPlot_fig.update_xaxes(type="log")
    #corrPlot_fig.update_yaxes(type="log")
    corrPlot_fig.update_layout(showlegend=False, hovermode=False)
    corrPlot_fig.show()
    corr = dfp.corr(method='pearson')
    mask = np.triu(np.ones_like(corr, dtype=np.bool)) #https://medium.com/@szabo.bibor/how-to-create-a-seaborn-correlation-heatmap-in-python-834c0686b88e
    plt.figure(figsize=(16, 6))
    sns.heatmap(corr,annot=True, mask = mask)
    plt.show()