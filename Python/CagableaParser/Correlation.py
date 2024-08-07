import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

def correlations(dataframe) :
    pd.plotting.scatter_matrix(dataframe,hist_kwds={'bins':50})
    scatterMatrix = px.scatter_matrix(dataframe, width=2500, height = 2500)


    corr = dataframe.corr(method='pearson')
    print(corr)

    #corr.style.background_gradient(cmap='coolwarm')



    layout = go.Layout(
        title_text="title",
        title_x=0.5,
        width=600,
        height=600,
        xaxis_showgrid=False,
        yaxis_showgrid=False,
        yaxis_autorange='reversed'
    )

    #fig=go.Figure(data=[dfpForHeatMaps],layout = layout)
    #fig.show()


    fig = px.imshow(corr,
                    x=corr.columns.values,
                    y=corr.columns.values,
                    color_continuous_scale='Viridis',
                    aspect="auto")
    fig.update_xaxes(side="top")
    fig.update_layout(layout)

    import seaborn as sns
    import matplotlib.pyplot as plt

    import matplotlib.pyplot as plt

    hm = sns.heatmap(corr, annot = True)

    plt.show()
    return scatterMatrix