import plotly.graph_objects as go
import plotly.io as plotly_io
import pandas as pd
import re
import settings
import variants.globals as vg


excludedMuts = []

def plotHeatMap(dfpForHeatMaps):
    index = [re.sub(r'[N]{7,}', "N7+", s) for s in dfpForHeatMaps.index]  # replace >= 7N by N7+
    index = [re.sub(r',.+', "", s) for s in index]  # remove everything after the ","
    df = dfpForHeatMaps.set_index(pd.Index(index))
    # newIndexlist = list()
    # for ind in dfA.index:
    #     for voc in VariantDefinitions.VOCsAllDefiningMuts:
    #         vocMutList = VariantDefinitions.VOCsAllDefiningMuts[voc]
    #         if any(sub in ind for sub in vocMutList):
    #             ind = voc + ":" + ind
    #     newIndexlist.append(ind)

    #dfA.set_index(pd.Index(newIndexlist), inplace=True)  # set new index

    #https://stackoverflow.com/questions/67700318/check-if-series-contains-any-element-from-a-list
    #dict of dataframes for VOCs
    vocDFs = dict()
    totRows = 0
    for voc  in vg.variantdict:
        #filt = df.apply(lambda x: any([i in x.index for i in VariantDefinitions.VOCsAllDefiningMuts[voc]]),axis = 1)
        filt = [any([i in x for i in  vg.variantdict[voc].getNucAcidMutList()]) for x in df.index]
        d = df[filt]
        if len(d) != 0 :
            vocDFs[voc] = d
            totRows += d.shape[0]
            #print("+++++++++++++++++++++\n",d)

    #MERGE voc dataframes
    dfVOCs = pd.concat(list(vocDFs.values())[::-1])
    #print("$$$$$$$$$$$$$$$$$$$$$$$\n",dfVOCs)
    dfNoVOC = df[~df.index.isin(dfVOCs.index)]
    #print("\n*************************\n",dfNoVOC)
    dfVOCs = pd.concat([dfNoVOC, dfVOCs])
    #print("\n*********$$$$$$$$$$**********\n", dfVOCs)
    #print("\n*********HHHHHHHHHHH**********\n", dfVOCs.index)

    heatmap_fig = go.Figure(data=go.Heatmap(
        z=dfVOCs,
        y=dfVOCs.index,
        x=dfVOCs.columns, colorscale=[(0, "blue"), (0.2, "red"), (1, "yellow")],
        xgap=2,
        ygap=2,
        colorbar=dict(title='Mut Freq')
      ),
        layout=go.Layout(plot_bgcolor='rgb(100, 100,100)', height=dfVOCs.shape[0] * 20 + 200, width=dfVOCs.shape[1] * 50 + 400
                         # title_text="Plotly heatmap"
                         )
    )


    # heatmap_fig.write_html(rootDir.as_posix() + "/heatmap.html")
    heatmap_fig.show()
    plotly_io.write_image(heatmap_fig, settings.rootDir.as_posix() + "/heatmapVariantsFirst.pdf", format='pdf')
    return heatmap_fig