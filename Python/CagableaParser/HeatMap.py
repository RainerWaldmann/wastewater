import plotly.graph_objects as go
import plotly.io as plotly_io
import pandas as pd
import re
import settings
import variants.globals as vg
import ReadData
import Data
import variants.MutInfo as mi
from plotly.subplots import make_subplots
import dash_bio


def __getHeatMapHoverText(df : pd.DataFrame, sample_list : list[Data.OneBCdata]) -> pd.DataFrame:
    """generates hover text for heatmap . Returns dataframe in same format as input with hoverdata"""
    def __format_intvalue(value):
        return str(round(value)) if pd.notna(value) else "NA"

    def __generateOneHoverLine(x : pd.Series, sample: Data.OneBCdata) -> pd.Series:
        """Does job for one line"""
        ret = f"<b>Sample: {sample.sample}</b><br>"
        ret += f"<b>Mutation: </b>{x['REF']}{x['POS']}{x['ALT']}"
        if (len(x["ALT"]) == 1):
            ret += mi.MutInfoSubst.getMutInfoFromPosMut(x["POS"], x["ALT"]).getAAmutstringForDFindex()
        ret += "<br>"
        ret += f"<b>Frequency: </b>{round(x['ALT_FREQ_'+ sample.sample],3) if pd.notna(x['ALT_FREQ_'+ sample.sample]) else 'NA'}<br>"
        ret += f"REF/ALT quality: {__format_intvalue(x['REF_QUAL_' + sample.sample])}/{__format_intvalue(x['ALT_QUAL_' + sample.sample])}"
        ret += f"<br>REF/ALT depth (FWD,REV): ({__format_intvalue(x['REF_DP_' + sample.sample] - x['REF_RV_' + sample.sample])},{__format_intvalue(x['REF_RV_' + sample.sample])})"
        ret += f")/({__format_intvalue(x['ALT_DP_' + sample.sample] - x['ALT_RV_' + sample.sample])},{__format_intvalue(x['ALT_RV_' + sample.sample])})"
        return ret

    dfp_hover = pd.DataFrame()
    for sample in sample_list:
        dfp_hover[sample.sample] = df.apply(lambda x: __generateOneHoverLine(x, sample), axis=1)
    return dfp_hover


def plotHeatMap( mergedData: pd.DataFrame, sample_list: list[Data.OneBCdata]) -> go.Figure:
    """gets heatmap figure""" 
     # get dataframes with just frequencies
    dfpForHeatMaps = ReadData.getFrequencyColumnsOnly(mergedData, sample_list, None)   
    dfp_hover = __getHeatMapHoverText(mergedData, sample_list)
    newIndexlist = list()
    for ind in dfpForHeatMaps.index:
        for voc in vg.variantdict:
            vocMutList = vg.variantdict[voc].getNucAcidMutList()
            if any(sub in ind for sub in vocMutList):
                ind = voc + ":" + ind
        newIndexlist.append(ind)

    dfA = dfpForHeatMaps.set_index(pd.Index(newIndexlist))  # set new index

    # index = [s.split(",")[0] for s in index] ,\w+
    # print([s.split(",")[0] for s in dfp.index])
    hoverinfo = '<b>%{customdata}</b>' if settings.hoverinfo_Onheatmaps else None
    hovermode = 'closest' if settings.hoverinfo_Onheatmaps else False #possible hovermode ['x', 'y', 'closest', False, 'x unified', 'y unified']
    heatmap_fig = go.Figure(data=go.Heatmap(
        z=dfA,
        y=dfA.index,
        x=dfA.columns, colorscale=[(0, "blue"), (0.2, "red"), (1, "yellow")],
        xgap=2,
        ygap=2,
        colorbar=dict(title='Mut Freq'),
        customdata=dfp_hover,
        #hoverinfo='skip',# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        hovertemplate= hoverinfo
        #            '<i>Site</i>: %{x}'+
        #            '<br>Mutation: %{y}'+
        #            '<br><b>Frequency= </b>: %{z}' +


    ),
        layout=go.Layout(plot_bgcolor='rgb(100, 100,100)', height=len(dfA) * 20, width=dfA.shape[1] * 50 + 400,
                hovermode=hovermode
                         # title_text="Plotly heatmap"
                         )
    )
    # custom hover:  https://plotly.com/python/hover-text-and-formatting/
    # custom hover: https://stackoverflow.com/questions/45569092/plotly-python-heatmap-change-hovertext-x-y-z
    # custom hover: https://community.plotly.com/t/heatmap-changing-x-y-and-z-label-on-tooltip/23588/4
    #    hovertemplate =
    #    '<i>Price</i>: x'+
    #    '<br><b>Frequency= </b>: y<br>',
    # + '<b>%{text}</b>',
    # text = ['Custom text {}'.format(i + 1) for i in range(5)],


    # ,
    # layout=go.Layout(plot_bgcolor=’#777777’)))
    #  layout=plot_bgcolor='rgb(230, 230,230)')))
    if settings.doPlot:
        heatmap_fig.show()

    # heatmap_fig.write_html(rootDir.as_posix() + "/heatmap.html")
    #plotly_io.write_image(heatmap_fig, settings.rootDir + "/heatmap.pdf", format='pdf')
    return heatmap_fig


def plotHeatMapVariantsOnly(mergedData: pd.DataFrame, sample_list: list[Data.OneBCdata] ): 
    dfpForHeatMaps = ReadData.getFrequencyColumnsOnly(mergedData, sample_list, None)
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
    nVariant=0
    for voc in vg.variantdict:
        filt = [any([i in x for i in vg.variantdict[voc].getNucAcidMutList()]) for x in df.index]
#    for voc in VariantDefinitions.VOCsLineageDef:
        #filt = df.apply(lambda x: any([i in x.index for i in VariantDefinitions.VOCsAllDefiningMuts[voc]]),axis = 1)
 #       filt = [any([i in x for i in  VariantDefinitions.VOCsLineageDef[voc]]) for x in df.index]
        d = df[filt]
        if len(d) != 0 :
            vocDFs[voc] = d
            totRows += d.shape[0]
            nVariant+=1

    row_heights = [vocDFs[v].shape[0]/totRows for v in vocDFs]

    #layout = go.Layout(plot_bgcolor='rgb(100, 100,100)', height=1000, width=1000)
    heatmap_fig = make_subplots(rows=len(vocDFs), cols=1, start_cell="top-left",
                                   subplot_titles=[v for v in vocDFs],row_heights = row_heights,
                               vertical_spacing = 0.02)
    i=1
    for voc in vocDFs:
        f = go.Heatmap(
        z=vocDFs[voc],
        y=vocDFs[voc].index,
        x=vocDFs[voc].columns, colorscale=[(0, "blue"), (0.2, "red"), (1, "yellow")],
        xgap=2,
        ygap=2,
        coloraxis="coloraxis"

        #,showlegend=(i == len(vocDFs))
        )
        #f.update_xaxes(visible= (i == len(vocDFs)))
        heatmap_fig.add_trace(f,row=i, col=1)
        i += 1

    colorbar_trace=go.Scatter(x=[None],y=[None],mode='markers',marker=dict(
                                 colorscale=[(0, "blue"), (0.2, "red"), (1, "yellow")],
                                 showscale=True,
                                 cmin=0,
                                 cmax=1,
                                 colorbar=dict(thickness=15, tickvals=[0,0.2,0.4,0.6,0.8,1], outlinewidth=0,orientation="h")
                             ),
                             hoverinfo='none'
                            )
    heatmap_fig.add_trace(colorbar_trace)

#    if settings.doPlot:
    heatmap_fig.update_layout(width=dfpForHeatMaps.shape[1]*80 + 200,height=totRows*100 + len(vocDFs)*80 + 200,coloraxis=dict(colorscale=[(0, "blue"), (0.2, "red"), (1, "yellow")]), showlegend=False)
    heatmap_fig.update_layout(width=dfpForHeatMaps.shape[1]*80 + 200, height=2000)
    heatmap_fig.update_layout(font=dict(size=10,color="RebeccaPurple"))
    heatmap_fig.update_xaxes(visible=False)
    heatmap_fig.update_xaxes(visible=True,col=1,row=nVariant)
    heatmap_fig.update_coloraxes(showscale=False)
    #heatmap_fig.show()

    # heatmap_fig.write_html(rootDir.as_posix() + "/heatmap.html")
    plotly_io.write_image(heatmap_fig, settings.rootDir + "/heatmapVariants.pdf", format='pdf')
    return heatmap_fig


transpose = True
def clusterMap(mergedData:pd.DataFrame, sample_list: list[Data.OneBCdata]):
    dfpForHeatMap = ReadData.getFrequencyColumnsOnly(mergedData, sample_list, None)
    dfz = dfpForHeatMap.transpose()  if  transpose else dfpForHeatMap
    columns = list(dfz.columns)
    rows = list(dfz.index)
    hovermode = 'closest' if settings.hoverinfo_Onheatmaps else False  # possible hovermode ['x', 'y', 'closest', False, 'x unified', 'y unified']
    cluster_fig = dash_bio.Clustergram(
        center_values= False,
        data=dfz.fillna(-0.0001).loc[rows].values,
        row_labels=rows,
        column_labels=columns,
        height= dfz.shape[0]*40 +350 if transpose else dfpForHeatMap.shape[0]*20 ,
        width=dfz.shape[1]*20 if transpose else dfpForHeatMap.shape[1]*100 + 150,
        color_threshold={
            'row': 0.3,
            'col': 0.3
        },
       color_list={
    #        'row': ['#636EFA', '#00CC96', '#19D3F3'],
    #        'col': ['#AB63FA', '#EF553B'],
            'bg': 'white'
        },
        color_map= [[0.0, 'grey'],[0.00011, 'rgb(0,0,100)'], [0.3, 'rgb(255,0,0)'], [1.0, 'rgb(255,255,0)']],
    #    color_map= [
    #        [0.0, "blue"],
    #        [1.0, "red"]
    #    ],
        #color_map= [(0.0, "black"),(0.05, "blue"), (0.2, "red"), (1, "yellow")],

        line_width=1.5,
        display_ratio=[0.1, 0.4] if transpose else [0.4, 0.1],
        paper_bg_color = 'rgb(240, 240,240)'
    )
    cluster_fig.update_traces(colorbar_orientation='h', selector=dict(type='heatmap'),
                              xgap = 2,
                              ygap = 2,
                     #customdata= hoverDictionary,
                     #hovertemplate ='<i>Site</i>: %{customdata[%{z}]}' + '<br><b>Frequency= </b>: %{z}'
                     hovertemplate = '<b>Site= </b>: %{y}<br>' +'<i>Mutation</i>: %{x}' +  '<br><b>Frequency= </b>: %{z}'
                     )
    cluster_fig.update_layout(hovermode = hovermode)
    #cluster_fig.update_layout(plot_bgcolor='rgb(100, 100,100)')
    if settings.doPlot:
        cluster_fig.show()
    #cluster_fig.write_html(rootDir.as_posix() + "/heatmap2.html")
    #cluster_fig.write_image(rootDir.as_posix() + "/heatmap2.pdf")
    #plotly_io.write_image(cluster_fig, settings.rootDir + "/clustermap.pdf", format='pdf')
    return cluster_fig