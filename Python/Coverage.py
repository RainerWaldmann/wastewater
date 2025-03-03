import pandas
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import math
import settings
import numpy as np
from joblib import Parallel, delayed
import plotly.io as pio
from multiprocessing import Pool
import time
import locale




def create_one_coverage_trace(onesample_data, i):
        def downsample_data(data, n):
            """Downsample a pandas DataFrame by taking every nth row."""
            return data.iloc[::n, :]

        def __cleanUpCoverageData(dat: pandas.DataFrame, column: str) -> pandas.DataFrame:
            """removes values from depth data that are too close to previous and next value
            --> no need to have lots of identical points in depth plot - will slow down html
            checks whether preceding and following data point are equal to current value"""
            x = dat if settings.downsampleFactorDepths is None or settings.downsampleFactorDepths <= 1 else downsample_data(
                dat, settings.downsampleFactorDepths)
            # after = x.shift(1)[column]
            # before = x.shift(-1)[column]
            # position = x[column]
            # r = [False if
            #      pos == aft == bef
            #      else True for aft, bef, pos in zip(after, before, position)]
            # return x[r]
            return x
        d = __cleanUpCoverageData(onesample_data.depths, 'COUNT')

        if settings.plotDepthsLogY:
            d['COUNT'] = np.where(d['COUNT'] == 1, np.nan, d['COUNT'])
            #d['COUNT'] = np.where(d['COUNT'] == 0, np.nan, d['COUNT'])
        trace = go.Scatter(y=d['COUNT'],
                           x=d['POS'],
                           line=dict(
                               color='Navy',  # 'rgb(160, 10, 10)',
                               width=1.5),
                           textposition="top center"
                           )
        return onesample_data.sample, trace

def plotCoverage(sample_data):
    """plot coverage for entire sequence and spike"""

    def __calcCoverage(d):
        """calculates % coverage with >= settings.minDepth reads from depth tsv data"""
        count = (d.loc[settings.maxAmplifiedRange[0]:settings.maxAmplifiedRange[1], 'COUNT'] >= settings.minDepth).sum()
        countSpike = (d.loc[settings.spikeRange[0]:settings.spikeRange[1], 'COUNT'] >= settings.minDepth).sum()
        # calculate relative SE of coverage
        # need depth with values for all positions
        total_bases = 29870
        spike_bases = 3822
        all_positions = pandas.DataFrame(
            {'POS': np.arange(settings.maxAmplifiedRange[0], settings.maxAmplifiedRange[1] + 1)}).set_index('POS')
        merged_df = pandas.merge(d, all_positions, left_index=True, right_index=True, how='outer')
        meancov = merged_df['COUNT'].mean()
        sdcov = merged_df['COUNT'].std()
        totalcov = count * 100 / total_bases
        spikecov = countSpike * 100 / spike_bases
        se = 100 * sdcov / meancov  if meancov > 0 else 100
        return totalcov, spikecov, se

    def __createCoverageDataTSV(sample_data):
        """creates a tsv file with coverage data for all samples"""
        df = sample_data[0].depths[['COUNT']].rename(columns={'COUNT': sample_data[0].sample})
        for sample in sample_data[1:]:
            # Join the 'TOTAL_DP' column from the sample's depth data to the new DataFrame
            # Use the sample's name as the column name in the new DataFrame
            df = df.join(sample.depths[['COUNT']].rename(columns={'COUNT': sample.sample}))
            df.to_csv(settings.rootDir + "/CoverageAllPositions.tsv", sep='\t',
                       decimal=locale.localeconv()['decimal_point'])
            df.loc[::10].to_csv(settings.rootDir + "/CoveragePositionsReduced.tsv", sep='\t',
                       decimal=locale.localeconv()['decimal_point'])

    start_time = time.time()
    __createCoverageDataTSV(sample_data)
    samples = [x.sample for x in sample_data]
    nCols = 3
    nRows = math.ceil(len(sample_data) / nCols)
    depthplots_fig = make_subplots(rows=nRows, cols=nCols, start_cell="top-left",
                                   subplot_titles=samples)

    coveragedf = pandas.DataFrame(columns=['Coverage', 'SpikeCoverage', 'CV%']) # holds coverage data for all samples

    with Pool() as p: # create a dict with sample as key
        results = {key: value for key, value in p.starmap(create_one_coverage_trace,
                                                          [(onesample_data, i) for i, onesample_data in
                                                           enumerate(sample_data)])}
    for i in range(len(sample_data)):
        trace = results[sample_data[i].sample]
        row = math.ceil((i + 1) / nCols)
        col = (i % nCols) + 1
        depthplots_fig.add_trace(trace, row=row, col=col)
        coverage = __calcCoverage(sample_data[i].depths)
        annotation_text = "<b>" + sample_data[i].sample + "</b>" + "<br>" + 'Coverage:' + '{:.2f}'.format(
            coverage[0]) + '%' + "<br>" + "CV: " + '{:.2f}'.format(
            coverage[2]) + '%' + "<br>" + 'Spike Coverage:' + '{:.2f}'.format(coverage[1]) + '%'
        depthplots_fig['layout']['annotations'][i].update(text=annotation_text)
        coveragedf.loc[sample_data[i].sample] = {'Coverage': coverage[0], 'SpikeCoverage': coverage[1],
                           'Coverage CV%': coverage[2] * 100}

    depthplots_fig.add_vrect(
        x0=settings.spikeRange[0], x1=settings.spikeRange[1],
        fillcolor="LightYellow", opacity=0.5,
        layer="below", line_width=0
    )

    # depthplots_fig.update_traces(marker=dict(color="RoyalBlue"))
    depthplots_fig.update_layout(width=800 * nCols,
                                 height=600 * nRows, title_text="Sequencing Depth (nReads per position)",
                                 yaxis=dict(
                                     showgrid=True
                                 ), title=""
                                 )
    depthplots_fig.update_layout(showlegend=False, hovermode=False)
    depthplots_fig.update_xaxes(tickfont=dict(size=18))
    depthplots_fig.update_yaxes(tickfont=dict(size=18))
    if settings.plotDepthsLogY:
        depthplots_fig.update_yaxes(type="log")
    print("Mean Coverage = ", coveragedf['Coverage'].mean() , "  SE = ",
          coveragedf['Coverage'].std() / np.sqrt(len(samples)), " N= ",
          len(samples))
    print("Mean Spike Coverage = ", coveragedf['SpikeCoverage'].mean(), "  SE = ",
          coveragedf['SpikeCoverage'].std() / np.sqrt(len(samples)), " N= ",
          len(samples))

    coveragedf.to_csv(settings.rootDir + "/coverageMeans.tsv", sep='\t', decimal=locale.localeconv()['decimal_point'])

    if settings.depthPDFfilename is not None:
        pio.write_image(depthplots_fig, settings.rootDir + "/" + settings.depthPDFfilename)

    print(f"--------- Sequencing depths plots took {time.time() - start_time} seconds to complete.")
    return depthplots_fig


def plotCoverageViolinPlots(sample_list):
    """plot coverage violin plots for entire sequence and spike"""
    # https://stackoverflow.com/questions/63082393/python-multiple-split-violine-plot-overlayed
    # https://plotly.com/python/subplots/

    violinsPerRow = 5
    nRows = math.ceil(len(sample_list) / violinsPerRow)
    violinsPlot_fig = make_subplots(rows=nRows, cols=violinsPerRow, start_cell="top-left")
    for i in range(0, len(sample_list)):
        row = math.ceil((i + 1) / violinsPerRow)
        col = i % violinsPerRow + 1
        sample = "<b>" + sample_list[i].sample + "</b>"
        violinsPlot_fig.add_trace(go.Violin(y=sample_list[i].depths['COUNT'], meanline_visible=True,
                                            x0= sample,
                                            box_visible=True,
                                            points=False,
                                            width=1,
                                            showlegend=i == 0,
                                            name='Sars-Cov-2',
                                            legendgroup='Sars-Cov-2', scalegroup='Sars-Cov-2',
                                            side='negative',
                                            line_color='blue',
                                            hoverinfo='skip'
                                            ),
                                  row=row, col=col
                                  )

        mask = (sample_list[i].depths['POS'] > settings.spikeRange[0]) & (
                    sample_list[i].depths['POS'] < settings.spikeRange[1])
        spikedepths = sample_list[i].depths.loc[mask, 'COUNT']

        violinsPlot_fig.add_trace(go.Violin(y=spikedepths, meanline_visible=True,
                                            showlegend=i == 0,
                                            legendgroup='Spike',
                                            scalegroup='Spike',
                                            x0=sample,
                                            box_visible=True,
                                            points=False,
                                            width=1,
                                            name='Spike',
                                            side='positive',
                                            line_color='orange',
                                            hoverinfo='skip'),
                                  row=row, col=col
                                  )

    violinsPlot_fig.update_traces(meanline_visible=True)
    violinsPlot_fig.update_layout(violingap=0, violinmode='overlay', width=violinsPerRow * 350,
                                  height=math.ceil(len(sample_list) / violinsPerRow) * 400,
                                  title_text="",
                                  yaxis=dict(
                                      showgrid=True
                                  )
                                  )
    # fig.update_yaxes(type="log", range=[-2,5]) # log range: 10^0=1, 10^5=100000
    print("DONE plotCoverageViolins")
    return violinsPlot_fig


def doShannon(sample_list):
    """Plot Shannon Entropy"""
    if not settings.doShannon:
        return None
    # def __calcMeanShannon(d):
    #     shannon = (d.loc[:, 'Shannon']).sum() / 29870
    #     shannonSpike = (d[21562:25384].loc[:, 'Shannon']).sum() / 3822
    #     return shannon, shannonSpike

    def __calcMeanShannon(d):
        shannon = (d.loc[:, 'Shannon']).sum() / 29870
        shannonSpike = (d[21562:25384].loc[:, 'Shannon']).sum() / 3822
        return shannon, shannonSpike

    def __cleanUpShannonData(x: pandas.DataFrame, column: str) -> pandas.DataFrame:
        """removes values from  data that are too close to previous and next value
        --> no need to have lots of identical points in plot """
        return x[[False if zero != zero  or math.fabs(zero - plusone) < 0.03 and math.fabs(zero - minusone) < 0.03 else True
             for zero,plusone,minusone in zip(x[column],x.shift(1)[column],x.shift(-1)[column])]]


    data = {}
    for s in sample_list:
        shannon, shannonSpike = __calcMeanShannon(s.depths)
        data[s.sample] = [shannon, shannonSpike]  # Use sample name as key

    # Create DataFrame from the dictionary
    meanShannons = pandas.DataFrame(data).transpose()
    meanShannons.columns = ['MeanShannon', 'MeanSpikeShannon']  # Rename columns
    if settings.useDateAxis:
        dates = [x.date for x in sample_list]
        meanShannons['Date'] = dates
    meanShannons.to_csv(settings.rootDir + "/Shannon.tsv", sep='\t', decimal=locale.localeconv()['decimal_point'])
    if settings.plotShannon:
        samples = [x.sample for x in sample_list] # truncate names at 20 chars ???
        nCols = 3
        nRows = math.ceil(len(sample_list) / nCols)
        # "subplot_titles=samples" in line below: titles are actually defined later. However If I put nothing here, shannonplots_fig['layout']['annotations'][i].update below throws an exception.
        shannonplots_fig = make_subplots(rows=nRows, cols=nCols, start_cell="top-left", subplot_titles=samples)
        for i, sample in enumerate(sample_list):
            row = math.ceil((i + 1) / nCols)
            # row = int((i)//nCols + 1)
            col = (i % nCols) + 1
            #print("Shannon Length df before filter:" + str(len(sample.depths)))
            df = __cleanUpShannonData(sample.depths, 'Shannon')
            #print("Shannon Length df after filter:" + str(len(df)))
            if settings.shannonPlotSpikeOnly:
                df =df[settings.spikeRange[0]:settings.spikeRange[1]]
            trace = go.Scatter(y=df['Shannon'],
                               x=df['POS'],
                               line=dict(
                                   color='Navy',  # 'rgb(160, 10, 10)',
                                   width=1.5
                               ),
                               hovertemplate="Pos: %{x}<br>Entropy: %{y}"
                               )

            shannonplots_fig.add_vrect(
                x0=settings.spikeRange[0], x1=settings.spikeRange[1],
                fillcolor="LightYellow", opacity=0.5,
                layer="below", line_width=0
            )
            shannonplots_fig.add_trace(trace, row=row, col=col)
            #shannon = __calcMeanShannon(sample.depths)
            #shannon_list.append(shannon[0])
            #spike_shannon_list.append(shannon[1])
            shannonplots_fig['layout']['annotations'][i].update(text= "<b>" + (samples[i] + "</b>" + "<br>" +
                                                                      'Mean shannon:' + '{:.5f}'.format(meanShannons.at[sample.sample, 'MeanShannon']) + "<br>" +
                                                                               'Mean spike shannon:' + '{:.5f}'.format(meanShannons.at[sample.sample, 'MeanSpikeShannon']) + "<br>" +
                                                                               'Mean shannon diff spike/tot:' + '{:.5f}'.format(meanShannons.at[sample.sample, 'MeanSpikeShannon'] - meanShannons.at[sample.sample, 'MeanShannon'])
                                                                               ))
        shannonplots_fig.update_layout(width=800 * nCols,
                                       height=600 * nRows, title_text="Shannon Entropy",
                                       yaxis=dict(
                                         showgrid=True),
                                       title="",
                                       hovermode='closest',
                                       showlegend=False)
        # shannonplots_fig.update_yaxes(type="log")
        return shannonplots_fig
    else:
        return None


def plotShannonViolins(sample_list):
    """NOT USED"""
    violinsPerRow = 5
    nRows = math.ceil(len(sample_list) / violinsPerRow)
    violinsPlot_fig = make_subplots(rows=nRows, cols=violinsPerRow, start_cell="top-left")
    for i in range(0, len(sample_list)):
        # row = int((i)//violinsPerRow + 1)
        row = math.ceil((i + 1) / violinsPerRow)
        col = i % violinsPerRow + 1
        # sample_list[i].depths.loc[:, 'Shannon']
        df = sample_list[i].depths
        violinsPlot_fig.add_trace(go.Violin(y=df['Shannon'],
                                            meanline_visible=True,
                                            x0=sample_list[i].sample,
                                            box_visible=True,
                                            points=False,
                                            width=1,
                                            showlegend=i == 0,
                                            name='Sars-Cov-2',
                                            legendgroup='Sars-Cov-2', scalegroup='Sars-Cov-2',
                                            side='negative',
                                            line_color='blue'
                                            #,customdata = df['POS']
                                            ,hovertemplate = "Position: %{y}: <br>Entropy: %{y}"
                                            # ,hoverinfo='skip'
                                            ),
                                  row=row, col=col
                                  )
        mask = (sample_list[i].depths.iloc[:, 1] > settings.spikeRange[0]) & (
                    sample_list[i].depths.iloc[:, 1] < settings.spikeRange[1])
        spikedepths = sample_list[i].depths.loc[mask, 'Shannon']
        
        violinsPlot_fig.add_trace(go.Violin(y=spikedepths, meanline_visible=True,
                                            showlegend=i == 0,
                                            legendgroup='Spike', scalegroup='Spike', x0=sample_list[i].sample,
                                            box_visible=True,
                                            points=False,
                                            width=1,
                                            name='Spike',
                                            side='positive',
                                            line_color='orange',
                                            hovertemplate="Position: %{y}: <br>Entropy: %{y}"
                                            ),
                                  row=row, col=col
                                  )

    violinsPlot_fig.update_traces(meanline_visible=True)
    violinsPlot_fig.update_layout(violingap=0, violinmode='overlay', width=violinsPerRow * 350,
                                  height=math.ceil(len(sample_list) / violinsPerRow) * 400,
                                  title_text="",
                                  yaxis=dict(
                                      showgrid=True
                                  )
                                  )
    violinsPlot_fig.update_yaxes(type="log")
    violinsPlot_fig.update_yaxes(range=[-3, -1])
    return violinsPlot_fig


def getShannonBoxPlots(sample_list):
    """Generate Shannon Entropy Box Plots"""
    if not settings.plotShannon:
        return None
    boxPlotsPerRow = 5
    nRows = math.ceil(len(sample_list) / boxPlotsPerRow)
    boxPlot_fig = make_subplots(rows=nRows, cols=boxPlotsPerRow, start_cell="top-left",
                                subplot_titles=["<b>" + s.sample + "</b>" for s in sample_list])

    for i in range(0, len(sample_list)):
        row = math.ceil((i + 1) / boxPlotsPerRow)
        col = i % boxPlotsPerRow + 1
        s1 = sample_list[i].depths
        s2 = sample_list[i].depths.iloc[settings.spikeRange[0]:settings.spikeRange[1]]
        boxPlot_fig.add_trace(go.Box(y=s1['Shannon'], marker_color='blue',
                                            name="All("+ "{:.5f}".format(np.mean(sample_list[i].depths['Shannon'])) +")",
                                            customdata = s1['POS'],
                                            hovertemplate = "Position: %{customdata}: <br>Entropy: %{y}" ),
                                        row=row,
                                        col=col)
        boxPlot_fig.add_trace(go.Box(y=s2['Shannon'], marker_color='orange', name ="Spike(" + "{:.5f}".format(np.mean(sample_list[i].depths['Shannon'].iloc[settings.spikeRange[0]:settings.spikeRange[1]])) +")",
                                            customdata = s2['POS'],
                                            hovertemplate = "Position: %{customdata}: <br>Entropy: %{y}"),
                                        row=row,
                                        col=col)

    boxPlot_fig.update_layout(width=boxPlotsPerRow * 350,
                                  height=math.ceil(len(sample_list) / boxPlotsPerRow) * 600,
                                  title_text="",
                                  showlegend=False,
                                  yaxis=dict(
                                      showgrid=True)
                                  )
    return boxPlot_fig

