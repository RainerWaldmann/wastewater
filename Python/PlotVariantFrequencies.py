
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as plotly_io
import settings
import variants
from variants import OneVariant
from itertools import groupby
from scipy.signal import savgol_filter
import scipy
import math
import numpy


class VariantsDataForHTML:
    """Just to group various data used in html. Avoids using too many seperate parameters for Jinja"""
    def __init__(self, variantPieList,variantsDetailedCounts, variantsdetailedCountsMask, variantDataMeans,variantHistoFigHTML, linechartHTML=None):
        self.variantPieList = variantPieList  # the list of pie charts
        self.variantsDetailedCounts = variantsDetailedCounts  # the Var frequencies for each mutation for each variant for the different samples
        self.variantsdetailedCountsMask = variantsdetailedCountsMask  # same structure as variantsDetailedCounts, contains just true or false to indicate whether value is an outlier
        self.variantDataMeans = variantDataMeans  # contains VariantData object for each variant and sample -> can retrieve means and se
        self. variantHistoFigHTML = variantHistoFigHTML  # boxplot html for each variant and site showing the median .. for all mutation defining the variant


def plotVariants(sampledata_list, detailedCounts:pd.DataFrame, detailedCountsMask : pd.DataFrame, variantDataMeansSD : pd.DataFrame, variantDataMeansSDUncorrected : pd.DataFrame) -> VariantsDataForHTML:
    """The function called from main"""

    # ****************    PIECHARTS ****************************************************************
    variantPieList =__generatePieChartsHTML(sampledata_list, variantDataMeansSD)

    #######   HISTOGRAM ####################################################################
    usedVariantDataMeansSD = variantDataMeansSD if settings.variantHistogram_sumChildsNotExceedParentCorrection else variantDataMeansSDUncorrected
    figVariantHisto = __generateHistograms(usedVariantDataMeansSD)
    variantDataTr = usedVariantDataMeansSD.transform(lambda x: [v.mean for v in x])
    if settings.plot_variant_freqs_line_chart:
        __generateStackedAreaChart(variantDataTr, sampledata_list)

    # add color dots for HTML to detailed counts
    for voc in detailedCounts:
        mutinfoList = [[d for d in variants.globals.variantdict[voc].data if d.getNAmutstring() in x][0] for x in
                       detailedCounts[voc].index]  # get corresponding MutIfo Objects for the mutation
        colorDotList = ["".join([x.getHTMLforColoredDot() for x in
                                 l.listOfOtherVariantsMutIsFound]) if l.listOfOtherVariantsMutIsFound is not None else ""
                        for l in mutinfoList]  # get the dots to print after mut
        detailedCounts[voc]["dots"] = colorDotList

    return VariantsDataForHTML(variantPieList, detailedCounts, detailedCountsMask, variantDataMeansSD,
                               figVariantHisto.to_html(full_html=False, include_plotlyjs='cdn'))


def __generateHistograms(usedVariantDataMeansSD: pd.DataFrame)-> go.Figure:
    """generates histograms for variant frequencies"""
    variantDataTr = usedVariantDataMeansSD.transform(lambda x: [v.mean for v in x])
    variantDataErr = usedVariantDataMeansSD.transform(lambda x: [v.se for v in x])
 # list of variants for histogram, only variants where print is enabled and where any of values is > 0
    variantsForHisto = [v for v in variants.globals.variantdict if
                        variants.globals.variantdict[v].histogramGroup and variantDataTr.loc[v].gt(0).any() and
                        variants.globals.variantdict[v].doPrint]
    # group variants by histogram bar, list must be sorted for groupby to work --> sort by variants.globals.variantdict[x].histogramOrderId
    groupedVariants = {key: list(value) for key, value in
                       groupby(sorted(variantsForHisto, key=lambda x: variants.globals.variantdict[x].histogramOrderId),
                               lambda x: variants.globals.variantdict[x].histogramOrderId)}
    data = []  # the parts of the histogram
    for offsetGroupIndex, oneGroup in groupedVariants.items():
        def defineOffset(variant: variants.OneVariant.OneVariant):  # get all childs and define offsets for histo bars in group
            """function that calculates the offset of childs within a bar, baseline if offset of parent"""
            childs = [v for v in oneGroup if variants.globals.variantdict[v].parent is not None and variants.globals.variantdict[v].parent.name is variant]
            if childs:
                childOffset = oneGroupWithOffsets[variant]  # first child offset equal parent offset
                for onechild in childs:
                    oneGroupWithOffsets[onechild] = childOffset
                    # groupWithFinishedStatus[onechild] = True
                    childOffset = (childOffset + variantDataTr.loc[
                        onechild])  # .copy() # increase child offset for next round by frequency of current child
                    defineOffset(onechild)  # treat childs of current child

        oneGroup = sorted(oneGroup, key=lambda x: variants.globals.variantdict[x].histogramOrderIdWithinBar)
        oneGroupWithOffsets = {v: 0 for v in oneGroup}  # holds variant name and offset of histogram bar

        # get top level variants (variants that don't have parent in current group)
        topNodes = [v for v in oneGroup if
                    (variants.globals.variantdict[v].parent is None or variants.globals.variantdict[
                        v].parent.name not in oneGroup)]
        if (len(topNodes) == 0):  # normally should not happen
            continue
        # set offset of first one to 0
        oneGroupWithOffsets[topNodes[0]] = 0
        for kk in range(1, len(topNodes)):
            oneGroupWithOffsets[topNodes[kk]] = (
                        oneGroupWithOffsets[topNodes[kk - 1]] + variantDataTr.loc[topNodes[kk - 1]])
        for variant in topNodes:
            # groupWithFinishedStatus[variant] = True
            defineOffset(variant)
        # legendgrouptitletext = "" if offsetGroupName.isdigit() else offsetGroupName

        for oneVariant in oneGroup:
            offsetGroupName = variants.globals.variantdict[oneVariant].histogramGroup
            data.append(go.Bar(
                # width = 0.5,
                # opacity=0.5,
                marker_line_color='black',
                legendgroup=variants.globals.variantdict[oneVariant].histogramGroup,
                legendgrouptitle_text=offsetGroupName,
                name=oneVariant,
                x=variantDataTr.columns,
                y=variantDataTr.loc[oneVariant],
                offsetgroup=offsetGroupIndex,
                base=oneGroupWithOffsets[oneVariant],  #
                marker_color=variants.globals.variantdict[oneVariant].color,
                error_y=dict(type='data', array=variantDataErr.loc[oneVariant]),
                # error_y_color=VariantDefinitions.VariantErrorbarColors[4]
                # text=oneVariant,
                # textposition='none',
                meta=[oneVariant],
                hovertemplate="<br>".join([
                    "Variant: %{meta[0]}",
                    "Frequency: %{y}"])
            ))

    figVariantHisto = go.Figure(data=data)
    figVariantHisto.update_layout(height=1000, width=len(variantDataTr.columns) * 100 + 300,
                                  paper_bgcolor='rgba(0,0,0,0)',
                                  plot_bgcolor='rgba(0,0,0,0)', barmode='group', bargap=0.4, bargroupgap=0.3,
                                  )
    # legendgrouptitle_font_size = 20,
    # legendgrouptitle_font_color="green", legendgrouptitle_font_family="Times New Roman")
    if settings.variant_hisogram_pdf:
        plotly_io.write_image(figVariantHisto, settings.rootDir + "/" + settings.variant_hisogram_pdf, format='pdf')
    return figVariantHisto


def __generateStackedAreaChart(variantFreqs: pd.DataFrame, samples):
    """generates stacked area chart for variant frequencies"""
    def get_ordered_variants(varForPlot: dict) -> list:
        """get list of variant nodes where each parent is followed by its childs - to define the order of traces in plot"""
        ordered_variants = []
        def dfs(variant: variants.OneVariant.OneVariant):
            ordered_variants.append(variant.name)
            #childs = variant.childs
            child_variants = [v for v in varForPlot.values() if v.parent is not None and v.parent.name == variant.name]
            for child in child_variants:
                dfs(child)

        top_nodes = [v for v in varForPlot.values() if v.parent is None or v.parent.doPrint == False or v.parent.name not in varForPlot]
        for v in top_nodes:
            dfs(v)
        return ordered_variants

    #__generateStackedAreaChartFromTSV()

    variants_for_plot = {k: v for k, v in variants.globals.variantdict.items() if
                         variantFreqs.loc[v.name].gt(0).any() and v.doPrint}
    ordered_variants = get_ordered_variants(variants_for_plot)
    variant_freqs = variantFreqs.reindex(ordered_variants)
    # Subtract frequencies of direct children from parents
    for variant in variants_for_plot.values():
        # Check if the variant has any child variants
        child_variants = [v.name for v in variants_for_plot.values() if v.parent == variant]
        if child_variants:
            # Calculate the sum of the frequencies of all child variants
            child_freqs_sum = variantFreqs.loc[child_variants].sum()
            # Subtract this sum from the frequency of the parent variant
            variant_freqs.loc[variant.name] -= child_freqs_sum
        # Following line Is Just for saving a dataframe we need for analysis!!!!!
        # variant_freqs.to_csv(settings.rootDir + "/VarFrequenciesLineChart.csv", index=True)

    data_per_variant = dict()
    for variant in variant_freqs.index:
        y = variant_freqs.loc[variant]
        # Apply Savitzky-Golay filter for smoothing
        if len(samples) > 7:
           y_smooth = y if not settings.smoothen_variantsLinePlot_with_filter else scipy.ndimage.gaussian_filter(y, sigma=settings.smoothen_variantsLinePlot_sigma)
           #y_smooth = numpy.clip(savgol_filter(y, window_length=4 if len(samples) > 30 else math.ceil(len(samples)/5), polyorder=2),0,1) if settings.smoothen_variantsLinePlot_with_filter else y
        else:
            y_smooth = y
        data_per_variant[variant] = y_smooth

        # Step 1: Loop through each x value (column) in data_per_variant
    for x in range(data_per_variant[next(iter(data_per_variant))].shape[0]):
        # Step 2: Calculate the sum of the smoothed values for all variants at this x value
        total_sum = sum(data_per_variant[variant][x] for variant in data_per_variant)
        if total_sum > 1:
            for variant in data_per_variant:
                data_per_variant[variant][x] /= total_sum

    fig = go.Figure()
    for variant in variant_freqs.index:
        x = variant_freqs.columns if not settings.useDateAxis else [sample.date for sample in samples]
        #y = variant_freqs.loc[variant]
        trace = go.Scatter(
            x=x,
            y=data_per_variant[variant],
            name=variant,
            mode="lines",
            stackgroup='stack',
            fill="tonexty",
            fillcolor=variants.globals.variantdict[variant].color,
            line_color=variants.globals.variantdict[variant].color
        )
        fig.add_trace(trace)

    fig.update_layout(
        width=4000,
        height=3000,
        title="SARS-CoV-2 Variant Frequencies Over Time",
        xaxis_title="Date",
        yaxis_title="Cumulative Frequency",
        legend_title="Variant",
    )

    if settings.useDateAxis:
        fig.update_xaxes(
            tickformat="%b\n%Y",  # Display month and year
            dtick="M1",  # Set tick every month
            tickangle=-45,  # Angle the ticks for better visibility
            tickmode="linear",  # Use linear mode for consistent tick intervals
            ticks="outside",  # Place the ticks outside the axis
            tickson="boundaries",  # Ensure ticks appear at month boundaries
            showgrid=True  # Optionally, show vertical grid lines for better alignment
        )
    plotly_io.write_image(fig, settings.rootDir + "/VariantFreqLineChart.pdf", format='pdf')


def __generatePieChartsHTML(samples, variantDataMeansSD ) -> dict:
    """generates pie charts for each variant and sample and returns dict with html for each sample"""
    if settings.plotVariantPies:
        variantsToPrint = [v for v in variants.globals.variantdict if variants.globals.variantdict[v].doPrint == True]
        pieData = variantDataMeansSD.loc[variantsToPrint]
        parents = [variants.globals.variantdict[v].parent.name if variants.globals.variantdict[v].parent and variants.globals.variantdict[v].parent.doPrint else "" for
                   v in variantsToPrint]
        variantcolors = [variants.globals.variantdict[v].color for v in variantsToPrint]
        variantPieDict = dict() # dict to hold pie charts for each sample
        for i in range(0, len(samples)):
            hoverdata = [[v.mean, v.se if v.se == v.se else "NA", v.n, v.nPosAboveZero] for v in
                         pieData[samples[i].sample]]

            trace = go.Sunburst(labels=variantsToPrint, parents=parents,
                                values=[v.mean for v in pieData[samples[i].sample]],
                                branchvalues="total"
                                ,customdata=hoverdata
                                ,hovertemplate='<b>Variant: %{label}<b><br>Frequency: %{value}<br>SE: %{customdata[1]:.3f}<br> %{customdata[3]} mutations of %{customdata[2]} above zero'
                                ,sort=settings.sortPieWedges
                                )
            figForList = go.Figure(trace)
            figForList.update_traces(
                marker=dict(colors=variantcolors, line=dict(color='#000000', width=0.5)))
            figForList.update_layout(height=350, width=350,
                                     margin=dict(t=10, l=0, r=0,
                                                 b=10))  # ,title=self.samples[i],margin = dict(t=0, l=0, r=0, b=0))
            variantPieDict[samples[i]] = figForList
        if settings.variantPiesPDFfilename is not None: # put variant pies into one figure to write them to PDF
            for sample, pie in variantPieDict.items():
                plotly_io.write_image(pie, settings.rootDir + "/" + "Pie_" + sample.sample + "_" + settings.variantPiesPDFfilename, format='pdf')
        return {sample:variantPie.to_html(full_html=False, include_plotlyjs='cdn') for sample,variantPie in variantPieDict.items()}
    else:
        return None


def __generateStackedAreaChartFromTSV():
    """generates stacked area chart for variant frequencies from TSV file (.e.g. patient data)"""

    def get_ordered_variants(varForPlot: dict) -> list:
        """get list of variant nodes where each parent is followed by its childs - to define the order of traces in plot"""
        ordered_variants = []

        def dfs(variant: variants.OneVariant.OneVariant):
            ordered_variants.append(variant.name)
            # childs = variant.childs
            child_variants = [v for v in varForPlot.values() if v.parent is not None and v.parent.name == variant.name]
            for child in child_variants:
                dfs(child)

        top_nodes = [v for v in varForPlot.values() if
                     v.parent is None or v.parent.doPrint == False or v.parent.name not in varForPlot]
        for v in top_nodes:
            dfs(v)
        return ordered_variants

    # Load the dataframe from a file
    file_path = r"C:\Users\raine\OneDrive - Université Côte d'Azur\Documents\LPUs\2024\Wastewater Nice\Epidemological Data\Patients\EMERGEN_06.csv"  # Replace with the actual file path
    variantFreqs = pd.read_csv(file_path, index_col=0)
    sample_ns = pd.read_csv(r"C:\Users\raine\OneDrive - Université Côte d'Azur\Documents\LPUs\2024\Wastewater Nice\Epidemological Data\Patients\EMERGEN_06_Ns.csv", index_col=0)
    # Extract data from sample_ns
    sample_ns = sample_ns.loc[:, ~sample_ns.columns.isna()]
    sample_ns = sample_ns[~sample_ns.index.isna()]
    # Drop columns where the name contains 'Unnamed' and all values are NaN
    sample_ns = sample_ns.loc[:, ~sample_ns.columns.str.contains('^Unnamed') | sample_ns.notna().any()]
    sample_ns.index = sample_ns.index.str.strip()
    if 'N' in sample_ns.index:
        sample_n_values = sample_ns.loc['N']
        sample_n_values.index = pd.to_datetime(sample_n_values.index, format="%m/%d/%Y")
    else:
        raise ValueError("No data with index 'N' found in sample_ns.")

    variants_for_plot = {k: v for k, v in variants.globals.variantdict.items() if
                         k in variantFreqs.index and variantFreqs.loc[k].gt(0).any() and v.doPrint}
   # variants_for_plot = {k: v for k, v in variants.globals.variantdict.items() if
   #                      variantFreqs.loc[v.name].gt(0).any() and v.doPrint}
    #protect against nan column names and indices
    variantFreqs = variantFreqs.loc[:, ~variantFreqs.columns.isna()]
    variantFreqs = variantFreqs[~variantFreqs.index.isna()]
    # Drop columns where the name contains 'Unnamed' and all values are NaN
    variantFreqs = variantFreqs.loc[:, ~variantFreqs.columns.str.contains('^Unnamed') | variantFreqs.notna().any()]

    ordered_variants = get_ordered_variants(variants_for_plot)
    variant_freqs = variantFreqs.reindex(ordered_variants)
    variant_freqs.columns = variantFreqs.columns.str.strip()
    variant_freqs.columns = pd.to_datetime(variant_freqs.columns, format="%m/%d/%Y")

    data_per_variant = dict()
    for variant in variant_freqs.index:
        y = variant_freqs.loc[variant]
        # Apply filter for smoothing
        y_smooth = y if not settings.smoothen_variantsLinePlot_with_filter else scipy.ndimage.gaussian_filter(y,sigma=settings.smoothen_variantsLinePlot_sigma)
        data_per_variant[variant] = y_smooth

        # Step 1: Loop through each x value (column) in data_per_variant
    for x in range(data_per_variant[next(iter(data_per_variant))].shape[0]):
        # Step 2: Calculate the sum of the smoothed values for all variants at this x value
        total_sum = sum(data_per_variant[variant][x] for variant in data_per_variant)
        if total_sum > 1:
            for variant in data_per_variant:
                data_per_variant[variant][x] /= total_sum

    # Create main stacked area chart
    fig = go.Figure()
    for variant in variant_freqs.index:
        x = variant_freqs.columns
        trace = go.Scatter(
            x=x,
            y=data_per_variant[variant],
            name=variant,
            mode="lines",
            stackgroup='stack',
            fill="tonexty",
            fillcolor=variants.globals.variantdict[variant].color,
            line_color=variants.globals.variantdict[variant].color
        )
        fig.add_trace(trace)

    # Add secondary y-axis for sample_ns data
    secondary_trace = go.Scatter(
        x=sample_n_values.index,
        y=sample_n_values,
        name="Sample N",
        mode="lines",
        line=dict(color="black", dash="dot"),
        yaxis="y2"  # Use secondary y-axis
    )
    fig.add_trace(secondary_trace)

    # Update layout for secondary axis
    fig.update_layout(
        width=4000,
        height=3000,
        title="SARS-CoV-2 Variant Frequencies Over Time",
        xaxis_title="Date",
        yaxis_title="Cumulative Frequency",
        legend_title="Variant",
        xaxis=dict(
            showgrid=True,  # Show vertical grid lines
            layer="above traces",  # Ensure grid is drawn above the traces
            gridwidth=1,  # Set grid line width
            gridcolor="LightGray",  # Color of vertical grid lines
            tickformat="%b\n%Y",  # Display month and year
            dtick="M1",  # Set tick every month
            tickangle=-45,  # Angle the ticks for better visibility
        ),
        yaxis=dict(
            showgrid=True,  # Show horizontal grid lines
            layer="above traces",  # Ensure grid is drawn above the traces
            gridwidth=1,  # Set grid line width
            gridcolor="LightGray",  # Color of horizontal grid lines
        ),
        yaxis2=dict(
            title="Sample N (log10)",
            overlaying="y",
            side="right",
            type="log",
            showgrid=False  # No grid on the secondary y-axis
        )
    )

    if settings.useDateAxis:
        fig.update_xaxes(
            tickformat="%b\n%Y",  # Display month and year
            dtick="M1",  # Set tick every month
            tickangle=-45,  # Angle the ticks for better visibility
            tickmode="linear",  # Use linear mode for consistent tick intervals
            ticks="outside",  # Place the ticks outside the axis
            tickson="boundaries",  # Ensure ticks appear at month boundaries
            showgrid=True  # Optionally, show vertical grid lines for better alignment
        )
    plotly_io.write_image(fig, settings.rootDir + "/VariantFreqLineChartTSV.pdf", format='pdf')