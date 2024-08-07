from builtins import sorted

import pandas as pd
import plotly.graph_objects as go
import plotly.io as plotly_io
import numpy as np
import math
import random
import copy
import settings
import Data
import variants as variants
import variants.globals as vg
#import variants.Node as nd
import operator
import sys
from itertools import groupby
import re

#calcSEsumOrDiff = lambda x,y : math.sqrt(x**2 + y**2) # calc se for sum or difference

class DataForHTML:
    """Just to group various data used in html. Avoids using too many seperate parameters for Jinja"""
    def __init__(self, variantPieList,variantsDetailedCounts, variantsdetailedCountsMask, variantDataMeans,variantHistoFigHTML):
        self.variantPieList = variantPieList # the list of pie charts
        self.variantsDetailedCounts = variantsDetailedCounts # the Var frequencies for each mutation for each variant for the different samples
        self.variantsdetailedCountsMask = variantsdetailedCountsMask # same structure as variantsDetailedCounts, contains just true or false to indicate whether value is an outlier
        self.variantDataMeans  = variantDataMeans # contains VariantData object for each variant and sample -> can retrieve means and se
        self. variantHistoFigHTML = variantHistoFigHTML # boxplot html for each variant and site showing the median .. for all mutation defining the variant


class VariantData:
    """used to hold the data for each Variant and site"""
    def __init__(self, mean: float, se: float, n: int,  nPosAboveZero : int, boxplot: str, meanSetToZero : bool , meanWasCorrectedToFitParents : bool ):
        self.mean = mean
        self.se = se
        self.n = n # number of mutations
        self.nPosAboveZero = nPosAboveZero # number of pos with above zero counts
        self.box = boxplot # boxplot html data for individual mutations for this variant
        self.meanSetToZero = meanSetToZero # True if mean was set to zero because not enough muts found
        self.meanWasCorrectedToFitParents = meanWasCorrectedToFitParents # whether mean was corrected to fit parent

    def __repr__(self):
        return "<VariantData mean:%s se:%s>" % (self.mean, self.se)

def getBoxPlot(x1, x2):
    """ used to get two boxplots for each variant and site : one for all data (x1) and one for the data where outliers are filtered out"""
    fig = go.Figure()
    fig.add_trace(go.Box(y=x1, boxpoints='all', name="Unfiltered (" + str(len(x1)) + ')', showlegend=False))
    if x2.empty == False:
        fig.add_trace(go.Box(y=x2, boxpoints='all', name="filtered (" + str(len(x2)) + ')', showlegend=False))
    fig.update_layout(height=400, width=200, hovermode=False)
    # fig.show()
    return fig


def calcSpecialVariants (variantData):
    """treat special situations where mean is calculated from other values. Currently only supports A operator B which is defined in <calcstrategy> of tsv with variants"""
    for voc in variants.globals.variantdict:
        var = variants.globals.variantdict[voc]
        if var.calcstrategy:
            cd = var.calcstrategy.split(':')
            if len(cd) != 3:  # TODO test also whether both variants are in dict
                sys.exit('calcstrategy split : should yield 3 elements')
            #op = operator.add if (Cd[1] == '+') else (lambda x,y: x-y if x>y else 0)
            op = lambda x,y: x-y if cd[1] == '-' else x+y
            #op = operator.add if (Cd[1] == '+') else operator.sub if Cd[1] == '-' else None
            if op:
                variantData.loc[var.name] = [
                    VariantData(op(i.mean, j.mean) if j.mean <= i.mean else 0, math.sqrt(i.se ** 2 + j.se ** 2), i.n, i.nPosAboveZero,"",
                                i.meanSetToZero, i.meanWasCorrectedToFitParents) if j != 1 else i for i, j in
                    zip(variantData.loc[cd[0]],
                        variantData.loc[cd[2]])]  # TODO protect against negative and other weird stuff

def substractChildFreqsFromFakeParents(data : pd.DataFrame):
    """substract  child freqs from fake parents (e.g. BA4 BA 5 from BA2"""
    for voc in variants.globals.variantdict:
        nd = variants.globals.variantdict[voc]
        if nd.parentforcalc and nd.removeFromParentForCalc:
            data.loc[nd.parentforcalc.name] = [
                VariantData(i.mean - j.mean if (i.mean - j.mean) > 0 else 0, math.sqrt(i.se ** 2 + j.se ** 2), i.n,i.nPosAboveZero,
                            i.box, i.meanSetToZero, i.meanWasCorrectedToFitParents) for i, j in # TODO correct box
                zip(data.loc[nd.parentforcalc.name], data.loc[voc])]



def plotVariants(sample_list : list[Data.OneBCdata] ,dfp : pd.DataFrame): #sample_list  only required to get info whether spike-in was used -> tests whether spikein data were used to correct to 100% in case of
    """The function called from main"""
    def getOutlierMask(l):
        """returns mask where outliers are set to false -> filter values that should not be used for mean"""
        mn = l.mean(skipna=True)
        qA, qB = np.percentile(l, [75, 25])
        intr_qr = qA - qB
        max = qA + (0.4 * intr_qr)
        min = qB - (0.4 * intr_qr)
        # protect from filtering values very close to mean if intr_qr is very small
        max = max if max > 1.05 * mn else 1.05 * mn
        min = min if min < 0.95 * mn else 0.95 * mn
        return x.between(min, max)  # mask with outliers

    samples = dfp.columns # just to get the sample list
    detailedCounts = dict() # holds data for each variant (key is variant)
    # Mask with TRUE false that indicate whether value was used for mean (True)  or is an outlier (False) , key of dict is variant name
    detailedCountsMask = dict() # key is variant name
    fold_sd_dev_from_mean = 1.5  # this fold sd is cutoff for deviation from mean

    # create dataframe for each variant that contains the mut frequencies for the individual mutations for this variant for each site
    for voc in variants.globals.variantdict:
        #get lines that match muts of variant
        df_variant = dfp[dfp.index.isin([st for st in dfp.index if
                                         any(sub.getNAmutstring().lower() in st.lower() for sub in variants.globals.variantdict[voc].data)])]
        # add lines for variant mutations that were not found in merged ivar
        df_variant_Indexaslist = [x.split('_')[0] for x in df_variant.index]
        # df_variant.index.get_level_values(0)
        resort = False
        for sub in variants.globals.variantdict[voc].data:
            if sub.getNAmutstring() not in df_variant_Indexaslist:
                resort = True
                df_variant.loc[sub.getNAmutstring() + sub.getAAmutstringForDFindex()] = \
                    [0 if sample.getCoverageAtPos(sub.position) is not None and sample.getCoverageAtPos(
                        sub.position) > settings.minDepth else float('nan') for sample in sample_list]
        if resort:
            df_variant.sort_index(key=lambda x: x.str.extract(r"^[^\d]*(\d+)", expand=False).astype(int), inplace=True)

        #TODO following line only checks whether it is present in index not whether it is non 0
        allrequiredMutsFound = all([any([y in x for x in df_variant.index]) for y in variants.globals.variantdict[voc].getRequiredMutations()])
        if df_variant.shape[0] >= variants.globals.variantdict[voc].minmutsforpass and allrequiredMutsFound: # if too few mutations found in dataset, ignore it COUNTS --> NO COUNTS CHECKED HERE
            detailedCounts[voc] = df_variant
            detailedCountsMask[voc] = pd.DataFrame(columns=df_variant.columns, index=df_variant.index) # just create an empty dataframe with the right dimension here. Will be populated later

    # empty dataframe for mean and SD, values in object of type VariantData that holds mean, SD, boxplots .....
    variantDataMeansSD = pd.DataFrame(columns=[s for s in samples], index=[x for x in
                                                                           variants.globals.variantdict])  # colums are samples, rows are variants
    # populate variantDataMeansSD and filter outliers

    for voc in variants.globals.variantdict:
        if voc not in detailedCounts:  # if no data generate 0 mean data
            for sample in samples:
                variantDataMeansSD[sample].loc[voc] = VariantData(0, np.nan, 0, 0,"", True,False)
        else:
            #print(voc)
            for sample in samples:
                #print("\t" + sample)
                x = detailedCounts[voc][sample] # Var freqs for one variant for one sample
                mask = [a == a for a in x]  # nans are marked false
                nPosAboveZero = (x[mask]>0).sum()
                nPositions = len(x)
                hasRequiredMuts = True
                #variants.globals.variantdict[voc].testWhetherBothPartsOfChimeraOK(x)
                minMutsAboveZero = True
                if variants.globals.variantdict[voc].minmutsforpass is not None and variants.globals.variantdict[voc].minmutsforpass > 0 and (x[mask]).sum()>0: # don't set flag if all are 0
                    minMutsAboveZero = (x[mask]!=0).sum() >= variants.globals.variantdict[voc].minmutsforpass
                if variants.globals.variantdict[voc].minstarredmutsforpass > 0 and (x[mask]).sum()>0:
                    #requiredMuts = [b for b in x.index if any([a in b for a in variants.globals.variantdict[voc].getRequiredMutations()])]
                    requiredMutsFound = x[mask][[b for b in x[mask].index if any ([a in b for a in variants.globals.variantdict[voc].getRequiredMutations()])]]
                    hasRequiredMuts = len(requiredMutsFound) >= variants.globals.variantdict[voc].minstarredmutsforpass and all([c > 0 for c in requiredMutsFound]) # Todo check this since it requires all starred to be above 0

                if sum(mask) >= 3: # mask outliers if at least 3 values
                        #print(x[mask])
                        mask = np.logical_and(getOutlierMask(x[mask]), mask)

                if sum(mask) > 1:
                        box = getBoxPlot(x, x[mask]).to_html(full_html=False, include_plotlyjs='cdn') if \
                            settings.plotBoxPlotsForDetailedVarCounts else ""
                else:
                        box = ""

                x = x[mask]
                detailedCountsMask[voc][sample] = mask
                if sum(mask) > 0 and hasRequiredMuts and minMutsAboveZero: #at least one left
                    variantDataMeansSD[sample].loc[voc] = VariantData(x.mean(skipna=True),
                                                       np.nan if sum(mask) <= 1 else x.std(skipna=True) / math.sqrt(sum(mask)),
                                                       nPositions, nPosAboveZero , box, False,False)
                else:
                    variantDataMeansSD[sample].loc[voc] = VariantData(0, np.nan, nPositions,nPosAboveZero , "", True,False)
                dummy = 0 # just for debugging stop

    def __correctChilds(node: variants.Node, correctionlambda):
        """sets children to zero if parent zero, sets sum of childs <= parent"""
        #dummy = sumchilds(node.childs)
        if len(node.childs) > 0:
            childcorrFactor = [correctionlambda(parent, child) for parent, child
                               in zip(variantDataMeansSD.loc[node.name], sumchilds(node.childs))]
            for c in node.childs:
                variantDataMeansSD.loc[c.name] = [
                    VariantData(i.mean * j, i.se * j, i.n,i.nPosAboveZero, i.box, i.meanSetToZero, j != 1 and i.mean !=0) if j != 1 else i for i, j in
                    zip(variantDataMeansSD.loc[c.name], childcorrFactor)]
                __correctChilds(c, correctionlambda)

    # ******* CORRECTION SPIKE-IN ****************

        # correct sum of all variants to close to 100% if spike-in was used  - TEMPORARY FIX , FIND permanent solution later
        # calculate correction factor for each sample to set Omicron + Delta to a random number between 0.985 and 1
        # find head nodes for variants
    headnodes = [variants.globals.variantdict[v] for v in variants.globals.variantdict if (variants.globals.variantdict[v].isHead() and variants.globals.variantdict[v].parentforcalc is None)]
    h = [variantDataMeansSD.loc[x.name] for x in headnodes]
    #correction factor to set sum of all root variants = 1
    rTotal = [ random.randint(985, 1000) / (1000 * sum([float(x.mean) for x in o if not (math.isnan(x.mean))])) if  sum([float(x.mean) for x in o if not (math.isnan(x.mean))]) !=0 else 1 for o in zip(*h)]
    rTotal = pd.Series(rTotal, index=samples)
    rTotal = [ rTotal[s.sample] if s.trimstatsSpikeIn != None else 1 for s in sample_list] #set rTotal to 1 when NO spike-in was used = don't correct


    def correctionSpikeIn(x):
        """Will correct individual values, if correction yiels value > 1, will set value to random between 0.990 and 1"""
        for i in range(x.size):
            if x.iloc[i] == x.iloc[i]:
                x.iloc[i] = x.iloc[i] if rTotal[i] == 1 else rTotal[i] * x.iloc[i] if rTotal[i] * x.iloc[i] <= 1 else random.randint(990, 1000) / 1000

    for voc in variants.globals.variantdict:
        #print("\nBefore correction: " + voc + "\n", detailedCounts[voc])
        if voc in detailedCounts:
            detailedCounts[voc].apply(lambda x: correctionSpikeIn(x),axis=1) #np.asarray(x) * np.asarray(rTotal), axis=1)
            # set means  to the corrected values
            variantDataMeansSD.loc[voc] = [VariantData(i.mean * j, i.se * j, i.n, i.nPosAboveZero, i.box, i.meanSetToZero,i.meanWasCorrectedToFitParents)  for i, j in
                                             zip(variantDataMeansSD.loc[voc], rTotal)]

    # ***************   SPIKE - IN correction end *************************************************

    #substract  child freqs from fake parents (e.g. BA4 BA 5 from BA2
    substractChildFreqsFromFakeParents(variantDataMeansSD)
    # keep copy of uncorrected data before corrections for PieChart
    variantDataMeansSDUncorrected = copy.deepcopy(variantDataMeansSD)
    calcSpecialVariants(variantDataMeansSDUncorrected)

    calcSpecialVariants(variantDataMeansSD)

    # correct childs
    # lambda to calc sum of child freqs
    sumchilds = lambda s: [0] * len(variantDataMeansSD.columns) if len(s) == 0 else [
        sum([float(x.mean) for x in o if not (math.isnan(x.mean))]) for o in
        zip(*[variantDataMeansSD.loc[v.name] for v in s if not v.excempFromChildsSum])]



    # avoid sublineages to depass parents,
    childcorrectionLambda = lambda parent,child: 0 if parent.mean == 0 else parent.mean / child if child > parent.mean else 1

    for h in headnodes:
        __correctChilds(h,childcorrectionLambda)

    #make sublineages sum equal to  parents if sum greater than parent,
    #childcorrectionLambda = lambda parent,child: 0 if parent.mean == 0 else 1 if child <= parent.mean  else parent.mean / child
    #for h in headnodes:
    #    __correctChilds(h,childcorrectionLambda)


    printVariantsTSV(variantDataMeansSD)
    # ****************    PIECHARTS ****************************************************************
    variantsToPrint = [v for v in variants.globals.variantdict if variants.globals.variantdict[v].doPrint == True]
    variantPieList = None

    pieData = variantDataMeansSD.loc[variantsToPrint]

    if settings.plotVariantPies:
        parents = [variants.globals.variantdict[v].parent.name if variants.globals.variantdict[v].parent else "" for v in variantsToPrint]
        variantcolors = [variants.globals.variantdict[v].color for v in variantsToPrint]
        variantPieList = dict()
        for i in range(0, len(samples)):
            #    row = math.ceil((i + 1) / nCols)
            #    col = i % nCols + 1

            hoverdata = [[v.mean, v.se if v.se == v.se else "NA", v.n, v.nPosAboveZero] for v in pieData[samples[i]]]
            trace = go.Sunburst(labels=variantsToPrint, parents=parents,
                                values=[v.mean for v in pieData[samples[i]]],
                                branchvalues="total",
                                customdata=hoverdata,
                                hovertemplate='<b>Variant: %{label}<b><br>Frequency: %{value}<br>SE: %{customdata[1]:.3f}<br> %{customdata[3]} mutations of %{customdata[2]} above zero')
            figForList = go.Figure(trace)
            figForList.update_traces(
                marker=dict(colors=variantcolors, line=dict(color='#000000', width=0.5)))
            figForList.update_layout(height=250, width=250,
                                     margin=dict(t=10, l=0, r=0,
                                                 b=10))  # ,title=samples[i],margin = dict(t=0, l=0, r=0, b=0))
            variantPieList[samples[i]] = figForList.to_html(full_html=False, include_plotlyjs='cdn')
            #figForList.show()


    #######   HISTOGRAM ####################################################################
    usedVariantDataMeansSD = variantDataMeansSD if settings.variantHistogram_sumChildsNotExceedParentCorrection else variantDataMeansSDUncorrected

    variantDataTr = usedVariantDataMeansSD.transform(lambda x: [v.mean for v in x])
    variantDataErr = usedVariantDataMeansSD.transform(lambda x: [v.se for v in x])
    #list of variants for histogram, only variants where print is enabled and where any of values is > 0
    variantsForHisto = [v for v in variants.globals.variantdict if variants.globals.variantdict[v].histogramGroup and variantDataTr.loc[v].gt(0).any() and variants.globals.variantdict[v].doPrint]
    # group variants by histogram bar, list must be sorted for groupby to work --> sort by variants.globals.variantdict[x].histogramOrderId
    groupedVariants = {key: list((value)) for key, value in groupby(sorted(variantsForHisto, key=lambda x: variants.globals.variantdict[x].histogramOrderId), lambda x: variants.globals.variantdict[x].histogramOrderId)}
    data = [] # the parts of the hustogram
    for offsetGroupIndex, oneGroup in groupedVariants.items():
        def defineOffset(variant: variants.Node):  # get all childs and define offsets
            """function that calculates the offset of childs within a bar, baseline if offset of parent"""
            childs = [v for v in oneGroup if variants.globals.variantdict[v].parent.name is variant]
            if childs:
                childOffset = oneGroupWithOffsets[variant] # first child offset equal parent offset
                for onechild in childs:
                    oneGroupWithOffsets[onechild] = childOffset
                    #groupWithFinishedStatus[onechild] = True
                    childOffset = (childOffset + variantDataTr.loc[onechild])#.copy() # increase child offset for next round by frequency of current child
                    defineOffset(onechild) # treat childs of current child

        oneGroup = sorted(oneGroup, key=lambda x: variants.globals.variantdict[x].histogramOrderIdWithinBar)
        oneGroupWithOffsets = {v:0 for v in oneGroup} # holds variant name and offset of histogram bar

        #get top level variants (variants that don't have parent in current group)
        topNodes = [v for v in oneGroup if
                    (variants.globals.variantdict[v].parent is None or variants.globals.variantdict[v].parent.name not in oneGroup)]
        if(len(topNodes)==0): # normally should not happen
            continue
        #set offset of first one to 0
        oneGroupWithOffsets[topNodes[0]] = 0;
        for kk in range(1,len(topNodes)):
            oneGroupWithOffsets[topNodes[kk]] = (oneGroupWithOffsets[topNodes[kk-1]] + variantDataTr.loc[topNodes[kk-1]])
        for variant in topNodes:
           #groupWithFinishedStatus[variant] = True
           defineOffset(variant)
        #legendgrouptitletext = "" if offsetGroupName.isdigit() else offsetGroupName

        for oneVariant in oneGroup:
                offsetGroupName = variants.globals.variantdict[oneVariant].histogramGroup
                data.append(go.Bar(
                    #width = 0.5,
                   # opacity=0.5,
                    marker_line_color='black',
                    legendgroup= variants.globals.variantdict[oneVariant].histogramGroup,
                    legendgrouptitle_text = offsetGroupName,
                    name=oneVariant,
                    x=variantDataTr.columns,
                    y=variantDataTr.loc[oneVariant],
                    offsetgroup=offsetGroupIndex,
                    base=oneGroupWithOffsets[oneVariant], #
                    marker_color=variants.globals.variantdict[oneVariant].color,
                    error_y=dict(type='data', array=variantDataErr.loc[oneVariant]),
                    #error_y_color=VariantDefinitions.VariantErrorbarColors[4]
                    #text=oneVariant,
                    #textposition='none',
                    meta=[oneVariant],
                    hovertemplate = "<br>".join([
                        "Variant: %{meta[0]}",
                        "Frequency: %{y}"])
                ))

    figVariantHisto = go.Figure(data=data)
    figVariantHisto.update_layout(height=1000, width=len(variantDataTr.columns) * 100 + 300, paper_bgcolor='rgba(0,0,0,0)',
                                  plot_bgcolor='rgba(0,0,0,0)', barmode='group', bargap=0.4,bargroupgap=0.3,
                                  )
                                  #legendgrouptitle_font_size = 20,
                                  #legendgrouptitle_font_color="green", legendgrouptitle_font_family="Times New Roman")
    if settings.doPlot :
        figVariantHisto.show()

     ##fig = px.bar(variantDataTr, x=variantDataTr.index, y=variantDataTr.columns, title="Variants",text_auto='.2f')
    # fig.update_traces(texttemplate=’%{percent:.3f}’)
    # fig.update_traces(marker=dict(line=dict(color='#000000', width=0.5)))
    # fig.update_layout(height=1500, width = 1200)
    plotly_io.write_image(figVariantHisto, settings.rootDir + "/Varianthisto.pdf", format='pdf')
    #add color dots for HTML to detailed counts
    for voc in detailedCounts:
        mutinfoList = [[d for d in variants.globals.variantdict[voc].data if d.getNAmutstring() in x][0] for x in
                   detailedCounts[voc].index]  # get corresponding MutIfo Objects for the mutation
        colorDotList = ["".join([x.getHTMLforColoredDot() for x in
                             l.listOfOtherVariantsMutIsFound]) if l.listOfOtherVariantsMutIsFound is not None else ""
                    for l in mutinfoList]  # get the dots to print after mut
        detailedCounts[voc]["dots"] = colorDotList

        # Add boxplots for VOCs
    variantsToPrint = [v for v in variants.globals.variantdict if
                       variants.globals.variantdict[v].showVOCplot == True]
    for voc in variantsToPrint:
        data = detailedCounts[voc]

    return DataForHTML(variantPieList,detailedCounts,detailedCountsMask, variantDataMeansSD, figVariantHisto.to_html(full_html=False, include_plotlyjs='cdn'))



def printVariantsTSV(variantData):
    """TODO: NOT CODED YET"""
    means = variantData.applymap(lambda x: x.mean)
    means.to_csv(settings.rootDir + "/" + settings.frequenciesOutFile, sep="\t")
     # writer = pd.ExcelWriter(settings.rootDir + "/" + settings.frequenciesOutFile, engine='xlsxwriter')
     # means = variantData.apply(lambda x : x.mean)
     # means.to_excel(writer, sheet_name='means', index=False)
     # se = variantData.apply(lambda x: x.se)
     # se.to_excel(writer, sheet_name='se', index=False)
     # writer.save()

