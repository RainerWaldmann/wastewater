from builtins import sorted

import pandas as pd
import plotly.graph_objects as go
import numpy as np
import math
import random
import copy
import settings
import Data
import variants as variants
from variants import OneVariant
import sys
from joblib import Parallel, delayed
import time



class VariantFrequencyData:
    """used to hold the data for each Variant and site"""
    def __init__(self, mean: float, se: float, n: int,  nPosAboveZero : int, boxplot: str, meanSetToZero: bool, meanWasCorrectedToFitParents : bool ):
        self.mean = mean
        self.se = se
        self.n = n # number of mutations
        self.nPosAboveZero = nPosAboveZero  # number of pos with above zero counts
        self.box = boxplot  # boxplot html data for individual mutations for this variant
        self.meanSetToZero = meanSetToZero  # True if mean was set to zero because not enough muts found
        self.meanWasCorrectedToFitParents = meanWasCorrectedToFitParents  # whether mean was corrected to fit parent

    def __repr__(self):
        return "<VariantData mean:%s se:%s>" % (self.mean, self.se)



def __getBoxPlot(x1, x2):
    """ used to get two boxplots for each variant frequency  : one for all data (x1) and one for the data where outliers are filtered out"""
    fig = go.Figure()
    fig.add_trace(go.Box(y=x1, boxpoints='all', name="Unfiltered (" + str(len(x1)) + ')', showlegend=False))
    if x2.empty == False:
        fig.add_trace(go.Box(y=x2, boxpoints='all', name="filtered (" + str(len(x2)) + ')', showlegend=False))
    fig.update_layout(height=400, width=200, hovermode=False)
    # fig.show()
    return fig


def __calcSpecialVariants (variantData):
    """treat special situations where mean is calculated from other values. Currently only supports A operator B which is defined in <calcstrategy> of tsv with variants"""
    for voc in variants.globals.variantdict:
        var = variants.globals.variantdict[voc]
        if var.calcstrategy:
            cd = var.calcstrategy.split(':')
            if len(cd) != 3:  # TODO test also whether both variants are in dict
                raise ValueError('calcstrategy split : should yield 3 elements')
            #op = operator.add if (Cd[1] == '+') else (lambda x,y: x-y if x>y else 0)
            op = lambda x,y: x-y if cd[1] == '-' else x+y
            #op = operator.add if (Cd[1] == '+') else operator.sub if Cd[1] == '-' else None
            if op:
                variantData.loc[var.name] = [
                    VariantFrequencyData(op(i.mean, j.mean) if j.mean <= i.mean else 0, math.sqrt(i.se ** 2 + j.se ** 2), i.n, i.nPosAboveZero, "",
                                         i.meanSetToZero, i.meanWasCorrectedToFitParents) if j != 1 else i for i, j in
                    zip(variantData.loc[cd[0]],
                        variantData.loc[cd[2]])]  # TODO protect against negative and other weird stuff


def __substractChildFreqsFromFakeParents(data : pd.DataFrame):
    """substract  child freqs from fake parents (e.g. BA4 BA 5 from BA2"""
    for voc in variants.globals.variantdict:
        nd = variants.globals.variantdict[voc]
        if nd.parentforcalc and nd.removeFromParentForCalc:
            data.loc[nd.parentforcalc.name] = [
                VariantFrequencyData(i.mean - j.mean if (i.mean - j.mean) > 0 else 0, math.sqrt(i.se ** 2 + j.se ** 2), i.n, i.nPosAboveZero,
                                     i.box, i.meanSetToZero, i.meanWasCorrectedToFitParents) for i, j in # TODO correct box
                zip(data.loc[nd.parentforcalc.name], data.loc[voc])]


def __correctChilds(varData : pd.DataFrame, node: OneVariant.OneVariant):
    """sets children to zero if parent zero, sets sum of childs <= parent. Will do it for given variant and all its childs recursively."""
    # lambda to calc sum of child freqs
    calculate_child_sum    = lambda s: [0] * len(varData.columns) if len(s) == 0 else [
        sum([float(x.mean) for x in o if not (math.isnan(x.mean))]) for o in
        zip(*[varData.loc[v.name] for v in s if not v.excempFromChildsSum])]
    calculate_correction_factor = lambda parent, child: 0 if parent.mean == 0 else parent.mean / child if child > parent.mean else 1
    if len(node.childs) > 0:
        childcorrFactor = [calculate_correction_factor(parent, child) for parent, child
                           in zip(varData.loc[node.name], calculate_child_sum(node.childs))]
        for c in node.childs:
            varData.loc[c.name] = [VariantFrequencyData(i.mean * j, i.se * j, i.n, i.nPosAboveZero, i.box, i.meanSetToZero, j != 1 and i.mean != 0) if j != 1 else i for i, j in zip(varData.loc[c.name], childcorrFactor)]
            __correctChilds(varData, c)


def __calcVariantFrequenciesOneVocThread(voc: str, dfp : pd.DataFrame, samples : list[str]):
        """function that generates the detailed counts for one variant. Will be called in parallel for each variant"""

        def getOutlierMask(data, mask):
            """returns mask where outliers are set to false -> filter values that should not be used for mean"""
            l = data[mask]
            mn = l.mean(skipna=True)
            qA, qB = np.percentile(l, [75, 25])
            intr_qr = qA - qB
            max = qA + (0.3 * intr_qr)
            min = qB - (0.3 * intr_qr)
            # protect from filtering values very close to mean if intr_qr is very small
            max = max if max > 1.05 * mn else 1.05 * mn
            min = min if min < 0.95 * mn else 0.95 * mn
            return data.between(min, max)  # mask with outliers


        # get lines that match muts of variant
        oneVOCdetailedCounts = dfp[dfp.index.isin([st for st in dfp.index if
                                                             any(sub.getNAmutstring().lower() in st.lower() for sub
                                                                 in
                                                                 variants.globals.variantdict[voc].data)])]
        # add lines for variant mutations that were not found in merged ivar
        df_variant_Indexaslist = [x.split('_')[0] for x in oneVOCdetailedCounts.index]
        # df_variant.index.get_level_values(0)
        resort = False
        for sub in variants.globals.variantdict[voc].data:
            if sub.getNAmutstring() not in df_variant_Indexaslist:
                resort = True
                # df_variant.loc[sub.getNAmutstring() + sub.getAAmutstringForDFindex()] = \
                #     [0 if sample.getCoverageAtPos(sub.position) is not None and sample.getCoverageAtPos(
                #         sub.position) > settings.minDepth else float('nan') for sample in sampledata_list]
                oneVOCdetailedCounts = dfp[dfp.index.isin([st for st in dfp.index if
                                                                     any(sub.getNAmutstring().lower() in st.lower()
                                                                         for
                                                                         sub in
                                                                         variants.globals.variantdict[
                                                                             voc].data)])].copy()
        if resort:
            oneVOCdetailedCounts.sort_index(key=lambda x: x.str.extract(r"^[^\d]*(\d+)", expand=False).astype(int),
                                            inplace=True)

        # TODO following line only checks whether it is present in index not whether it is non 0
        allrequiredMutsFound = all([any([y in x for x in oneVOCdetailedCounts.index]) for y in
                                    variants.globals.variantdict[voc].getRequiredMutations()])
        # if not enough mutations found in dataset, bail out, return None
        if oneVOCdetailedCounts.shape[0] < variants.globals.variantdict[
            voc].minmutsforpass or not allrequiredMutsFound:
            return voc, None, None, [VariantFrequencyData(0, np.nan, 0, 0, "", True, False) for sample in samples]

        oneVocdetailedCountsMask = pd.DataFrame(columns=oneVOCdetailedCounts.columns,
                                                index=oneVOCdetailedCounts.index)  # just create an empty dataframe with the right dimension here. Will be populated later

        variant_data_dictOneVOC = {}  # holds the VariantData objects for one voc for each sample
        for sample in samples:
            # print("\t" + sample)
            x = oneVOCdetailedCounts[sample]  # Var freqs for one variant for one sample
            mask = [a == a for a in x]  # nans are marked false
            nPosAboveZero = (x[mask] > 0).sum()
            nPositions = len(x)
            hasRequiredMuts = True
            # variants.globals.variantdict[voc].testWhetherBothPartsOfChimeraOK(x)
            minMutsAboveZero = True
            if variants.globals.variantdict[voc].minmutsforpass is not None and variants.globals.variantdict[
                voc].minmutsforpass > 0 and (x[mask]).sum() > 0:  # don't set flag if all are 0
                minMutsAboveZero = (x[mask] != 0).sum() >= variants.globals.variantdict[voc].minmutsforpass
            if variants.globals.variantdict[voc].minstarredmutsforpass > 0 and (x[mask]).sum() > 0:
                # requiredMuts = [b for b in x.index if any([a in b for a in variants.globals.variantdict[voc].getRequiredMutations()])]
                requiredMutsFound = x[mask][[b for b in x[mask].index if any([a in b for a in
                                                                              variants.globals.variantdict[
                                                                                  voc].getRequiredMutations()])]]
                hasRequiredMuts = len(requiredMutsFound) >= variants.globals.variantdict[
                    voc].minstarredmutsforpass and all([c > 0 for c in
                                                        requiredMutsFound])  # Todo check this since it requires all starred to be above 0

            if sum(mask) >= 3:
                # print(x[mask])
                mask = np.logical_and(getOutlierMask(x, mask), mask)
                box = __getBoxPlot(x, x[mask]).to_html(full_html=False, include_plotlyjs='cdn') if \
                    settings.plotBoxPlotsForDetailedVarCounts else ""
            elif sum(mask) > 1:
                    box = __getBoxPlot(x, x[mask]).to_html(full_html=False, include_plotlyjs='cdn') if \
                        settings.plotBoxPlotsForDetailedVarCounts else ""
            else:
                box = ""

            x = x[mask]
            oneVocdetailedCountsMask[sample] = mask
            if sum(mask) > 0 and hasRequiredMuts and minMutsAboveZero:  # at least one left
                variant_data_dictOneVOC[sample] = VariantFrequencyData(x.mean(skipna=True),
                                                                       np.nan if sum(mask) <= 1 else x.std(
                                                                  skipna=True) / math.sqrt(sum(mask)),
                                                                       nPositions, nPosAboveZero, box, False,
                                                                       False)
            else:
                variant_data_dictOneVOC[sample] = VariantFrequencyData(0, np.nan, nPositions, nPosAboveZero, "",
                                                                       True, False)

        return voc, oneVOCdetailedCounts, oneVocdetailedCountsMask, variant_data_dictOneVOC


def calcVariantFrequencies(sample_list : list[Data.OneBCdata] ,dfp : pd.DataFrame):
    """Calculates the variant frequencies for each variant and sample.
    Will also generate the detailed counts (frequencies for each mutation) for each variant and sample and a similar df with a mask indicating which mutations were used to calculate mean. """

    start_time = time.time()
    samples = dfp.columns
    # empty dataframe for mean and SD, values in object of type VariantData that holds mean, SD, boxplots .....
    variantDataMeansSD = pd.DataFrame(columns=[s for s in samples],
                                      index=[x for x in  # dfp.colums -> series of samples
                                             variants.globals.variantdict])  # colums are samples, rows are variants
    variantDataMeansSDUncorrected = None  # will hold copy of variantDataMeansSD before corrections
    detailedCounts = dict()  # holds data for each variant (key is variant)
    detailedCountsMask = dict()  # Mask with TRUE false that indicate whether value was used for mean (True)  or is an outlier (False) , key of dict is variant name

    #with concurrent.futures.ProcessPoolExecutor() as executor:
    #    executor.map(generateDetailedCounts, variants.globals.variantdict)
    vd = Parallel(n_jobs=settings.num_cores)(
        delayed(__calcVariantFrequenciesOneVocThread)(voc, dfp, samples) for voc in variants.globals.variantdict)
    print(f"DefineLineagesThreads finished {time.time() - start_time} seconds to complete.")
    for voc, df_variant, df_variant_mask, onevocvariantDataMeansSD in vd:
        variantDataMeansSD.loc[voc] = onevocvariantDataMeansSD
        if df_variant is not None: # if not enough mutations found, df_variant will be None
            detailedCounts[voc]= df_variant
            detailedCountsMask[voc] = df_variant_mask
            variantDataMeansSD.loc[voc] = pd.Series(onevocvariantDataMeansSD)

    # ******* CORRECTION SPIKE-IN ****************

    # correct sum of all variants to close to 100% if spike-in was used  - TEMPORARY FIX , FIND permanent solution later
    # calculate correction factor for each sample to set Omicron + Delta to a random number between 0.985 and 1
    # find head nodes for variants
    headnodes = [variants.globals.variantdict[v] for v in variants.globals.variantdict if (
                variants.globals.variantdict[v].isHead() and variants.globals.variantdict[v].parentforcalc is None)]
    headnodesdata = [variantDataMeansSD.loc[x.name] for x in headnodes]
    # correction factor to set sum of all root variants = 1
    rTotal = [1 / (sum([float(x.mean) for x in o if not (math.isnan(x.mean))])) if sum(
            [float(x.mean) for x in o if not (math.isnan(x.mean))]) != 0 else 1 for o in zip(*headnodesdata)]
    rTotal = pd.Series(rTotal, index=samples)
    rTotal = [rTotal[s.sample] if s.trimstatsSpikeIn != None else 1 for s in
              sample_list]  # set rTotal to 1 when NO spike-in was used = don't correct

    def __correctionSpikeIn(x):
        """Will correct individual values, if correction yiels value > 1, will set value to random between 0.990 and 1"""
        for i in range(x.size):
            if x.iloc[i] == x.iloc[i]:
                x.iloc[i] = x.iloc[i] if rTotal[i] == 1 else rTotal[i] * x.iloc[i] if rTotal[i] * x.iloc[
                    i] <= 1 else random.randint(990, 1000) / 1000

    for voc in variants.globals.variantdict:
        # print("\nBefore correction: " + voc + "\n", detailedCounts[voc])
        if voc in detailedCounts:
            detailedCounts[voc].apply(lambda x: __correctionSpikeIn(x),
                                      axis=1)  # np.asarray(x) * np.asarray(rTotal), axis=1)
            # set means  to the corrected values
            variantDataMeansSD.loc[voc] = [
                VariantFrequencyData(i.mean * j, i.se * j, i.n, i.nPosAboveZero, i.box, i.meanSetToZero,
                                     i.meanWasCorrectedToFitParents) for i, j in
                zip(variantDataMeansSD.loc[voc], rTotal)]

    # ***************   SPIKE - IN correction end *************************************************

    # substract  child freqs from fake parents (e.g. BA4 BA 5 from BA2
    __substractChildFreqsFromFakeParents(variantDataMeansSD)
    __calcSpecialVariants(variantDataMeansSD)
    # keep copy of uncorrected data before corrections for PieChart
    variantDataMeansSDUncorrected = copy.deepcopy(variantDataMeansSD)

    # correct childs
    # avoid sum of sublineages to depass parents, typically a problem for very low frequency variants
    for headnodesdata in headnodes:
        __correctChilds(variantDataMeansSD, headnodesdata)
    if settings.showprocesstime:
        print(f"-----------------calcVariantFrequencie finished took {time.time() - start_time} seconds to complete.-------------------------------")
    return detailedCounts, detailedCountsMask, variantDataMeansSD, variantDataMeansSDUncorrected



def printVariantsTSV(variantData : pd.DataFrame):
    """Save TSV with variant frequencies to file."""
    # # Step 1: Transform pieData to contain only the required attributes
    # df = variantData.applymap(lambda x: [x.mean, x.se, x.n])
    #
    # # Step 2: Split the list in each cell into multiple columns
    # df = pd.concat([df[col].apply(pd.Series) for col in df], axis=1)
    #
    # # Step 3: Rename the columns
    # df.columns = pd.MultiIndex.from_product([variantData.columns, ['mean', 'se', 'n']])
    #
    # # Step 4: Write df to an Excel file
    # output_path = os.path.join(settings.rootDir, settings.frequenciesOutFile)
    # print(f"Output file path: {output_path}")
    # df.to_excel(output_path)
    # pass
    # Step 1: Transform pieData to contain only the required attributes
    df = variantData.applymap(lambda x: x.mean)
    df = df.transpose()
    # Step 2: Reshape df from wide format to long format
    df = df.unstack().reset_index()
    df.columns = ['variable', 'attribute', 'value']

    # Step 3: Reshape df back to wide format with multi-level columns
    df = df.pivot(index='variable', columns='attribute', values='value')

    # Step 4: Write df to an Excel file
    output_path = settings.rootDir + '\\' + settings.frequenciesOutFile
    print(f"Output file path: {output_path}")
    df.to_excel(output_path)
    #means = variantDataMeansSD.applymap(lambda x: x.mean)
    #means.to_csv(settings.rootDir + "/" + settings.frequenciesOutFile, sep="\t")
     # writer = pd.ExcelWriter(settings.rootDir + "/" + settings.frequenciesOutFile, engine='xlsxwriter')
     # means = variantData.apply(lambda x : x.mean)
     # means.to_excel(writer, sheet_name='means', index=False)
     # se = variantData.apply(lambda x: x.se)
     # se.to_excel(writer, sheet_name='se', index=False)
     # writer.save()

