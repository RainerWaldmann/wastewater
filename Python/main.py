import Globals
import pandas as pd
from jinja2 import Template
import math

import DefineLineages
import PlotVariantFrequencies
import ReadData
import settings
import HeatMap
import Coverage
import Pepper
import JinjaTemplate
import argparse
import MutlLstProcessing.ProcessMutlist
import variants.globals
import seqtools.CodonTables as Cd
import time
import sys
from concurrent.futures import ThreadPoolExecutor

start_time = time.time()

pd.options.display.float_format = '{:,.3f}'.format

import variants.MutInfo

if settings.heatmapDisplayOnlyVariants:
    settings.hoverinfo_Onheatmaps = False


def analyzeSamplesMain():
    sample_list = ReadData.readData()
    functions = [ReadData.mergeAndFilterData, Coverage.plotCoverage, Coverage.plotCoverageViolinPlots, Coverage.doShannon, Coverage.getShannonBoxPlots]
    start_time = time.time()

    with ThreadPoolExecutor() as executor:
        mergedData, coveragePlots_fig, violinCoveragePlots_fig, shannonPlots_fig, shannonBoxFig = list(executor.map(lambda func: func(sample_list), functions))

    if settings.showprocesstime:
        print(f"----------------- mergedData, coveragePlots_fig, violinCoveragePlots_fig, shannonBoxFig, shannonPlots_fig finished took {time.time() - start_time} seconds to complete.-------------------------------")
    #write TSV with merged data
    mergedData.to_csv(settings.rootDir + "/VarFrequencies.tsv", sep="\t", index=False)
    # get dataframes with just frequencies
    frequencyDF = ReadData.getFrequencyColumnsOnly(mergedData, sample_list, None)
    frequencyDF.to_csv(settings.rootDir + "/VarFrequenciesLight.tsv", sep="\t")
    settings.plotPepper = settings.plotPepper & all([s.ppmov is not None for s in sample_list])
    if settings.plotPepper and any([s.ppmov is None or s.articPool1Count is None for s in sample_list]): #if data for ppmov counting not present, disable it
        settings.plotPepper = False
        Globals.warningList.put("settings.plotPepper set to False because no artic pool 1 count found in one or more of the .trimstats files")
    if settings.plotPepper:
        pepperHistogram_fig, pepperDF = Pepper.countPepper(sample_list)
    else:
        pepperHistogram_fig = None
        pepperDF = None

    mergedDataForHeatMaps = ReadData.getMergedDataForHeatMap(mergedData, sample_list)
    if mergedDataForHeatMaps.empty:
        print("Error: No mutation passed filter for heatmaps(min reads). ---- EXITING")
        sys.exit(1)
    start_time = time.time()
    functions = [HeatMap.plotHeatMap, HeatMap.plotHeatMapVariantsOnly, HeatMap.clusterMap]

    with ThreadPoolExecutor() as executor:
        heatmap_fig, heatmapVariantsOnly, clusterMap_fig = list(executor.map(lambda func: func(mergedDataForHeatMaps, sample_list), functions))
        
    if settings.showprocesstime:
        print(f"----------------- heatmap_fig, heatmapVariantsOnly, clusterMap_fig finished took {time.time() - start_time} seconds to complete.-------------------------------")

    #variant_dataForHtml = Var.plotVariants(sample_list, frequencyDF)
    detailedCounts, detailedCountsMask, variantDataMeansSD, variantDataMeansSDUncorrected = DefineLineages.calcVariantFrequencies(sample_list, frequencyDF)
    DefineLineages.printVariantsTSV(variantDataMeansSD)
    variant_dataForHtml = PlotVariantFrequencies.plotVariants(sample_list,  detailedCounts, detailedCountsMask, variantDataMeansSD , variantDataMeansSDUncorrected )
    #variant_dataForHtml = DefineLineages.LineageCalculator(sample_list, frequencyDF).plotVariants()

    plotrg = "all" if (settings.plotRange is None or (settings.plotRange[0] < 2) & (settings.plotRange[1] > 29902)) else str(settings.plotRange[0]) + "-" + str(settings.plotRange[1])
    minIndelFreq = settings.minFreq if settings.minIndelFrequency < settings.minFreq else settings.minIndelFrequency
    # samples_gridded = np.array_split([s.sample for s in sample_list],3)
    x = [s for s in sample_list]
    samples_gridded = [x[4 * i:4 * (i + 1)] for i in range(math.ceil(len(x) / 4))]

    t = Template(JinjaTemplate.template)
    htmldata = t.render(samples=samples_gridded,
                 nodes=variants.globals.variantdict,
                 plotrg=plotrg,
                 heatmap=heatmap_fig.to_html(full_html=False, include_plotlyjs='cdn'),
                 clustermap=clusterMap_fig.to_html(full_html=False, include_plotlyjs='cdn'),
                 variantsHeatmap=heatmapVariantsOnly.to_html(full_html=False, include_plotlyjs='cdn'),
                 depth=coveragePlots_fig.to_html(full_html=False, include_plotlyjs='cdn') if coveragePlots_fig is not None else None,
                 violins=violinCoveragePlots_fig.to_html(full_html=False, include_plotlyjs='cdn') if violinCoveragePlots_fig != None else None,
                 shannon=shannonPlots_fig.to_html(full_html=False,
                                                  include_plotlyjs='cdn') if shannonPlots_fig is not None else None,
                 shannonBoxFig=shannonBoxFig.to_html(full_html=False,
                                                     include_plotlyjs='cdn') if shannonBoxFig is not None else None,
                 pepperfig=pepperHistogram_fig ,
                 pepperDF=pepperDF ,
                 variantData = variant_dataForHtml,
                 settings = settings,
                 )

    with open(settings.rootDir + "/" + settings.htmlOutFile, 'w') as f:
        f.write(htmldata)
#######################################  MAIN ####################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputDir", default=False, help="Input Directory", type=str)
    parser.add_argument("--analyzemutList", default=False, action="store_true")
    parser.add_argument("--parseProtMutationFile", default=False,
                        help="parse File with amino acid muts and generate file with nucleic acid muts", type=str)
    args = parser.parse_args()
    if args.inputDir:
        settings.rootDir = args.inputDir
    if args.analyzemutList:
        MutlLstProcessing.ProcessMutlist.processMutList(settings.rootDir)
    elif args.parseProtMutationFile:
        Cd.getNucMutsForAAmutsFile(args.parseProtMutationFile)
    else:
        analyzeSamplesMain()
    print(f"The program took {time.time() - start_time} seconds to complete.")
    if not Globals.warningList.empty() or not Globals.errorList.empty():
        print("\033[93m")
        if not Globals.warningList.empty():
            print("+++++++++++++++++++++++++++++++++++++  Warnings:")
            while not Globals.warningList.empty():
                item = Globals.warningList.get()
                print(item)
        if  not Globals.errorList.empty():
            print("+++++++++++++++++++++++++++++++++++++  Errors:")
            while not Globals.errorList.empty():
                item = Globals.errorList.get()
                print(item)
        print("\033[0m")

