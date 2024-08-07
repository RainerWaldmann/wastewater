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
import threading
import time

start_time = time.time()

pd.options.display.float_format = '{:,.3f}'.format

import variants.MutInfo


def analyzeSamplesMain():
    def execute_in_threads(functions , args):
        """executes the functions (identical args) in separate threads, waits until all done """
        # Initialize results list
        results = [None] * len(functions)

        def store_result(index, func, func_args):
            results[index] = func(*func_args)

        # Create threads
        threads = [threading.Thread(target=store_result, args=(i, func, args), name=f"Thread-{func.__name__}") for i, func in enumerate(functions)]

        # Start threads
        for thread in threads:
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        return results


    sample_list = ReadData.readData()
    settings.useDateAxis = settings.useDateAxis and not any(
        pd.isna(sample.date) for sample in sample_list)  # all samples contain valid date
    functions = [ReadData.mergeAndFilterData, Coverage.plotCoverage, Coverage.plotCoverageViolinPlots]
    if settings.plotShannon:
        functions.append(Coverage.plotShannon)
        functions.append(Coverage.getShannonBoxPlots)
    start_time = time.time()
    results = execute_in_threads(functions, (sample_list,))
    mergedData = results[0]
    coveragePlots_fig = results[1]
    violinCoveragePlots_fig = results[2]
    shannonPlots_fig = results[3] if settings.plotShannon else None
    shannonBoxFig = results[4] if settings.plotShannon else None
    if settings.showprocesstime:
        print(f"----------------- mergedData, coveragePlots_fig, violinCoveragePlots_fig, shannonBoxFig, shannonPlots_fig finished took {time.time() - start_time} seconds to complete.-------------------------------")
    #write TSV with merged data
    mergedData.to_csv(settings.rootDir + "/Ivar.tsv", sep="\t",index=False)
    # get dataframes with just frequencies
    frequencyDF = ReadData.getFrequencyColumnsOnly(mergedData, sample_list, None)
    frequencyDF.to_csv(settings.rootDir + "/IvarLight.tsv", sep="\t")
    mergedDataForHeatMaps = ReadData.getMergedDataForHeatMap(mergedData, sample_list)
    start_time = time.time()
    functions = [HeatMap.plotHeatMap, HeatMap.plotHeatMapVariantsOnly, HeatMap.clusterMap]
    heatmap_fig, heatmapVariantsOnly, clusterMap_fig = execute_in_threads(functions, (mergedDataForHeatMaps, sample_list))
    if settings.showprocesstime:
        print(f"----------------- heatmap_fig, heatmapVariantsOnly, clusterMap_fig finished took {time.time() - start_time} seconds to complete.-------------------------------")

    if settings.plotPepper:
        pepperHistogram_fig, pepperDF = Pepper.countPepper(sample_list)
    else:
        pepperHistogram_fig = None
        pepperDF = None

    #variant_dataForHtml = Var.plotVariants(sample_list, frequencyDF)
    detailedCounts, detailedCountsMask, variantDataMeansSD, variantDataMeansSDUncorrected = DefineLineages.calcVariantFrequencies(sample_list, frequencyDF)
    DefineLineages.printVariantsTSV(variantDataMeansSD)
    variant_dataForHtml = PlotVariantFrequencies.plotVariants(sample_list,  detailedCounts, detailedCountsMask, variantDataMeansSD , variantDataMeansSDUncorrected )
    #variant_dataForHtml = DefineLineages.LineageCalculator(sample_list, frequencyDF).plotVariants()

    plotrg = "all" if (settings.plotRange is None or (settings.plotRange[0] < 2) & (settings.plotRange[1] > 29902)) else str(settings.plotRange[0]) + "-" + str(settings.plotRange[1])
    minIndelFreq = settings.minFreq if settings.minIndelFrequency < settings.minFreq else settings.minIndelFrequency
    # samples_gridded = np.array_split([s.sample for s in sample_list],3)
    x = [s.sample for s in sample_list]
    samples_gridded = [x[4 * i:4 * (i + 1)] for i in range(math.ceil(len(x) / 4))]

    t = Template(JinjaTemplate.template)
    htmldata = t.render(samples=samples_gridded,
                 nodes=variants.globals.variantdict,
                 minDepth=settings.minDepth,
                 minFreq=settings.minFreq,
                 minFreqForHeatMaps=settings.minFreqForHeatMaps,
                 minIndelFreq=settings.minIndelFrequency,
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

