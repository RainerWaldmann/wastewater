

#import seqtools.GTF as gtf
import multiprocessing

#SUBFOLDERS (SAMPLES) STARTING WITH '_' WILL BE IGNORED
#rootDir = r"D:\Delete\comparaison_Tris_ipmc_21Nov2023\analysisTest"
#rootDir = r"D:\Delete\CagableaParserTest"
#rootDir = r"D:\Delete\PAPER\20221711_cagablea METHODS COMP"
#rootDir = r"D:\Delete\PAPER\20220106_cagablea_VEOLIA"
#rootDir = r"D:\Delete\PAPER\20211216_cagablea_Veolia"
#rootDir = r"D:\Delete\PAPER\Test 20240204"
rootDir = r"D:\Wastewater\SAMPLES"

useDateAxis = True # will try to parse date in format YYYY-MM-DD from sample names and use it as x-axis in certain plots

plotRange = None #(0, 99999999)
#plotRange = (21563,25384)
minDepth = 10
minFreq = 0.003
minFreqForHeatMaps = 0.1
minFreqForNonMultiplesOfThree = 0.3
minIndelFrequency = 0.05 # indels below this freq will be filtered out

depthExtension = 'depth.tsv'
ivarExtension = 'ivar.tsv'
pepperCountExtension = 'ppmov.txt'
trimstatsExtension = 'TrimStats.txt'
spikeRange = (21563, 25384)
maxAmplifiedRange = (47, 29826) #  max range that is potentially amplified with the used primers  after trimming
sarsCov2length = 29903
doPlot = False # if True will pop up individual plots
andersonLabIvarUsed = False # if True, ivar files were generated with Ivar from the Anderson lab, if False own software was used
#READ QUALITY FILTER
conditionALT_QUAL = lambda iv_data: (iv_data['ALT_QUAL'] < 17) & (iv_data['REF_QUAL'] - iv_data['ALT_QUAL'] > 7) # condition for filtering out lo

#FWD/REV IMBALANCE Condition for filtering out
do_filter_FWD_REV_balance = False
maxALT_FwdRevImbalance = 5 # maximum fold difference between ALT forward and ALT reverse reads
condition_FWD_REV_balance = lambda iv_data: (((iv_data['ALT_DP'] - iv_data['ALT_RV']) / iv_data['ALT_RV'] > maxALT_FwdRevImbalance) |
                                             (iv_data['ALT_RV'] / (iv_data['ALT_DP'] - iv_data['ALT_RV']) > maxALT_FwdRevImbalance)) # condition for filtering out

#sarscovgff = None #GFF data for SarsCov2
#sarscov2seq = None #SarsCov2 sequence

htmlOutFile = "summary.html"
frequenciesOutFile = "VariantFrequencies.xlsx"
removeFirstCharsFromSampleName = True
#HEATMAPS
hoverinfo_Onheatmaps = True
#HISTOGRAM
variantHistogram_sumChildsNotExceedParentCorrection = True  # if true, assures that sum of childs does not exceed parent
#DEPTHS
plotDepths = True
downsampleFactorDepths = 10 # scatterplots are slow to generate. downsampling factor for depth speeds up scatter plots. 1
depthPDFfilename = None # if None, no pdf will be generated
plotDepthsLogY = True
plotViolinDepths = True
plotPepper = True
#PIECHARTS
plotVariantPies = True
sortPieWedges = True # if True, will sort wedges in pie chart
variantPiesPDFfilename = None # "variantPies.pdf" # creates one PDF per sample if None, no pdf will be generated
plotBoxPlotsForDetailedVarCounts = False
#ENTROPY
plotShannon = False
#Variant frequencies stacked line chart
plot_variant_freqs_line_chart = True
smoothen_with_savgol_filter = True

verbose = True
showprocesstime = True
num_cores = multiprocessing.cpu_count()
