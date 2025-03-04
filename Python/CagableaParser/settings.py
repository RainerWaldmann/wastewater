

#import seqtools.GTF as gtf
import multiprocessing
import threading

#INPUT DIRECTORY SUBFOLDERS (SAMPLES) STARTING WITH '_' WILL BE IGNORED
#can be supplied here or as command line argument
rootDir = r"D:\Wastewater\DELETE\juejun"

#if set to True, will remove the first characters from the sample name until the first underscore at max pos 5 in the name. Used to define sorting order
# for example a sample (folder name) AAA_..... will be shown before AAB_ .
# AAA_ and AAB_ will be removed in graphs, report .. if this is set to True
removeFirstCharsFromSampleName = True

# if True, will try to parse date in format YYYY-MM-DD from sample names and use it as x-axis in certain plots
useDateAxis = False

# names of input files. The files in the sample folders
depthExtension = 'depth.tsv'
ivarExtension = 'ivar.tsv'
pepperCountExtension = 'ppmov.txt' # File only needs to be present when ppMov levels should be analyzed
trimstatsExtension = 'TrimStats.txt' # file only required if ppMov amplicon was added.



spikeRange = (21563, 25384)
plotRange = None #If none will plot all, e.g. (21563, 25384) will just plot data for the spike range
#plotRange = (21563,25384)
minDepth = 10 # minimum total sequencing depth at a position for a variation to be considered
minFreq = 0.1 # variations that do not have at least  this frequency in any of the samples will be filtered out
minFreqForHeatMaps = 0.1 # variations that do not have at least  this frequency in any of the samples won't be displayed on heat maps
# Nanopore reads have a higher error rate in particular indels. Those filters allow to eliminate some noise
minFreqForNonMultiplesOfThree = 0.1
filterIndelsNotMultipleOfThree = False # if true will completely filter non multiples of three
minIndelFrequency = 0.1 # indels below this freq will be filtered out
filterInsertions = False # completely ignore insertions


maxAmplifiedRange = (47, 29826) #  max range that is potentially amplified with the used primers  after trimming
sarsCov2length = 29903

#READ QUALITY FILTER
conditionALT_QUAL = lambda iv_data: (iv_data['ALT_QUAL'] < 18) & (iv_data['REF_QUAL'] - iv_data['ALT_QUAL'] > 7) # condition for filtering out low quality variant calls

#FWD/REV IMBALANCE Condition for filtering out
andersonLabIvarUsed = False # if True, ivar files were generated with Ivar from the Anderson lab, if False own software was used. Some quality filters (like indel QVs) cannot be used when ivar was used to define variation frequencies
do_filter_FWD_REV_balance = True
maxALT_FwdRevImbalance = 4.0 # maximum fold difference between ALT forward and ALT reverse reads
condition_FWD_REV_balance = lambda iv_data: (
    (iv_data['ALT_DP'] - iv_data['ALT_RV'] == 0) | # both forward and reverse reads should back variation
    (iv_data['ALT_RV'] == 0) |
    ((iv_data['ALT_DP'] - iv_data['ALT_RV']) / iv_data['ALT_RV'] > maxALT_FwdRevImbalance) |
    (iv_data['ALT_RV'] / (iv_data['ALT_DP'] - iv_data['ALT_RV']) > maxALT_FwdRevImbalance)
)

# various output files
htmlOutFile = "summary.html"
frequenciesOutFile = "VariantFrequencies.xlsx" # excel file with lineage frequencies


#HEATMAPS
hoverinfo_Onheatmaps = True # if True, will show the sample name and the frequency in the hover info
heatmap_pdf = "heatmap.pdf" # if set, will write a pdf with the heatmap, if None doesn't write
heatmap_variantsonly_pdf = "heatmapVariants.pdf" # if set, will write a pdf with the heatmap of variants only, if None doesn't write


#HISTOGRAM
variantHistogram_sumChildsNotExceedParentCorrection = True  # if true, assures that sum of childs does not exceed parent. Pie charts and histograms don't show if child sum is very slightly exceeding parent
variant_hisogram_pdf = "Varianthisto.pdf" # if set, will write a pdf with the histograms, if none won't write

#DEPTHS
#plotDepths = True
downsampleFactorDepths = 10 # scatterplots are slow to generate. downsampling factor for depth speeds up scatter plots. None -> no downsample
depthPDFfilename = None # If set, will write a pdf with the coverage plots.
plotDepthsLogY = True
plotViolinDepths = True

#PEPPER, only relevant if ppMov amplicon was added to Artic pool 1
plotPepper =  False

#PIECHARTS
plotVariantPies = True # generates variant pie charts and addds them to html. Needs also to be True for pie chart PDF writing
sortPieWedges = True # if True, will sort wedges in pie chart
variantPiesPDFfilename = None # "variantPies.pdf" # creates one PDF per sample if None, no pdf will be generated

#If set html will contain box plots with frequencies for individual mutations defining a variant for each variant found in each sample. Helps to visualize outliers
plotBoxPlotsForDetailedVarCounts = False


#ENTROPY
doShannon = False # calculates mean shannon entropy and generates file Shannon.tsv
plotShannon = False # if True plots shannon entropy for each sample in html
shannonPlotSpikeOnly = True # if True, will only plot shannon entropy for the spike region

#Variant frequencies stacked line chart
plot_variant_freqs_line_chart = True # if True, will plot a stacked line chart with variant frequencies (VariantFreqLineChart.pdf)
smoothen_variantsLinePlot_with_filter = True    # if True, will smoothen the line plot with gaussian filter
smoothen_variantsLinePlot_sigma = 2 # sigma for gaussian filter

"""Outlier filtering 
Filters outliers of variation frequencies defining a lineage
filters out variants that are outside of the interquartile range (above max and below min)
max = qA + (interq_range_outliers * intr_qr)
min = qB - (interq_range_outliers * intr_qr)
"""
interq_range_outliers = 0.1


#verbose = True
showprocesstime = True
num_cores = multiprocessing.cpu_count()


# if this is not None, will only display the listed variations on the heatmap
#heatmapDisplayOnlyVariants = None

#
#JUSTTOKEEPIT_

heatmapDisplayOnlyVariants = \
[
    'C21618T',
    'C21618G',
    'C21622T',
    'G21624C',
    '21633del9',
    '21652del3',
    'C21711T',
    'G21718T',
    'C21721T',
    'C21762T',
    '21765del6',
    'T21810C',
    'C21846T',
    'G21941T',
    'G21987A',
    '21987del9',
    '21991del3',
    'C22000A',
    'A22001G',
    '22029del6',
    'A22101T',
    'C22109G',
    '22194del3',
    'T22200G',
    'T22200A',
    'C22208T',
    'C22295A',
    'G22317T',
    'A22320G',
    'G22331A',
    'C22353A',
    'C22480T',
    'A22556G',
    'C22570T',
    'G22577C',
    'G22578A',
    'G22599C',
    'G22599A',
    'A22629C',
    'C22664A',
    'T22673C',
    'C22674T',
    'T22679C',
    'C22686T',
    'A22688G',
    'G22770A',
    'G22775A',
    'A22786C',
    'C22792T',
    'T22795G',
    'G22813T',
    'T22882G',
    'A22893C',
    'G22895C',
    'T22896A',
    'T22896C',
    'G22898A',
    'A22910G',
    'C22916T',
    'T22917G',
    'T22926C',
    'G22927T',
    'T22928C',
    'T22930A',
    'T22942G',
    'T22942A',
    'G22992A',
    'C22995G',
    'C22995A',
    'T23005A',
    '23008del3',
    'G23012A',
    'A23013C',
    'T23018G',
    'T23018C',
    'T23019C',
    'T23031C',
    'C23039G',
    'A23040G',
    'G23048A',
    'A23055G',
    'A23063T',
    'T23075C',
    'C23123T',
    'C23202A',
    'G23222A',
    'C23271T',
    'C23277T',
    'A23403G',
    'C23423T',
    'C23525T',
    'T23599G',
    'C23604A',
    'C23604G',
    'C23854A',
    'G23948T',
    'C24130A',
    'C24378T',
    'G24410A',
    'A24424T',
    'T24469A',
    'C24503T',
    'G24872T',
    'C24990T',
    'C25000T',
    'C25207T'
]



heatmapDisplayOnlyVariantsOLD =  \
[
    '21608insTTATGCCGCTGT', '21608insTCATGCCGCTGTT', 'C21618G', 'C21618T', 'C21622T',
    'G21624C', '21633del9', '21652del3', 'C21711T',
    'G21718T',  'C21762T', '21765del6', 'T21810C',
    'C21846T', 'G21941T', 'G21987A', '21987del9', '21991del3',
    'T22200A', 'T22200G', 'C22208T', 'C22295A', 'G22317T',
    'A22320G',  'C22353A', 'C22480T',
     'A22556G', 'C22570T', 'G22577C', 'G22578A',
    'G22599A', 'G22599C', 'A22629C', 'C22664A', 'T22673C',
    'C22674T', 'T22679C', 'C22686T', 'A22688G', 'G22770A',
    'G22775A', 'A22786C', 'C22792T', 'T22795G', 'G22813T',
    'T22882G', 'A22893C', 'G22895C', 'T22896C', 'T22896A',
    'G22898A', 'A22910T', 'C22916T', 'T22917G', 'T22926C',
    'G22927T', 'T22928C', 'T22930A', 'T22942A', 'T22942G',
    'G22992A', 'A22994G', 'C22995A', 'T23005A', '23008del3',
    'G23012A', 'A23013C', 'T23018G', 'T23018C', 'T23019C',
    'T23031C', 'A23040G', 'G23048A', 'C23202A', 'G23222A',
    'C23271T', 'C23423T', 'C23525T', 'T23599G', 'C23604A',
    'C23604G', 'C23854A', 'G23948T', 'C24130A', 'C24378T',
    'G24410A', 'A24424T', 'T24469A', 'C24503T', 'G24872T',
    'C24990T', 'C25000T', 'C25207T'
]



