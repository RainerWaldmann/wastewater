import Globals
import pandas
import os
import numpy as np
import settings
import Data
from pathlib import Path
import math
import re
import variants.MutInfo as mi
import variants
import pandas as pd
import gc
from joblib import Parallel, delayed
import time
import dateutil.parser as dparser
import sys


def readData():
    start_time = time.time()

    def readOneDirectory(f):
        fileList = [p.as_posix() for p in f.iterdir() if p.is_file()]
        print("Parsing: " + str(f))

        if settings.ivarExtension in str(fileList) and settings.depthExtension in str(fileList):

            depth_data = pd.read_csv(next(name for name in fileList if settings.depthExtension in name), sep='\t',
                                     names=["REF", "POS", "COUNT"],
                                     dtype={"POS": "int", "COUNT": "int"}, usecols=["POS", "COUNT"])
            #fill in 0 depth for lacking positions 4 steps, TODO find easier way
            all_positions = pd.DataFrame({
                'POS': range(1, settings.sarsCov2length),
                'COUNT': 0
            })

            # Step 2: Set 'POS' as the index in both DataFrames
            depth_data.index = depth_data['POS']
            all_positions.index = all_positions['POS']

            # Step 3: Combine the DataFrames
            depth_data = depth_data.combine_first(all_positions)

            iv_data = pd.read_table([i for i in fileList if settings.ivarExtension in i][0],
                                    dtype={"REGION": "string", "POS": "int32", "REF": 'category',
                                           "ALT": 'string', "REF_DP": "int32", "REF_RV": "int32",
                                           "REF_QUAL": 'float32', "ALT_DP": "int32", "ALT_RV": "int32",
                                           "ALT_QUAL": 'float32', "ALT_FREQ": "float32", "TOTAL_DP": "int32",
                                           "PVAL": "float16", "PASS": "bool",
                                           "GFF_FEATURE": "string", "REF_CODON": "string", "REF_AA": "string",
                                           "ALT_CODON": "string", "ALT_AA": "string",
                                           "PANGOSTR": "string",
                                           "LINEAGES": "string"})
            # remove unneeeded columns
            iv_data = iv_data.drop(
                ['REGION', "PANGOSTR", "LINEAGES", "GFF_FEATURE", "REF_CODON", "ALT_CODON", "REF_AA", "ALT_AA", "PASS",
                 "PVAL"], errors='ignore', axis=1)


            # in Anderson lab ivar some lines are duplicated when GTF is supplied and CDS in GTF overlap
            dupFilter = iv_data.duplicated(["POS", "REF", "ALT"])
            iv_data = iv_data[~dupFilter]
            initial_iv_data_len = len(iv_data)
            # DO SOME FILTERING
            #remove rows with bad quality or bad forward reverse balance or ALT qual
            conditionALT_QUAL = settings.conditionALT_QUAL(iv_data)
            conditionBal = pd.Series([False]*len(iv_data)) if not settings.do_filter_FWD_REV_balance else settings.condition_FWD_REV_balance(iv_data)
            condition4 = settings.andersonLabIvarUsed or ~iv_data['ALT'].str.contains('del')  # deletions are always forward in ivar -> don't filter rows with deletions
            conditionQualAndBalance = (conditionALT_QUAL | conditionBal) & condition4
            #for now filter deletions that are not multiples of three. They are typically artefacts
            conditionDelNotMultipleOfThree =  ~iv_data['ALT'].str.match(r'(^[+-](\w{3})*$|^[^+-])') if settings.filterIndelsNotMultipleOfThree else False # indels that are not multiples of 3
            conditionFilterInsertions = iv_data['ALT'].str.startswith('+') if settings.filterInsertions else False# remove insertions # TODO treat insertions current xorkaround remove ALT lines that contain '+'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #filter low depth rows . Use depths from the depth dataframe since depths from the Anderson lab Ivar are inconsistent
            depthhighenough_pos = depth_data[depth_data["COUNT"] >= settings.minDepth].index
            conditionSeqDepth = ~iv_data['POS'].isin(depthhighenough_pos)
            iv_data = iv_data.drop(iv_data[conditionQualAndBalance |
                                           conditionDelNotMultipleOfThree |
                                           conditionSeqDepth |
                                           conditionFilterInsertions].index)
            print(" ---------- Read data for SAMPLE: ", os.path.basename(f.as_posix()) + "----------\n" +
                  "Frequency table number of rows: "+ str(initial_iv_data_len) + "\n" +
                  "Filtered low ALT_QUAL: " + conditionALT_QUAL.sum().astype(str) + "\n" +
                  "Filtered low ALT FWD/REV balance: " + conditionBal.sum().astype(str) + "\n" +
                  "Filtered Depth too low: " + conditionSeqDepth.sum().astype(str) + "\n" +
                    "Frequency table number of rows after filters : " + str(len(iv_data)) + "\n"
                  )
            #correct del positions. Ivar prints pos before del -> use pos where del starts (increase by one)
            iv_data.loc[iv_data['ALT'].str.startswith('-'), 'POS'] += 1
            #replace deletions such as -NNNNNN by del6
            iv_data.ALT = iv_data.ALT.apply(lambda z: 'del' + str(len(z[1:])) if z[0] == '-'  else z).astype('category')

              # add column with shannon entropy to depth dataframe
            if settings.doShannon:
                dfForShannon = __addWTfrequenciesToDataframe(iv_data) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                depth_data = __computeshannonfromivar(dfForShannon, depth_data)

            # read amplicon data
            file_name = next((name for name in fileList if 'TrimStats.' in name), None)
            if file_name is not None:
                amplicon_data = pd.read_table(file_name)
                trimstats = Data.AmpliconData(amplicon_data)
                # get Artic pool 1 count for pepper ref counts
                artic_data = next((item for item in trimstats.data.items() if 'artic' in item[0].lower()), None)
            else:
                trimstats = None
                artic_data = None

            if artic_data is not None:
                articPool1Count = artic_data[1].pool1Count
            else:
                settings.plotPepper = False
                articPool1Count = None
                print("\033[91m ++++++++   settings.plotPepper set to False because no artic pool 1 count found in: " + f.name +"\033[0m")

            spikeInTrimStatsFileName = [name for name in fileList if 'TrimStats.SpikeIn.' in name]
            if len(spikeInTrimStatsFileName) != 0: # spikeIn Trimstat file exists
                t = pd.read_table(spikeInTrimStatsFileName[0])
                trimstatsSpikeIn = Data.AmpliconData(t)
            else:
                trimstatsSpikeIn = None

            # read file with pepper counts
            ppmov = None
            if settings.plotPepper:
                pepper_file_path = next(f.glob(f"*{settings.pepperCountExtension}"), None)
                if pepper_file_path and os.path.isfile(pepper_file_path):
                    with open(pepper_file_path) as pep:
                        ppmov = int(pep.readline()) #file just contains a number, the ppmov read count
                else:
                    e = "settings.plotPepper set to False because no pepper data found in: " + f.name
                    Globals.warningList.put(e)
                    print("\033[91m+++++++++++++++++++++ " + e  + "+++++++\033[0m")
            #search whether date in format YYYY-MM-DD is in filename

            match = re.search(r'\d{4}-\d{2}-\d{2}', f.name)

            if match:
                date = dparser.parse(match.group(), fuzzy=True)
            else:
                if settings.useDateAxis:
                    e = "Sample: " + f.name + " DATE in format YYYY-MM-DD not found, settings.useDateAxis was true, settings.useDateAxis set to False"
                    Globals.warningList.put(e)
                    print("\033[91m" + e + "\033[0m")  # print error message in red
                    settings.useDateAxis = False
                date = None

            # return data for one sample
            return Data.OneBCdata(f.name, iv_data, depth_data, trimstats, trimstatsSpikeIn,articPool1Count, ppmov, date)


    directories = [f for f in Path(settings.rootDir).iterdir() if f.is_dir() and not f.name.startswith("_")]
    sample_list = Parallel(n_jobs=settings.num_cores)(
            delayed(readOneDirectory)(f) for f in directories)
    sample_list = [item for item in sample_list if item is not None] # filter out None values which are folders that did not have ivar and depth files

    if settings.showprocesstime:
        print(
            f"-----------------Read files finished took {time.time() - start_time} seconds to complete.-------------------------------")
    # Remove first two chars when samples are named A_, B_ ... for sorting
    if settings.removeFirstCharsFromSampleName:
        sample_list.sort(key=lambda x: x.sample)
        for x in range(len(sample_list)):
            if '_' in sample_list[x].sample[0:4]:
                sample_list[x].sample = sample_list[x].sample.split("_", 1)[1]

    settings.useDateAxis = settings.useDateAxis and not any(
        pd.isna(sample.date) for sample in sample_list)  # all samples contain valid date

    if settings.useDateAxis:
        for s in sample_list:
            s.sample = s.date.date().isoformat()
    return sample_list


def __mergeCloseDeletions (data:pd.DataFrame) -> pd.DataFrame:
    """NOT USED AND IMPLEMENTED YET merge deletions that are in close positions and have similar lengths UNFINISHED !!!!!!!!"""
    df = data[data['ALT'].str.contains('del')]
    finishedIndices = {}
    obsoleteIndices = {} # indices that were merged
    previousRow = None
    #df.loc[df_a.index, :] = df_a[:] #
    for i in df.index:
        row = df.loc[i]
    return data


def mergeAndFilterDataTEST(sample_list):
    """merge ivar dataframes"""
    columns_tomerge = ["POS", "REF", "ALT"]

    othercolumns = ['REF_DP', 'REF_RV', 'REF_QUAL', 'ALT_DP', 'ALT_RV', 'ALT_QUAL', 'ALT_FREQ', 'TOTAL_DP']
    merged_df_list = [sample_list[0].ivar.rename(columns={col: f"{col}_{sample_list[0].sample}" for col in othercolumns if col in sample_list[0].ivar.columns})]

    if len(sample_list) == 1:
        return merged_df_list[0]

    # Prepare the suffixes and list for merging
    suffixes = [f"_{sample.sample}" for sample in sample_list]

    for i in range(1, len(sample_list)):
        print("Merging sample ", i + 1)

        suffix = suffixes[i]
        temp_df = sample_list[i].ivar

        # Convert temp_df to efficient data types
        temp_df = temp_df.astype({
            "POS": "int32", "REF": "category", "ALT": "category",
            "REF_DP": "int32", "REF_RV": "int32", "REF_QUAL": "float32",
            "ALT_DP": "int32", "ALT_RV": "int32", "ALT_QUAL": "float32",
            "ALT_FREQ": "float32", "TOTAL_DP": "int32"
        })

        # Free RAM used by the ivar df
        sample_list[i].ivar = None

        # Add suffix to all columns except those in columns_tomerge
        temp_df = temp_df.rename(columns={col: col + suffix for col in temp_df.columns if col not in columns_tomerge})

        merged_df_list.append(temp_df)

    # Perform the merge
    mer = pd.concat(merged_df_list, axis=0, join='outer', ignore_index=True).groupby(columns_tomerge,
                                                                                     as_index=False).first()
    #filter columns with too low freqs
    alt_freq_cols = [col for col in mer.columns if col.startswith('ALT_FREQ')]
    mer = mer[(mer[alt_freq_cols] >= settings.minFreq).any(axis=1)]
    # Reduce RAM usage - convert columns to int32
    c = mer.filter(regex='^(ALT_DP_|REF_DP_|REF_REV_|TOTAL_DP)').columns
    #print("**********MEM",mer.dtypes,mer.memory_usage())
    mer[c] = mer[c].astype(pd.Int32Dtype())
    #mer[c] = mer[c].apply(pd.to_numeric, errors='coerce').astype(pd.Int32Dtype()) # does not change still remains float64
    mer["ALT"] = mer["ALT"].astype("category") # alt became object after merge -> reset to category to save RAM
    #print("**********MEM",mer.dtypes,mer.memory_usage())
    # Sort by pos
    mer.sort_values(by=['POS'], inplace=True)
    mer.reset_index(drop=True, inplace=True)
    mer = mer.copy() # copy to free memory
    if mer.empty:
        print("Error: No mutation passed filter (min reads, qual ....). ---- EXITING")
        sys.exit(1)
    print("DONE mergeData")
    return __filterMergedData(mer, sample_list)


def mergeAndFilterData(sample_list):
    """merge ivar dataframes"""
    mer = sample_list[0].ivar
    columns_tomerge = ["POS", "REF", "ALT"]
    if len(sample_list) == 1:  # if just one sample merge won't add sample name to columns -> need to rename columns
        othercolumns = ['REF_DP', 'REF_RV', 'REF_QUAL', 'ALT_DP', 'ALT_RV', 'ALT_QUAL', 'ALT_FREQ', 'TOTAL_DP']
        for col in othercolumns:
            if col in mer.columns:
                mer = mer.rename(columns={col: col + '_' + sample_list[0].sample})

    for i in range(1, len(sample_list)):
        print("Merging sample ", i + 1)
        bc1 = "_" + sample_list[i - 1].sample
        bc2 = ("_" + sample_list[i].sample) if i == len(sample_list) - 1 else None
        mer = mer.merge(sample_list[i].ivar, on=columns_tomerge, how='outer'
                        , suffixes=(bc1, bc2))
        # the following two lines free RAM used by the ivar df
        sample_list[i].ivar = None
        gc.collect()
    # filter columns with too low freqs
    alt_freq_cols = [col for col in mer.columns if col.startswith('ALT_FREQ')]
    mer = mer[(mer[alt_freq_cols] >= settings.minFreq).any(axis=1)]
    # Reduce RAM usage - convert columns to int32
    c = mer.filter(regex='^(ALT_DP_|REF_DP_|REF_REV_|TOTAL_DP)').columns
    # print("**********MEM",mer.dtypes,mer.memory_usage())
    mer[c] = mer[c].astype(pd.Int32Dtype())
    # mer[c] = mer[c].apply(pd.to_numeric, errors='coerce').astype(pd.Int32Dtype()) # does not change still remains float64
    mer["ALT"] = mer["ALT"].astype("category")  # alt became object after merge -> reset to category to save RAM
    # print("**********MEM",mer.dtypes,mer.memory_usage())
    # Sort by pos
    mer.sort_values(by=['POS'], inplace=True)
    mer.reset_index(drop=True, inplace=True)

    print("DONE mergeData")
    return __filterMergedData(mer, sample_list)


def __filterMergedData(mergedData, sample_list) -> pandas.DataFrame:
    """Filters merged dataframe"""

    # def __filterNAsAndBelowCutoff(x: pd.Series, cutoff: float, sample_list: list) -> bool:
    #     """used to Filter lines with just NAs and lines where no freq is above cutoff freqs. Might not be necessary"""
    #     return any(
    #         (not pd.isna(x["ALT_FREQ_" + sample.sample]) and float(x["ALT_FREQ_" + sample.sample]) >= cutoff) for sample
    #         in sample_list)

    # TODO treat insertions current xorkaround remove ALT lines that contain '+'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #mergedData = mergedData[~mergedData['ALT'].str.contains(r'\+', na=False)]
    #passFilter = mergedData.apply(lambda x: __filterNAsAndBelowCutoff(x, settings.minFreq, sample_list), axis=1)
    #mergedData = mergedData[passFilter]
    # fill up missing values (NA in ALT_FREQ_) set to 0 if NA and enough sequencing depth -> was not in varfreq file because ALT_FREQ low
    # first replace depth values in varfreq file with values from depth file  
    for i in range(0, len(sample_list)):
        cov = sample_list[i].depths.copy() # depth plots require still initial depth data
        # rename columns --> depth column in coverage and depth column in merge need to have same name
        col_toupdate = f'TOTAL_DP_{sample_list[i].sample}'
        cov.rename(columns={'COUNT': col_toupdate}, inplace=True)
        cov.set_index('POS',inplace=True )
        # should have coverage dataframe with same columns / rows as ivar dataframe
        # replace NA in merged ivar with data from depth dataframe  see https://stackoverflow.com/questions/55093574/fill-nans-within-1-column-of-a-df-via-lookup-to-another-df-via-pandas , https://stackoverflow.com/questions/29177498/python-pandas-replace-nan-in-one-column-with-value-from-corresponding-row-of-sec
        mergedData[col_toupdate] = np.where(mergedData[col_toupdate].isnull(),
                                           mergedData['POS'].map(cov[col_toupdate])  , # val to use if NA
                                           mergedData[col_toupdate]) # value to use if not NA
        # update missing frequency values
        # 0 if FREQ for pos is NA and depth > mindepth otherwise FREQ
        sample = sample_list[i].sample
       
        mergedData[f'ALT_FREQ_{sample}'] = mergedData.apply(
            lambda x: 0 if np.isnan(x[f'ALT_FREQ_{sample}']) and x[f'TOTAL_DP_{sample}'] > settings.minDepth else x[f'ALT_FREQ_{sample}'],
            axis=1)
    #filter bad forward reverse balance
    #mergedData = mergedData.apply(__correctFrequenciesWithBadFwdRevBalance, axis=1, args=(sample_list,))
    print("DONE mergeFilterData")
    return mergedData


def getMergedDataForHeatMap(dfp:pandas.DataFrame, sample_list: list[Data.OneBCdata]) -> pandas.DataFrame :
    """filters merged dataframe for heatmaps"""
    merForHeatMaps = dfp[dfp.apply(
        lambda x: any(float(x["ALT_FREQ_" + sample.sample]) >= settings.minFreqForHeatMaps for sample in sample_list),
        axis=1)]
    # passFilterForHeatMaps = dfp.apply(
    #     lambda x: __filterNAsAndBelowCutoff(x, settings.minFreqForHeatMaps, sample_list), axis=1)
    #merForHeatMaps = dfp[passFilterForHeatMaps]
    print("length dataframe after filter ", len(dfp))
    print("length dataframe for heatmaps after filter ", len(merForHeatMaps))
    return merForHeatMaps


def getFrequencyColumnsOnly(df : pandas.DataFrame, sample_list : list, range : tuple) -> pandas.DataFrame:
    """extract requency columns into new dataframe, range is a tuple with a 1-based begin and end of range. If range is None the entire df is used
    modifies index and adds info on amino acid substitution to index"""
    retval = df[["ALT_FREQ_" + x.sample for x in sample_list]]
    retval.columns = retval.columns.str.split('FREQ_').str[-1]  # remove ALT_FREQ from column names
    # TODO treat insertions current workaround remove ALT lines that contain '+'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    def __generate_index(ref, pos, alt):
        if 'del' in alt.lower(): # deletion
            return f"{pos}{alt}{mi.MutInfoDeletion.getMutInfoFromPosMut(pos, alt).getAAmutstringForDFindex()}"
        elif len(alt) == 1:# substitution
            return f"{ref}{pos}{alt}{mi.MutInfoSubst.getMutInfoFromPosMut(pos, alt).getAAmutstringForDFindex()}"
        elif '+' in alt: # insertion
            return f"{pos}ins{alt[1:]}"
        else:
            return ""

    retval.index = [__generate_index(ref, pos, alt) for ref, pos, alt in zip(df['REF'], df["POS"].astype(int), df["ALT"])]
    # if range e.g. Spike was specified
    if range is not None:
        positions = retval.index.str.extract(r'(\d+)', expand=False).astype(int)
        retval = retval[positions.between(range[0], range[1])]

    return retval


def __computeshannonfromivar(dfivar: pd.DataFrame, dfdepth: pd.DataFrame) -> pd.DataFrame:
    """uses dataframe where wt freqs are retained and depth data to calculate shannon entropy"""
    #ignore deletions because mapping issues always lead to high entropy
    d = dfivar[~dfivar['ALT'].str.match(r'(^del)')].reset_index()
    iv = d.groupby('POS')
    #use log4 for shannon entropy
    shannon = iv.agg(shannon = ('ALT_FREQ', lambda x : -sum([i*math.log(i,2) for i in x if i > 0 ])))
    #nature paper
    # def shannonFct(x):
    #     p= sum([i  for i in x ]) # can't be right this way
    #     return 0 if p <= 0 or p >= 1 else -p* math.log(p, 2) - (1-p) * math.log(1-p, 2) # according to san diego nature paper
    # shannon = iv.agg(shannon=('ALT_FREQ', shannonFct)).reset_index()
    dfdepth ['Shannon'] = dfdepth.index.map(shannon['shannon'])
    dfdepth['Shannon'] = dfdepth ['Shannon'].fillna(0)
    return dfdepth




def __addWTfrequenciesToDataframe(df1: pd.DataFrame) -> pd.DataFrame:
    """Adds rows with wildtype frequencies to ivar dataframe
    currently only used for Shannon Entropy"""
    def __generateWTfreqRow(prevRow, frequencysum : float):
        """uses info from row and generates row for wt frequency"""
        return {"POS": int(prevRow.POS), "REF": prevRow.REF, "ALT": prevRow.REF,
             "ALT_FREQ": (1 - frequencysum) if frequencysum <= 1 else 0}


    # filter indels that are not multiples of three and have a freq < 2%
    df = df1[(df1['ALT'].str.startswith('del') == False) ]
            #df1['ALT'].str.match(r'(^del(\w{3})*$)') |
            #pd.to_numeric(df1["ALT_FREQ"]) > 0.02]
    df = df[['POS', 'REF','ALT', 'ALT_FREQ']]
    previousRow = None
    rows_list = []
    freqsum = 0

    for row in df.itertuples():
        pos = row.POS
        if previousRow is not None and int(previousRow.POS) != pos:  # new pos
            rows_list.append(__generateWTfreqRow(previousRow, freqsum))
            freqsum = 0
        rows_list.append(row._asdict())
        previousRow = row
        freqsum += float(row.ALT_FREQ)
    rows_list.append(__generateWTfreqRow(previousRow, freqsum))

    new_df = pd.DataFrame(rows_list).drop(columns=['Index'])

    return new_df


#currently NOT USED !!!!!!!!!!!!!!!!!!!
def __addWTfrequenciesForAllPositionsToDataframeNEW(df: pd.DataFrame, depths: pd.DataFrame) -> pd.DataFrame:
    """Adds rows with wildtype frequencies to ivar dataframe
    Also adds WT frequencies of 1 for all pos with coverage even if they are not in the ivar file !!!!!"""
    # todo add depth data
    def generateWTfreqRow(prevRow: pd.Series, pos:int, alt:int, ref:int, frequencysumothers : float):
        """ generates row for wt frequency"""
        wt_row = pd.Series(index=df.columns,dtype=prevRow.dtype)
        wt_row["POS"] = pos
        wt_row["REF"] = ref
        wt_row["ALT"] = alt
        wt_row["ALT_FREQ"] = (1 - frequencysumothers) if frequencysumothers <= 1 else 0
        return wt_row
    # # filter indels that are not multiples of three and have a freq < 2%
    # df = df[(df['ALT'].str.startswith('del') is False) or
    #     df['ALT'].str.match(r'(^del(\w{3})*$)') | pd.to_numeric(df["ALT_FREQ"]) > 0.01]
    #currentpos = None
    lineList = list()
    previousRow = None
    freqsum = 0
    for r in df.iterrows():
        row = r[1]
        pos = int(row["POS"])
        #generate line with wt freqs
        if previousRow is not None and int(previousRow["POS"]) != pos: # new pos
            r = generateWTfreqRow(previousRow,previousRow["POS"], previousRow["REF"], previousRow["REF"], freqsum)
            lineList.append(r)
            freqsum = 0
        # add wt line rows for positions without ivar lines
        if previousRow is not None and pos > previousRow['POS'] + 1:
            for x in range(1, pos - previousRow['POS']):
                p = previousRow['POS'] + x
                if p in depths.index and depths.loc[p]['COUNT'] > 0:
                    lineList.append(generateWTfreqRow(previousRow, p, variants.globals.sarscov2seq[p -1], variants.globals.sarscov2seq[p -1],0))
        lineList.append(row)
        previousRow = row
        freqsum = freqsum + float(row["ALT_FREQ"])
    lineList.append(generateWTfreqRow(previousRow,previousRow["POS"], previousRow["REF"], previousRow["REF"], freqsum))
    dfnew = pd.DataFrame(data=lineList,columns=df.columns)
    return dfnew








