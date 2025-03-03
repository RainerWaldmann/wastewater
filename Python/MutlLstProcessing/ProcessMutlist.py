
import pandas as pd
import MutlLstProcessing.settings as settings
from pathlib import Path
import os
import variants.MutInfo as mi
import re



regexformutsplit = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")

def processMutList(rootdir: str):
    """takes tsvs generated from https://covariants.org/ to define muts for variant.
    Only variants in the root directory are listed. Variants in subdirectories are treated
    as if they were the variant in the first level dir"""
    mutant_list = list()
    dummy = os.listdir(rootdir.encode())
    directories = [f for f in Path(rootdir).glob('**/*') if f.is_dir()]
    for dirCounter , currentDir in enumerate(directories):
        #base_name = os.path.basename(d.as_posix())
        #fileCounter = -1
        files = [f for f in os.listdir(currentDir) if os.path.isfile(str(currentDir)+'/'+f)]
        for fileCounter, currentFile in enumerate(files):
                #fileCounter +=1
                print("Reading: " + currentFile)
                if "XBB" in currentFile:
                    dummy = 1
                try:
                    p = pd.read_table(os.path.join(currentDir, currentFile), header=None)
                except:
                    print (currentFile + " threw exception !!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                    os._exit(1)

                p = p.set_index(p.iloc[:, 0] + p.iloc[:, 1])  # set sum of first and second column as index
                substlist = p.loc['nucchanges', 2].strip().split(",") if ('nucchanges' in p.index) else []
                reversionlist = p.loc['nucreversionsToRoot', 2].strip().split(",") if (
                            'nucreversionsToRoot' in p.index) else []
                dellist = p.loc['nucgaps', 2].strip().split(",") if ('nucgaps' in p.index) else None
                dels = []
                if dellist:
                    previousDelPos = -1
                    currentDelLength = 0
                    # delstart = 0
                    for m in dellist:
                        delpos = int(re.findall(r'\d+', m)[0])
                        if delpos != previousDelPos + 1:
                            if previousDelPos > 0:
                                dels.append(str(delstart) + "del" + str(currentDelLength))
                            delstart = delpos
                            currentDelLength = 1
                        else:
                            currentDelLength += 1
                        previousDelPos = delpos
                    if currentDelLength != 0:  # add last del to list
                        dels.append(str(delstart) + "del" + str(currentDelLength))

                x = pd.DataFrame(substlist +
                                # reversionlist + # Todo treat reversion list seperately
                                 dels)
                #x.columns = ["Mut" + str(fileCounter)]
                x.columns = [currentFile]
                x.index = x[currentFile]
                if fileCounter == 0: #first file in directory
                    merge = x
                else:
                    merge = pd.concat([merge, x], axis=1)
        #merge["AAchange"] = [mut if "del" in mut else (mi.MutInfoSubst.getMutInfoFromMutString(mut)).getAAmutstring() for mut
        merge["AAchange"] = [mut if "del" in mut else (mi.getMutInfo(mut)).getAAmutstring() for mut
                             in merge.index]
        # sort by pos

        merge['indexNumber'] = [int(re.findall(r'\d+', i)[0]) for i in merge.index]  # add helper column
        merge.sort_values(by=['indexNumber'], inplace=True)
        merge.drop(columns=['indexNumber'], axis='columns', inplace=True)  # delete helper column
        # def dummy(x):
        #     u = x.to_series().str.findall(r'\d+')
        #     print('AAAAAAAAAAAAAAAA',u)
        #     return x
        # merge.sort_index(key=lambda x: dummy(x),
        #                          inplace=True)
        # merge = merge.sort_index(key=lambda x: x.to_series().str.extract(r'[a-zA-Z]*(\d+)[a-zA-Z]*').astype(int), inplace=True)
        # merge.sort_index(key=lambda x: (regexformutsplit.match(x.to_series().str).groups()[1].astype(int)), inplace=True)
        if len(files) > 1:
            merge.to_csv(rootdir + "\\" + os.path.basename(currentDir)  + ".tsv", sep="\t")
        # filter muts that are not in majority of files for this variant  -> accept 20% NAs if at least 10, one na if at least 5, if less than 5 accept only mutations found in all
        naThreshold = int(len(merge.columns)*0.8) if len(merge.columns) >= 10  else len(merge.columns) -1 if len(merge.columns) >= 5 else len(merge.columns)
        merge = merge.dropna(thresh= naThreshold)
        merge = merge[["AAchange"]]
        merge = merge.rename({'AAchange': os.path.basename(currentDir)}, axis='columns')
        if dirCounter == 0:
            summarydf = merge
        else:
            summarydf = pd.concat([summarydf, merge], axis=1)
    summarydf['indexNumber'] = [int(re.findall(r'\d+', i)[0]) for i in
                                    summarydf.index]  # add helper column with just positions
    summarydf.sort_values(by=['indexNumber'], inplace=True)
    summarydf.drop(columns='indexNumber', axis='columns', inplace=True)  # delete helper column
    summarydf = summarydf[sorted(summarydf.columns)]
    summarydf.to_csv(rootdir + "\\" + "summary.tsv", sep="\t")
    # write csv with muts found only in one variant

    mask = summarydf.apply(lambda x: x.count()==1, axis=1)
    summarydf[mask].to_csv(rootdir + "\\" + "summaryUnique.tsv", sep="\t")



