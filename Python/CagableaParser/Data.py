import pandas as pd

class AmpliconData:
    "Holds trimstats data"
    def __initialize(self, ts):
        data = dict()
        for i in range(0, len(ts.columns)-3, 6):
            name = ts.columns.values[i]
            pool1Count = int(ts.iloc[0, i+1])
            pool2Count = int(ts.iloc[1, i+1])
            #print(pool1Count,pool2Count)
            #print(name,"****")
            #print(list(ts.iloc[2, i:i+5]))
            df = pd.DataFrame(ts.iloc[3:, i:i+5])
            df.columns = list(ts.iloc[2, i:i+5])
            df. dropna(inplace=True)
            df.index = df["SetNo"]
            df.drop('SetNo', axis =1, inplace = True)
            #print(name)
            #print(df)
            data[name] = OneAmpliconData(name,pool1Count,pool2Count,df)
            #print(ts.iloc[i,0])
        return data
    def __init__(self, trmstats):
        self.data = self.__initialize(trmstats)

class OneAmpliconData:
    def __init__(self, name : str, pool1Count : int, pool2Count : int, dataframe):
        self.name = name
        self.pool1Count = pool1Count
        self.pool2Count = pool2Count
        self.dataframe = dataframe


class OneBCdata:
    "Holds the data for each sample"

    def __init__(self, sample : str,ivar,depths, trimstats, trimstatsSpikeIn, articPool1Count, ppmov, date):
        self.sample = sample
        self.ivar = ivar
        self.depths = depths
        self.ppmov = ppmov
        self.date = date
        self.articPool1Count = articPool1Count
        self.trimstats = trimstats
        self.trimstatsSpikeIn = trimstatsSpikeIn
        self.spikeInData = None # dictionary with amplicon names, sarsCov/spikeIn counts

    def getCoverageAtPos(self, pos):
        """returns coverage at given position"""
        if pos in self.depths.index.get_level_values(0):
            return self.depths.loc[pos].iloc[1]
        else:
            return None



