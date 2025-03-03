import pandas as pd
import plotly.express as px
import settings

def countPepper(sample_list):
    #peppercounts = list()
    peppercountsTest =  list()
    peppersamples = list()

    for s in sample_list:
        #print ([name for name in s.trimstats.data.keys() if name in pepperPool][0])
        if s.ppmov==s.ppmov and s.ppmov != 0:
            ratioSarsPpmov = (int(s.articPool1Count) / int(s.ppmov))
        else:
            ratioSarsPpmov = s.ppmov
        peppercountsTest.append(ratioSarsPpmov)
        peppersamples.append(s.sample)

        #if len([name for name in s.trimstats.data.keys() if name in settings.pepperPool]) > 0 :
            #peppersamples.append(s.sample)
            #peppercounts.append(s.ppmov if s.ppmov == 0 else s.trimstats.data[[name for name in s.trimstats.data.keys() if name in settings.pepperPool][0]].pool1Count / s.ppmov)
    #print(peppercountsTest)

    pepperDataFrameTest = pd.DataFrame(peppercountsTest,index=peppersamples, columns=["SarsCov/ppmov ratio"])
    #pepperDF = pd.DataFrame(peppercounts)
    #if pepperDF.empty:
    #    return None, None
    #pepperDF.index = peppersamples
    #pepperDF.columns = ['SarsCov/ppmov ratio']
    #pepperDF['SarsCov/ppmov ratio'] = pepperDF['SarsCov/ppmov ratio'].round(decimals=2).astype(str)
    #print(peppercounts)
    #pepperDF.info()
    #display(pepperDF)
    #import plotly.graph_objects as go
    #fig = go.Figure(data=[go.Table(#header=dict(values=['Sample', 'SARSCov/ppmov']),
    #                 cells=dict(zip(pepperDF.index,pepperDF[0])))
    #                     ])

    pepperfig = px.histogram(pepperDataFrameTest, x=pepperDataFrameTest.index, y='SarsCov/ppmov ratio', text_auto=".3s" #, title="Sars-Cov-2/pepper count ratios"
        )
    pepperfig.update_layout(
        width=200 + len(pepperDataFrameTest.index)*50,
        height=800,
        title="",
        xaxis_title="Sample",
        yaxis_title='SarsCov/ppmov ratio',
        hovermode = False
    )

    return pepperfig, pepperDataFrameTest