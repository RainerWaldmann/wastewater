import pandas
import pandas as pd
import variants.MutInfo as mi
import re


#static variables (not used here): https://stackoverflow.com/questions/68645/static-class-variables-and-methods-in-python
from variants.MutInfo import MutInfoSubst

class ChimeraInfo:
    """Holds info on junction, parent variants, min mutations that need to be found before and after junction
    Info is supplied in TSV in the format parent1/parent2...//junction1/junction2...//minMuts1/minMuts2
    """
    def __init__(self,parentlineages : list[str], junctionpos : list, minMutsForEachPart : list):
        self.parentlineages = parentlineages
        self.junctions = [int(i) for i in junctionpos]
        self.minMutsForEachPart = [int(i) for i in minMutsForEachPart]

    @classmethod
    def getInstance(cls, s : str):
        """factory method to generate ChimeraInfo from string in tsv"""
        p = s.split("//")
        return cls(p[0].split("/"),p[1].split("/"),p[2].split("/") )

class OneVariant:
    """ Holds data for one Sars-Cov2 variant. the list of mutations ..... """
    def __init__(self, name:str, comment:str, mutInfoList :  list[MutInfoSubst], minmutsforpass : int = 0, minstarredmutsforpass : int = 0, parent  = None , parentForCalc  = None, removeFromParentForCalc : bool = False, calcstrategy = None, doPrint : bool = True, childs  = [], histogramGroup = None,
                 histogramOrderId = None, histogramOrderIdWithinBar = None,color : str = None, hatched : bool= False, excempFromChildsSum : bool = False, chimeraInfo : ChimeraInfo = None, showVOCplot : bool =False) :
        self.name = name
        self.comment = comment
        self.minmutsforpass = minmutsforpass
        self.minstarredmutsforpass = minstarredmutsforpass
        self.parent = parent
        self.parentforcalc = parentForCalc # parent used for calc (e.g. BA4/BA5 for BA5)
        self.removeFromParentForCalc = removeFromParentForCalc # indicates whether freq should be removed from parentforcalc(e.g. BA4/BA5 from BA.2
        self.calcstrategy = calcstrategy # calculation strategy for this variant e.g. for BA5 BA.4_BA.5:-:BA.4 means BA4/BA5 minus BA.4
        self.doPrint = doPrint
        self.childs = childs
        self.data = mutInfoList # the MutInfo objects for the mutations
        self.histogramGroup = histogramGroup # string indicating the name of the bar in which variant should be plotted in hist
        self.histogramOrderId = histogramOrderId # indicates bar number this should go in
        self.histogramOrderIdWithinBar = histogramOrderIdWithinBar  # indicates position within bar
        self.color = color
        self.hatched = hatched
        self.excempFromChildsSum = excempFromChildsSum # e.g. BA.4 and BA.5 should not be counted as Omicron childs since BA4_BA5 is already used as Omicron child
        self.naMutList = None # contains list of mutations as string will be initialized at firs call to getNucAcidMutList. NOT USED YET
        self.requiredMutations = None
        self.chimeraInfo = chimeraInfo
        self.showVOCplot = showVOCplot # whether to show special boxplot for this variant

    def addChild(self,c):
        if not c in self.childs:
            self.childs.append(c)

    def setParent(self, parent):
        """adds parent and adds this as child to parent"""
        self.parent = parent
        self.parent.addChild(self)

    def getNucAcidMutList(self):
        """returns list of nuc acid mutations for this variant with mutations as string"""
        if self.naMutList is None: # avoids generating this list everytime the function is called
            self.naMutList = [d.getNAmutstring() for d in self.data]
        return self.naMutList

    def __repr__(self):
        retval = "Node: \n"
        retval += "Name: " + self.name + "\n"
        retval += "Parent: " + ("None" if self.parent is None else str(self.parent)) + "\n"
        retval += "Parent for calc: " + ("None" if self.parentforcalc is None else self.parentforcalc) + "\n"
        retval += "Remove counts from Parent for calc: " + str(self.removeFromParentForCalc) + "\n"
        retval += "Children: "
        retval += ''.join([c.name + " " for c in self.childs])
        retval += "\n"
        retval += "Muts: "
        retval += ''.join([str(d) + " , " for d in self.data])
        return retval

    def isHead(self):
        return self.parent is None

    def testIfMutIsMandatory(self, mutString:str):
        """will test if the supplied mut is mandatary mutation. Will test if supplied str contains the mutstring"""
        if self.data is None:
            return False
        else:
            return any([x.getNAmutstring() in mutString for x in self.data if x.requiredMut])

    def getRequiredMutations(self):
        if self.requiredMutations is None:
            self.requiredMutations = [x.getNAmutstring() for x in self.data if x.requiredMut] if self.data is not None else []
        return self.requiredMutations

    def containsMutation(self, mutInfo : mi.MutInfo) -> bool:
        """returns true if contains the mutation"""
        return mutInfo in self.data
        #return any([x == mutInfo for x in self.data])

    def testWhetherBothPartsOfChimeraOK(self, ser : pandas.Series) -> pandas.Series:
        """will test whether all parts of potential chimeric variant were found and the required minmuts for each part are there
        returns a mask for this series. Big differences in freq between parts leads to masking of part with bigger freq
        returns pos mask for mean"""
        if self.chimeraInfo is None:
            return [True]*len(ser)
        mask = [a == a for a in ser] # NA mask
        ms = ser[mask]
        mutsForRegions = list()
        minPos = 1
        for maxPos in self.chimeraInfo.junctions:
            mutsForRegions.add(ms.loc[lambda x: minPos <= re.search(r'\d+', x).group() <= maxPos ])
            minPos = maxPos + 1
        mutsForRegions.add(ms.loc[lambda x: minPos <= re.search(r'\d+', x).group() <= 40000])
        pass

        #dummy = [re.search(r'\d+', z).group() for z in x.index]

    # def getCSSforColouredDot(self) -> str:
    #     """gets CSS style for coloured dot using color for this variant"""
    #     return '.'+ self.color + """dot {
    #         height:13px;
    #         width:13px;
    #         font-size:13px;
    #         background:""" + self.color + """;
    #         border-radius: 50%;
    #         display:inline-block;
    #         }\n""" if self.color is not None else ""

    def getHTMLforColoredDot(self)-> str:
        """get html string to insert colored dot"""
        #return "<div class=\"" + self.color + "dot\"></div>"
        return "<div class=\"dotForTables\"style=\"background:" + self.color + "\"></div>"

    @classmethod
    def getInstance(cls, name: str, data : pd.DataFrame) :
        """ returns Node object from tsv column"""
        x = data.loc['parent']
        parent = None if x != x else x
        x = data.loc['comment']
        comment = None if x != x else x
        x = data.loc['remove counts from dummy parent']
        removeFromParentForCalc = False if x != x else True
        x = data.loc['dummy parent for calc']
        parentforcalc = None if x != x else x
        x = data.loc['min muts for pass']
        minmutsforpass = 0 if x != x else int(x)
        x = data.loc['min * muts for pass']
        minstarredmutsforpass = 0 if x != x else int(x)
        x = data.loc['calcstrategy']
        calcstrategy = None if x != x else x
        x = data.loc['print']
        doPrint = False if x != x or x != 'TRUE' else True
        x = data.loc['histogram group']
        # histogram info are supplied in the format <histogramGroupName>/<histogramOrderId> - <histogramOrderIdWithinBar>, - <histogramOrderIdWithinBar> is optional
        #just a name for the histogram group
        histogramGroupName = "" if x != x else x.split('/')[0] #[int(v) for v in x.split(',')]
        # the  histogram group, lower will be displayed first
        histogramOrderId = None if x !=x or len(x.split('/'))<2 else int(x.split('/')[1].split('-')[0])
        # the  order within a bar, lower will be displayed first
        histogramOrderIdWithinBar = 9999 if x !=x or len(x.split('/'))<2 or len(x.split('/')[1]) < 2 else int(x.split('/')[1].split('-')[1])
        x = data.loc['color']
        color  = None if x != x else x
        x = data.loc['hatched']
        hatched = False if x != x or x != 'TRUE' else True
        x = data.loc['excempt from child sum']
        excempFromChildsSum = False if x != x or x != 'TRUE' else True
        x = data.loc['chimeraInfo']
        chimeraInfo = None if x != x  else ChimeraInfo.getInstance(x)
        showVOCplot = False # TODO remove this and uncomment two following lines
        #x = data.loc['showVOCplot']
        #showVOCplot = False if x != x or x != 'TRUE' else True
        l = data.loc["mutations":].dropna()
        data = [mi.getMutInfo(x) for x in l if  not x.startswith('#')]                                                       #, (None if len(x.split('_')) < 2 else  x.split('_')[1].split(',')))
           #     for x in l if 'del' not in x and not x.startswith('#')] #  ignore cells starting with # !!!!!!!!! FOR NOW EXCLUDE DELETIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return cls(name,comment,data,minmutsforpass,minstarredmutsforpass,parent,parentforcalc,removeFromParentForCalc,calcstrategy,doPrint,[],histogramGroup =  histogramGroupName, histogramOrderId=histogramOrderId,
                   histogramOrderIdWithinBar = histogramOrderIdWithinBar, color = color, hatched = hatched, excempFromChildsSum = excempFromChildsSum, chimeraInfo = chimeraInfo, showVOCplot = showVOCplot)