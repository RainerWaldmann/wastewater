import variants.globals
import math
import seqtools.CodonTables as ct
import re
from enum import Enum
import variants.OneVariant as OneVariant


regexformutsplitSubst = re.compile("([a-zA-Z]+)([0-9]+)([a-zA-Z]+)")
regexformutsplitDel = re.compile("([0-9]+)(del)([0-9]+)", re.IGNORECASE)
#regexformutsplitDel = re.compile("([0-9]+)((?i)del)([0-9]+)")

class MutInfo:
    """Just a base class that holds common info on mutations"""
    class MutType(Enum):
        INSERTION = 1
        DELETION = 2
        SUBSTITUTION = 3


    def __init__(self, mt : MutType, position: int,  wuhanseq: str, mutseq: str, gene: str, requiredMut:bool):
        self.mutType = mt
        self.position = position  # position in sarscov
        self.wuhanseq = wuhanseq  # wt seq for substitution
        self.mutseq = mutseq  # mutated seq
        self.gene = gene
        self.requiredMut = requiredMut
        self.listOfOtherVariantsMutIsFound : list[OneVariant.OneVariant] = None # can hold list of variants that have this mutation

    def getMutType(self)  -> MutType:
        return self.mutType

    def __eq__(self,other):
        """simply tests whether pos and mutseq are equal"""
        return self.position == other.position and self.mutseq == other.mutseq

    def addOtherNodeWhereMutIsFound(self, nde ):
        """adds Node to list of nodes where mut is found"""
        if self.listOfOtherVariantsMutIsFound is None:
            self.listOfOtherVariantsMutIsFound = [nde]
        elif nde not in self.listOfOtherVariantsMutIsFound:
            self.listOfOtherVariantsMutIsFound.append(nde)

    def getAAmutstringForDFindex(self) -> str:
        """return string such as _spike:K452L or empty string if None. Used for index of dataframe"""
        s = self.getAAmutstring()
        aachange = s.split(':')[1] if s is not None else ""
        return "_" + s if s is not None and (aachange[0] != aachange[
            -1] or  isinstance(self, MutInfoDeletion)) else ""  # aachange[0] != aachange[-1] to avoid adding AA mutation for silent nuc mutations


###################################################################################################################################
class MutInfoDeletion(MutInfo):
    """holds info on codon and mutation for one deletion
    mutseq is e.g. del9 for 9 nt deletion"""
    def __init__(self, position: int,   delLength: int, wuhancodons: str, mutcodons: str, gene: str,
                 aaPosInGene: int,requiredMut:bool = False):
        super().__init__(MutInfo.MutType.DELETION,position, "", "DEL" + str(delLength), gene, requiredMut)
        self.wuhancodons = wuhancodons  # wuhan codon
        self.mutcodons = mutcodons  # mutated codon
        self.aaPosInGene = aaPosInGene  # 1 - based
        self.delLength = delLength
        self.oldpeptide = None # initialized on request
        self.newpeptide = None # initialized on request

    def __repr__(self):
        return str(self.position) + 'del' + str(self.delLength)

    def getNAmutstring(self) -> str:
        """return string such as C2000T"""
        return self.__repr__()

    def getAAmutstring(self) -> str:
        """return string such as spike:KM452RN. Pos is pos of first AA"""
        retval = None
        if self.wuhancodons is not None and self.mutcodons is not None:
            if self.oldpeptide is None:
                self.oldpeptide = ct.translateprot(self.wuhancodons)
            if self.newpeptide is None:
                if not self.mutcodons: #empty string
                    self.newpeptide = 'del' # add del instead of empty string to print e.g. V32 del instead of V32 for mutstring in heatmap....
                else:
                    self.newpeptide = ct.translateprot(self.mutcodons)
            retval =  self.gene + ':' + self.oldpeptide + str(self.aaPosInGene) + \
               self.newpeptide
        return retval

    @classmethod
    def getMutInfoFromPosMut(cls, pos, delStr:int, required=False) -> MutInfo:
        """accepts string with multiple nucleotides"""
        retval = None
        delLength = int(re.findall('\d+', delStr)[0])
        matchingGffData = None
        for row in variants.globals.sarscovgff.itertuples():
            if pos >= int(row.start) and pos <= int(row.end):
                matchingGffData = row
                break
        if matchingGffData is not None:
            cds = variants.globals.sarscov2seq[int(matchingGffData.start) - 1: int(matchingGffData.end)]
            codonnumberDelStart = math.floor((pos  - int(matchingGffData.start)) / 3)  # 0 -based
            codonnumberEndBeforeDel = math.floor((pos + delLength -1 - int(matchingGffData.start)) / 3)
            offsetDelOnCodon = pos  - int(matchingGffData.start) - codonnumberDelStart * 3
            oldCodons = cds[codonnumberDelStart * 3:codonnumberEndBeforeDel * 3 + 3]
            newCodons = "" if offsetDelOnCodon == 0 else oldCodons[:offsetDelOnCodon] + oldCodons[-3 + offsetDelOnCodon :]
            #x = len(oldCodons) - delLength
            #newCodons = oldCodons[0: x + (0 if x % 3 == 0 else 3 - x % 3)]
            #newCodons = oldCodons[0: len(oldCodons) - delLength + (0 if x % 3 == 0 else 3 - x % 3)] --> good example that Python sucks !!!!!!!!!!!
            retval = cls(pos, delLength,oldCodons,newCodons,matchingGffData.ID,codonnumberDelStart + 1, required) # codonnumberDelStart + 1 : O-based -> 1-based
        else:
            retval = cls(pos, delLength, None, None, None,None,required)
        return retval
#################################################################################################################################
class MutInfoSubst(MutInfo):
    """holds info on codon and mutation for one variation"""

    def __init__(self, position: int,  wuhanseq: str, mutseq: str, wuhancodon: str, mutcodon: str, gene: str,
                 aaPosInGene: int,requiredMut:bool = False, listOfVariantsMutIsFound : list[str] = None):
        super().__init__(MutInfo.MutType.SUBSTITUTION,position, wuhanseq, mutseq, gene, requiredMut)
        self.wuhancodon = wuhancodon  # wuhan codon
        self.mutcodon = mutcodon  # mutated codon
        self.aaPosInGene = aaPosInGene  # 1 - based

    def __repr__(self):
        return str(self.wuhanseq) + str(self.position) + str(self.mutseq)

    def getNAmutstring(self) -> str:
        """return string such as C2000T"""
        return self.wuhanseq.upper() + str(self.position) + self.mutseq.upper()

    def getAAmutstring(self) -> str:
        """return string such as spike:K452L"""
        return self.gene + ':' + ct.codonTable[self.wuhancodon] + str(self.aaPosInGene) + \
               ct.codonTable[self.mutcodon] if self.mutcodon else None

    @classmethod
    def getMutInfoFromPosMut(cls,pos,  newnuc,  required = False) -> MutInfo :
        """accepts string with multiple nucleotides"""
        retval = None
        matchingGffData = None
        for row in variants.globals.sarscovgff.itertuples():
            if pos >= int(row.start) and pos <= int(row.end):
                matchingGffData = row
                break
        if matchingGffData is not None:
            cds = variants.globals.sarscov2seq[int(matchingGffData.start) - 1: int(matchingGffData.end)]
            codonnumber = math.floor((pos - int(matchingGffData.start)) / 3)  # 0 -based
            codon = cds[codonnumber * 3:codonnumber * 3 + 3]
            posincodon = (pos - int(matchingGffData.start)) % 3  # 0-based
            newcodon = list(codon)  # because str is imutable
            for i in range(len(newnuc)):
                newcodon[posincodon + i] = newnuc[i]
            retval = cls(pos, codon[posincodon:posincodon + len(newnuc)], newnuc,
                         codon.upper(), "".join(newcodon).upper(), matchingGffData.ID, codonnumber + 1, required)
        else:
            retval = cls(pos, variants.globals.sarscov2seq[pos - 1: pos - 1 + len(newnuc)],
                         newnuc, None, None, None, required)
        return retval



def getMutInfoFromPosMut(pos,  newnuc,  required = False) -> MutInfo :
        """accepts string with multiple nucleotides"""
        retval = None
        matchingGffData = None
        for row in variants.globals.sarscovgff.itertuples():
            if pos >= int(row.start) and pos <= int(row.end):
                matchingGffData = row
                break
        if matchingGffData is not None:
            cds = variants.globals.sarscov2seq[int(matchingGffData.start) - 1: int(matchingGffData.end)]
            codonnumber = math.floor((pos - int(matchingGffData.start)) / 3)  # 0 -based
            codon = cds[codonnumber * 3:codonnumber * 3 + 3]
            posincodon = (pos - int(matchingGffData.start)) % 3  # 0-based
            newcodon = list(codon)  # because str is imutable
            for i in range(len(newnuc)):
                newcodon[posincodon + i] = newnuc[i]
            retval = MutInfoSubst(pos, codon[posincodon:posincodon + len(newnuc)], newnuc,
                         codon.upper(), "".join(newcodon).upper(), matchingGffData.ID, codonnumber + 1, required)
        else:
            retval = MutInfoSubst(pos, variants.globals.sarscov2seq[pos - 1: pos - 1 + len(newnuc)],
                         newnuc, None, None, None, required)
        return retval



def getMutInfo (s:str):
    required = s[0] == '*'
    x = s[1:].split('/')[0] if s.startswith('*') else s.split('/')[0] # get REF POS ALT, kept the '/' split in case I add something behind a '/'
    d = regexformutsplitSubst.match(x.strip())
    if d is not None:
        z = d.groups()
        return MutInfoSubst.getMutInfoFromPosMut(int(z[1]), z[2], required)
    else:
        d = regexformutsplitDel.match(x.strip())
        if d is not None:
            z = d.groups()
            r = MutInfoDeletion.getMutInfoFromPosMut(int(z[0]), 'del' + z[2], required)
            return r
        else:
            raise Exception("argument does not have nucsposnucs format for subst or posDELlength for deletion")





