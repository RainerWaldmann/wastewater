from itertools import groupby
import variants.globals
import pandas as pd

codonTable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

reverseTranslationTable = {'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                           'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'],
                           'E': ['GAA', 'GAG'], 'D': ['GAT', 'GAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                           'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'N': ['AAT', 'AAC'], 'M': ['ATG'], 'K': ['AAA', 'AAG'],
                           'Y': ['TAT', 'TAC'], 'I': ['ATT', 'ATC', 'ATA'], 'Q': ['CAA', 'CAG'], 'F': ['TTT', 'TTC'],
                           'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
                           '*': ['TAA', 'TAG', 'TGA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                           'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC']}


def translateprot(seq):
    """Translates string of NAs"""
    if not seq:
        return ""
    protein = ""
    seq = seq.upper()
    if len(seq) % 3 != 0:
        print("SEQ is:" + seq)
        raise Exception(" Translated sequence length % 3 must be 0")
    else:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += codonTable[codon]
    return protein

def translateAA(seq : str):
    """ Translates one codon"""
    aa = None
    if len(seq.upper()) == 3 :
        aa = codonTable[seq]
    return aa



def getNucMutsForAAmutsFile(inFile):
    """takes TSV with AA muts in columns and creates TSV with Nuc muts
    data is in format gene:POS:AA"""
    p = pd.read_table(inFile)
    result = {}
    for(colname,colval) in p.iteritems():
        oneVarRslt = []
        for v in colval:
            if v == v:
                s = v.split(':')
                r = getPossibleNucMutsForOneAAmut(s[0], s[1], s[2])
                #r = [x + "/" + s[0] + ':' + str(s[1]) + s[2] for x in m]
                oneVarRslt += r
        result[colname] = oneVarRslt
    #https://stackoverflow.com/questions/19736080/creating-dataframe-from-a-dictionary-where-entries-have-different-lengths
    new_df = pd.DataFrame({key: pd.Series(value) for key, value in result.items()})
    new_df.to_csv(inFile.rsplit( ".", 1 )[ 0 ]  + "_NAmuts.tsv", sep="\t")



def getPossibleNucMutsForOneAAmut(gene, pos, newAminoAcid):
    """takes gene pos of AA in gene and returns list of NA muts format: OldaaPosNewaa/OldNucPosNewnuc"""
    start = None
    for i, row in variants.globals.sarscovgff.iterrows():
        if row['ID'] == gene:
            start = int(row.start)
            break
    if start is not None :
        startOfCodon = int(pos) * 3 - 2 # yields 1,4,7 for 1,2,3 -> start of codon in ORF; 1-based
        startOfCodonOnSars = start + startOfCodon - 1 # start of codon in SarsCov-2, 1 - based
        wt_codon = variants.globals.sarscov2seq[startOfCodonOnSars-1:startOfCodonOnSars+2]# the wild type codon
        wtAA = codonTable[wt_codon]

        poss_codons = reverseTranslationTable[newAminoAcid]  # list with all possible codons for new amino acid
        # find the codon that requires the least mutations
        codons_muts = []  #list of tuples with n mutations , codon
        for mut_codon in poss_codons:
            mm = sum([mut_codon[i] != wt_codon[i] for i in range(3)])  # count number of mutations required
            if mm != 0:
                codons_muts.append((mm, mut_codon))
        codons_muts.sort(key=lambda x: x[0])  # sort by number of muts lowest first
        g = groupby(codons_muts, lambda x: x[0])  # get the first group (least mutations)
        # get list of codons that require the least changes
        x= next(g) # gets the first element -> the one with the least mutations, x is here a tuple with nmuts,itertool.grouper object with data
        x = list(x[1]) # transforms itertool.grouper object into list, still list of tuples with nmuts,codon
        possMutCodons = [v[1] for v in x] # extract the codons
        # other option would be to transform groupby output into dictionary ; groups = {key: list(v) for (key, v) in g}
        # get list of muts in format NucPosNuc
        mutStrings = []
        # for accepted codon changes get list of mutations
        for cod in possMutCodons:
            for i in range(len(cod)):
                if cod[i] != wt_codon[i]:
                    mutStrings.append(str(wt_codon[i]) + str(startOfCodonOnSars + i) + str(cod[i]))
        mutStrings = [m + "/" + gene + ':' + wtAA + str(pos) + newAminoAcid for m in mutStrings]
    return mutStrings