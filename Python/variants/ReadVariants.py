import pandas as pd
import variants.OneVariant as OneVariant
import variants.globals

def readVariantTSV() -> dict():
    """returns dictionary of Nodes"""
    p = pd.read_table("./Data/Variants.tsv.txt", index_col=0, skip_blank_lines=True)
    p.dropna(how='all', axis='columns',inplace=True) # drop columns that only contain NaNs
    retval = dict()
    for column in p.columns :
        if not column.startswith('#'): #ignore columns where first cell starts with '#'
            x = OneVariant.OneVariant.getInstance(column, p[column])
            retval[column] = x
    for z in retval.values():
        if z.parent:
            z.setParent(retval[z.parent]) # parent was just string before -> set to node, will also set z as child for parent
        if z.parentforcalc :
            z.parentforcalc =retval[z.parentforcalc]
    for z in retval.values(): # iterate over variants
        for y in retval.values(): # iterate over mutations
            if y is not z: # avoid self comparison
                for mut in z.data:
                    if y.containsMutation(mut):
                        mut.addOtherNodeWhereMutIsFound(y)

    return retval

