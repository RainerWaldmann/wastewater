import seqtools.GTF as gtf
import variants.globals
import variants.ReadVariants


variants.globals.sarscovgff = gtf.dataframe("./Data/SarsCov2_modified.gff")
with open(r"./Data/sarscov2.fa", 'r') as fp:
    variants.globals.sarscov2seq = fp.readlines()[1].upper()  # read second line
variants.globals.variantdict = variants.ReadVariants.readVariantTSV()
pass
