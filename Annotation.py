from xplib import DBI

def overlap(bed1,bed2):
    '''
    compare overlap for same chromosome
    '''
    try:
        return (bed1.stop>bed2.start) and (bed1.start<bed2.stop)
    except: # in case for "NonType" of bed2
        return False

def Subtype(bed1,genebed):
    '''
    determine intron or exon or utr
    '''
    subtype="intron"
    if overlap(bed1,genebed.utr3()):
        for i in genebed.utr3().Exons():
            if overlap(bed1,i):
                subtype="utr3"
    elif overlap(bed1,genebed.utr5()):
        for i in genebed.utr5().Exons():
            if overlap(bed1,i):
                subtype="utr5"
    else:
        for i in genebed.Exons():
            if overlap(bed1,i):
                subtype="exon"
                break
    return subtype
    
    
def annotation(bed,ref_allRNA,ref_detail):
    '''
    ref_allRNA is the DBI init file for bed6 file of all kinds of RNA
    ref_detail is the DBI init file for bed12 file of lincRNA and mRNA with intron, exon, UTR 
    '''
    
    flag=0
    typ="non"
    name="."
    subtype="."
    for hit in ref_allRNA.query(bed):
        if flag==0:
            name=hit.id.split(".")[1]
            typ=hit.id.split(".")[0]
            flag=1
    if (typ=="lincRNA" or typ=="protein_coding"):
        flag=0
        for hit in ref_detail.query(bed):
            if flag==0:
                subtype=Subtype(bed,hit)
                flag=1
    return [typ,name,subtype]
                    
           