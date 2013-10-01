"""
For the purpose of annotating RNA types for genomic regions.
"""


#from xplib import DBI
#from cogent.db.ensembl import HostAccount, Genome

def overlap(bed1,bed2):
    """
    This function compares overlap of two Bed object from same chromosome
    
    :param bed1: A Bed object from xplib.Annotation.Bed (BAM2X)
    :param bed2: A Bed object from xplib.Annotation.Bed (BAM2X)
    :returns: boolean -- True or False

    Example:
    
    >>> from xplib.Annotation import Bed
    >>> bed1=Bed(["chr1",10000,12000])
    >>> bed2=Bed(["chr1",9000,13000])
    >>> print overlap(bed1,bed2)
    True

    """
    try:
        return (bed1.stop>bed2.start) and (bed1.start<bed2.stop)
    except: # in case for "NonType" of bed2
        return False

def Subtype(bed1,genebed):
    """
    This function determines intron or exon or utr from a BED12 file.
    
    :param bed1: A Bed object defined by xplib.Annotation.Bed (BAM2X)
    :param genebed: A Bed12 object representing a transcript defined by xplib Annotaton.Bed with information of exon/intron/utr from an BED12 file
    :returns: str -- RNA subtype. "intron"/"exon"/"utr3"/"utr5"/"."

    Example:
    
    >>> from xplib.Annotation import Bed
    >>> from xplib import DBI
    >>> bed1=Bed(["chr13",40975747,40975770])
    >>> a=DBI.init("../../Data/Ensembl_mm9.genebed.gz","bed")
    >>> genebed=a.query(bed1).next()
    >>> print Subtype(bed1,genebed)
    "Intron"
        
    """
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


def annotation(bed,ref_allRNA,ref_detail,genome):
    """
    This function is based on :func:`overlap` and :func:`Subtype` functions to annotate RNA type/name/subtype for any genomic region.

    :param bed: A Bed object defined by xplib.Annotation.Bed (in BAM2X).
    :param ref_allRNA: the DBI.init object (from BAM2X) for bed6 file of all kinds of RNA
    :param ref_detail: the DBI.init object for bed12 file of lincRNA and mRNA with intron, exon, UTR
    :param genome: the Genome object from PyCogent.
    :returns: list of str -- [type,name,subtype]

    Example:

    >>> from xplib.Annotation import Bed
    >>> from xplib import DBI
    >>> from cogent.db.ensembl import Genome
    >>> bed=Bed(["chr13",40975747,40975770])
    >>> ref_allRNA=DBI.init("../../Data/all_RNAs-rRNA_repeat.txt.gz","bed")
    >>> ref_detail=DBI.init("../../Data/Ensembl_mm9.genebed.gz","bed")
    >>> genome=Genome('mouse', Release=67, account=None)
    >>> print annotation(bed,ref_allRNA,ref_detail,genome)
    ["protein_coding","gcnt2","intron"]

    """
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
                tempname=hit.id
                subtype=Subtype(bed,hit)
                if subtype!="intron":
                    flag=1
        try:
            tran=genome.getTranscriptByStableId(StableId=tempname).Gene
            typ=tran.BioType
            name=tran.Symbol
        except: pass
    if typ=="non":
        try:
            repeats=genome.getFeatures(CoordName=bed.chr[3:], Start=bed.start, End=bed.stop, feature_types='repeat')
            for r in repeats:
               if r.RepeatClass!='dust':
                   typ=r.RepeatType
                   name=r.Symbol
                   break 
        except: pass
    return [typ,name,subtype]
                    
           
