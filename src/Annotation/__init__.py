"""
For the purpose of annotating RNA types for genomic regions.
"""


#from xplib import DBI
#from cogent.db.ensembl import HostAccount, Genome

def overlap(bed1,bed2):
    """
    This function compares overlap of two Bed object from same chromosome
    
    :param bed1: A Bed object from `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :param bed2: A Bed object from `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :returns: boolean -- True or False

    Example:
    
    >>> from xplib.Annotation import Bed
    >>> from Annotation import overlap
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
    
    :param bed1: A Bed object defined by `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (BAM2X)
    :param genebed: A Bed12 object representing a transcript defined by xplib Annotaton.Bed with information of exon/intron/utr from an BED12 file
    :returns: str -- RNA subtype. "intron"/"exon"/"utr3"/"utr5"/"."

    Example:
    
    >>> from xplib.Annotation import Bed
    >>> from xplib import DBI
    >>> from Annotation import Subtype
    >>> bed1=Bed(["chr13",40975747,40975770])
    >>> a=DBI.init("../../Data/Ensembl_mm9.genebed.gz","bed")
    >>> genebed=a.query(bed1).next()
    >>> print Subtype(bed1,genebed)
    "intron"
        
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


def annotation(bed,ref_allRNA,ref_detail,ref_repeat):
    """
    This function is based on :func:`overlap` and :func:`Subtype` functions to annotate RNA type/name/subtype for any genomic region.

    :param bed: A Bed object defined by `xplib.Annotation.Bed <http://bam2xwiki.appspot.com/bed>`_ (in BAM2X).
    :param ref_allRNA: the `DBI.init <http://bam2xwiki.appspot.com/DBI>`_ object (from BAM2X) for bed6 file of all kinds of RNA
    :param ref_detail: the `DBI.init <http://bam2xwiki.appspot.com/DBI>`_ object for bed12 file of lincRNA and mRNA with intron, exon, UTR
    :param ref_detail: the `DBI.init <http://bam2xwiki.appspot.com/DBI>`_ object for bed6 file of mouse repeat
    :returns: list of str -- [type,name,subtype]

    Example:

    >>> from xplib.Annotation import Bed
    >>> from xplib import DBI
    >>> from Annotation import annotation
    >>> bed=Bed(["chr13",40975747,40975770])
    >>> ref_allRNA=DBI.init("../../Data/all_RNAs-rRNA_repeat.txt.gz","bed")
    >>> ref_detail=DBI.init("../../Data/Ensembl_mm9.genebed.gz","bed")
    >>> ref_repeat=DBI.init("../../Data/mouse.repeat.txt.gz","bed")
    >>> print annotation(bed,ref_allRNA,ref_detail,ref_repeat)
    ["protein_coding","gcnt2","intron"]

    """
    flag=0
    typ="non"
    name="."
    subtype="."
    max_overlap = 0  # find annotation with largest overlap
    for hit in ref_allRNA.query(bed):
        overlap = min(hit.stop,bed.stop)-max(hit.start,bed.start)
        if overlap>max_overlap:
            name=hit.id.split(".")[1]
            typ=hit.id.split(".")[0]
            max_overlap = overlap
        if overlap == hit.stop-hit.start and typ=="snoRNA":  #annotated as snoRNA if region covers whole snoRNA
            break
    if (typ=="lincRNA" or typ=="protein_coding"):
        flag=0
        for hit in ref_detail.query(bed):
            if flag==0:
                tempname=hit.id.split("&")
                subtype=Subtype(bed,hit)
                try:
                  typ = tempname[1]
                except:
                  pass
                name = tempname[0]
                if subtype!="intron":
                    flag=1
        '''
        try:
            tran=genome.getTranscriptByStableId(StableId=tempname).Gene
            typ=tran.BioType
            name=tran.Symbol
        except: pass
        '''
    if typ=="non":
        for hit in ref_repeat.query(bed):
            tempname=hit.id.split("&")
            name = tempname[0]
            typ = tempname[1]
            subtype = tempname[2]
            break
        '''
        try:
            repeats=genome.getFeatures(CoordName=bed.chr[3:], Start=bed.start, End=bed.stop, feature_types='repeat')
            for r in repeats:
               if r.RepeatClass!='dust':
                   typ=r.RepeatType
                   name=r.Symbol
                   break 
        except: pass
        '''
    if typ=="lincRNA" and subtype!="intron":
        subtype="utr"
    if typ in ["snoRNA","snRNA","miRNA","miscRNA"]:
        subtype='.'
    if typ=="pseudogene":
        subtype="."
    return [typ,name,subtype]
                    
           
