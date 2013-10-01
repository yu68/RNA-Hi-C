#!/usr/bin/env python
# Pengfei Yu yupf05@gmail.com

import argparse,sys
from xplib import TableIO
from xplib import DBI
from Annotation import *
from cogent.db.ensembl import HostAccount, Genome

def ParseArg():
    p=argparse.ArgumentParser(description = "RNA composition for aligned sample using ensembl annotation",epilog="Library dependency: bam2x")
    p.add_argument("-i","--input",type=str,help="input aligned file")
    p.add_argument("-f","--format",type=str,default="bam",help="input file format, default: bam")
    p.add_argument("-o","--output",type=str,default="stdout",help="output txt file")
    p.add_argument("-a","--annotation",dest="db",type=str,help="annotation bed6 file for all kinds of RNA")
    p.add_argument("-A","--annotationGenebed",dest="db_detail",type=str,help="annotation bed12 file for lincRNA and mRNA with intron and exon")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def Main():
    global args,out
    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",args.output,"to write. Using stdout instead"
            out=sys.stdout

    count={}
    dbi1=DBI.init(args.db,"bed") # the DBI init file for bed6 file of all kinds of RNA
    dbi2=DBI.init(args.db_detail,"bed") # the DBI init file for bed12 file of lincRNA and mRNA with intron, exon, UTR
    genome=Genome('mouse', Release=67, account=None)
    for bed in TableIO.parse(args.input,args.format):
        [typ,name,subtype]=annotation(bed,dbi1,dbi2,genome)
        if count.has_key(typ):
            count[typ]+=1
        else:
            count[typ]=1
        print >>out, "\t".join (str(f) for f in [bed.chr,bed.start,bed.stop,bed.id,name,bed.strand,typ, subtype])

    print >>out, "\n".join ("#"+typ+"\t%d"%(count[typ]) for typ in count.keys())


if __name__=="__main__":
    Main()

