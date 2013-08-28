#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:        Stitch-seq_Aligner
# Purpose:
#
# Author:      Pengfei
#
# Created:     21/08/2012
# Copyright:   (c) Pengfei 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import sys, os, argparse
import pysam
import itertools
from Bio import SeqIO
from xplib import TableIO
from xplib import DBI
from xplib.Annotation import Bed


def ParseArg():
    p=argparse.ArgumentParser( description = 'Align miRNA-mRNA pairs for Stitch-seq. print the alignable miRNA-mRNA pairs with coordinates', epilog = 'Library dependency: Bio, pysam, itertools')
    p.add_argument('input1',type=str,metavar='miRNA_reads',help='paired miRNA fastq file')
    p.add_argument('input2',type=str,metavar='mRNA_reads',help='paired mRNA fastq file')
    p.add_argument('bowtie_path',type=str,metavar='bowtie_path',help="path for the bowtie program")
    p.add_argument('-s','--samtool_path',dest='spath', type=str,metavar='samtool_path',help="path for the samtool program",default='samtools')
    p.add_argument('miRNA_ref',type=str,metavar='miRNA_ref',default="hairpin_mouse.fa",help="reference hairpin miRNA seq from miRBase")
    p.add_argument('mRNA_ref',type=str,metavar='mRNA_ref',help="reference genomic seq from mm9")
    p.add_argument('-a','--annotation',type=str,help='If specified, include the RNA type annotation for each aligned pair, need to give bed annotation RNA file')


    if len(sys.argv)==1:
        #print (p.print_help(),file=sys.stderr)
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def bowtie_align(b_path,read,ref,s_path):
    # b_path: bowtie path;
    # s_path: samtools path;
    sam=read.split("/")[-1].split(".")[0]+".sam"
    if ref.split(".")[-1] in ["fa","fasta"]:
        base=ref.split("/")[-1].split(".")[0]
        os.system("rm "+read+".log")
        os.system(b_path+"-build "+ref+" "+base+" >> "+read+".log 2>&1")
        os.system(b_path+ " -f -n 1 -l 15 -e 200 -p 6 -S "+base+" "+read+" "+sam+" >> "+read+".log 2>&1")
    else:
        os.system("rm "+read+".log")
        os.system(b_path+ " -f -n 1 -l 15 -e 200 -p 6 -S "+ref+" "+read+" "+sam+" >> "+read+".log 2>&1")
    bam=read.split("/")[-1].split(".")[0]+".bam"
    os.system(s_path+ " view -Sb -o "+bam +" "+sam)
    os.system("rm "+sam)
    pysam.sort("-n",bam,"sort_"+read.split("/")[-1].split(".")[0])
    align=pysam.Samfile("sort_"+bam,"rb")
    os.system("rm "+bam)
    return align

def annotate(bed,dbi):
    flag=0
    typ="other"
    name="."
    for hit in dbi.query(bed):
        if flag==0:
            name=hit.id.split(".")[1]
            typ=hit.id.split(".")[0]
            flag=1
    return [name,typ]


def Main():
    args=ParseArg()

    miRNA_align=bowtie_align(args.bowtie_path,args.input1,args.miRNA_ref,args.spath)
    mRNA_align=bowtie_align(args.bowtie_path,args.input2,args.mRNA_ref,args.spath)
    
    if args.annotation:
        dbi=DBI.init(args.annotation,"bed")

    for record1, record2 in itertools.izip(miRNA_align, mRNA_align):
        #print record1.qname, record2.qname
        if ( not record1.is_unmapped) & ( not record2.is_unmapped):
            if args.annotation:
                bed1=Bed([miRNA_align.getrname(record1.tid),record1.pos,record1.aend])
                bed2=Bed([miRNA_align.getrname(record2.tid),record2.pos,record2.aend])
                [name1,typ1]=annotate(bed1,dbi)
                [name2,typ2]=annotate(bed2,dbi)
                print '\t'.join(str(f) for f in [miRNA_align.getrname(record1.tid),record1.pos,record1.aend,record1.seq,name1,typ1,record1.qname,mRNA_align.getrname(record2.tid),record2.pos,record2.aend,record2.seq,name2,typ2])
            else:
                print '\t'.join(str(f) for f in [miRNA_align.getrname(record1.tid),record1.aend-record1.alen+1,record1.aend,record1.seq,record1.qname,mRNA_align.getrname(record2.tid),record2.aend-record2.alen+1,record2.aend,record2.seq])
    miRNA_align.close()
    mRNA_align.close()

if __name__ == '__main__':
    Main()







