#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:        Stitch-seq_Aligner
# Purpose:
#
# Author:      Pengfei
#
# Created:     21/08/2012
# Modified:    03/09/2013
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
from Annotation import *
from cogent.db.ensembl import HostAccount, Genome


def ParseArg():
    p=argparse.ArgumentParser( description = 'Align miRNA-mRNA pairs for Stitch-seq. print the alignable miRNA-mRNA pairs with coordinates', epilog = 'Library dependency: Bio, pysam, itertools')
    p.add_argument('input1',type=str,metavar='part1_reads',help='paired part1 fasta file')
    p.add_argument('input2',type=str,metavar='part2_reads',help='paired part2 fasta file')
    p.add_argument('bowtie_path',type=str,metavar='bowtie_path',help="path for the bowtie program")
    p.add_argument('-b','--bowtie2',action="store_true",help="set to use bowtie2 (--sensitive-local) for alignment")
    p.add_argument('-u','--unique',action="store_true",help="set to only allow unique alignment")
    p.add_argument('-s','--samtool_path',dest='spath', type=str,metavar='samtool_path',help="path for the samtool program",default='samtools')
    p.add_argument('miRNA_ref',type=str,metavar='part1_ref',default="mm9",help="reference genomic seq for part1")
    p.add_argument('mRNA_ref',type=str,metavar='part2_ref',help="reference genomic seq for part2")
    p.add_argument('-a','--annotation',type=str,help='If specified, include the RNA type annotation for each aligned pair, need to give bed annotation RNA file')
    p.add_argument("-A","--annotationGenebed",dest="db_detail",type=str,help="annotation bed12 file for lincRNA and mRNA with intron and exon")


    if len(sys.argv)==1:
        #print (p.print_help(),file=sys.stderr)
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def bowtie_align(b_path,read,ref,s_path,bowtie2):
    # b_path: bowtie path;
    # s_path: samtools path;
    # bowtie2: logic, true/false
    
    sam=read.split("/")[-1].split(".")[0]+".sam"
    if ref.split(".")[-1] in ["fa","fasta"]:
        base=ref.split("/")[-1].split(".")[0]
        os.system("rm "+read+".log")
        os.system(b_path+"-build "+ref+" "+base+" >> "+read+".log 2>&1")
        if not bowtie2:
            os.system(b_path+ " -f -n 1 -l 15 -e 200 -p 9 -S "+base+" "+read+" "+sam+" >> "+read+".log 2>&1")
        else:
            os.system(b_path+ " -x "+base+" -f -U "+read+" --sensitive-local -p 8 --reorder -t -S "+sam+" >> "+read+".log 2>&1")
    else:
        os.system("rm "+read+".log")
        if not bowtie2:
            os.system(b_path+ " -f -n 1 -l 15 -e 200 -p 9 -S "+ref+" "+read+" "+sam+" >> "+read+".log 2>&1")
        else:
            os.system(b_path+ " -x "+ref+" -f -U "+read+" --sensitive-local -p 8 --reorder -t -S "+sam+" >> "+read+".log 2>&1")
    bam=read.split("/")[-1].split(".")[0]+".bam"
    os.system(s_path+ " view -Sb -o "+bam +" "+sam)
    os.system("rm "+sam)
    pysam.sort("-n",bam,"temp")
    align=pysam.Samfile("temp.bam","rb")
    os.system("rm temp.bam")
    os.system(s_path+ " sort "+bam+ " "+"sort_"+read.split("/")[-1].split(".")[0])
    os.system("rm "+bam)
    return align

def Included(record,RequireUnique):
    # record is a pysam read
    # non-unique alignment in Bowtie2 has 'XS' tag: https://www.biostars.org/p/18965/ 
    if RequireUnique:
        try:
            record.opt('XS')
            unique=False
        except:
            unique=True
    else:
        unique=True # not consider unique
    return (not record1.is_unmapped)&unique


def Main():
    args=ParseArg()

    miRNA_align=bowtie_align(args.bowtie_path,args.input1,args.miRNA_ref,args.spath,args.bowtie2)
    mRNA_align=bowtie_align(args.bowtie_path,args.input2,args.mRNA_ref,args.spath,args.bowtie2)
    
    if args.annotation:
        dbi1=DBI.init(args.annotation,"bed")
        dbi2=DBI.init(args.db_detail,"bed")
        dbi3=DBI.init("/home/yu68/bharat-interaction/new_lincRNA_data/mouse.repeat.txt","bed")
          
    for record1, record2 in itertools.izip(miRNA_align, mRNA_align):
        
        if record1.qname.split(" ")[0]!=record2.qname.split(" ")[0]:
            print record1.qname.split(" ")[0]
            print record2.qname.split(" ")[0]
            print >> sys.stderr, "Not match!!"
            sys.exit(0)

        if Included(record1,args.unique) & Included(record2,args.unique):
            strand1="+"
            strand2="+"
            if record1.is_reverse:
                strand1="-"
            if record2.is_reverse:
                strand2="-"
            if args.annotation:
                bed1=Bed([miRNA_align.getrname(record1.tid),record1.pos,record1.aend])
                bed2=Bed([miRNA_align.getrname(record2.tid),record2.pos,record2.aend])
                [name1,typ1,subtype1]=annotation(bed1,dbi1,dbi2,dbi3)
                [name2,typ2,subtype2]=annotation(bed2,dbi1,dbi2,dbi3)
                print '\t'.join(str(f) for f in [miRNA_align.getrname(record1.tid),record1.pos,record1.aend,strand1,record1.seq,name1,typ1,subtype1,record1.qname,mRNA_align.getrname(record2.tid),record2.pos,record2.aend,strand2,record2.seq,name2,typ2,subtype2])
            else:
                print '\t'.join(str(f) for f in [miRNA_align.getrname(record1.tid),record1.aend-record1.alen+1,record1.aend,strand1,record1.seq,record1.qname,mRNA_align.getrname(record2.tid),record2.aend-record2.alen+1,record2.aend,strand2,record2.seq])
    miRNA_align.close()
    mRNA_align.close()

if __name__ == '__main__':
    Main()




