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



def ParseArg():
    p=argparse.ArgumentParser( description = 'Align miRNA-mRNA pairs for Stitch-seq. print the alignable miRNA-mRNA pairs with coordinates', epilog = 'Library dependency: Bio, pysam, itertools')
    p.add_argument('input1',type=str,metavar='miRNA_reads',help='paired miRNA fastq file')
    p.add_argument('input2',type=str,metavar='mRNA_reads',help='paired mRNA fastq file')
    p.add_argument('bowtie_path',type=str,metavar='bowtie_path',help="path for the bowtie program")
    p.add_argument('-s','--samtool_path',dest='spath', type=str,metavar='samtool_path',help="path for the samtool program",default='samtools')
    p.add_argument('miRNA_ref',type=str,metavar='miRNA_ref',default="hairpin_mouse.fa",help="reference hairpin miRNA seq from miRBase")
    p.add_argument('mRNA_ref',type=str,metavar='mRNA_ref',help="reference genomic seq from mm9")

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
        os.system(b_path+ " -n 1 -l 15 -e 200 -p 6 -q -S "+base+" "+read+" "+sam+" >> "+read+".log 2>&1")
    else:
        os.system("rm "+read+".log")
        os.system(b_path+ " -n 1 -l 15 -e 200 -p 6 -q -S "+ref+" "+read+" "+sam+" >> "+read+".log 2>&1")
    bam=read.split("/")[-1].split(".")[0]+".bam"
    os.system(s_path+ " view -Sb -o "+bam +" "+sam)
    os.system("rm "+sam)
    pysam.sort("-n",bam,"sort_"+read.split("/")[-1].split(".")[0])
    align=pysam.Samfile("sort_"+bam,"rb")
    os.system("rm "+bam)
    return align

def Main():
    args=ParseArg()

    hairpin=SeqIO.parse(args.miRNA_ref,"fasta")
    hairpin_mouse=open("hairpin_mouse.fa","w")
    for record in hairpin:
        if record.id[:3] == 'mmu':
            SeqIO.write(record,hairpin_mouse,"fasta")
    miRNA_align=bowtie_align(args.bowtie_path,args.input1,"hairpin_mouse.fa",args.spath)
    mRNA_align=bowtie_align(args.bowtie_path,args.input2,args.mRNA_ref,args.spath)
    for record1, record2 in itertools.izip(miRNA_align, mRNA_align):
        #print record1.qname, record2.qname
        if ( not record1.is_unmapped) & ( not record2.is_unmapped):
            print miRNA_align.getrname(record1.tid),record1.aend-record1.alen+1,record1.aend,record1.seq
            print mRNA_align.getrname(record2.tid),record2.aend-record2.alen+1,record2.aend,record2.seq
    miRNA_align.close()
    mRNA_align.close()

if __name__ == '__main__':
    Main()







