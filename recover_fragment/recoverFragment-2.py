#-------------------------------------------------------------------------------
# Name:        recoverFragment-2
# Purpose:     classify stitch-seq reads based on fragment lengths. contain wier
#              d read pairs with one end aligned to primer, but not the other
#
# Author:      Pengfei
#
# Created:     25/09/2012
# Updated:     06/08/2013
# Copyright:   (c) Pengfei 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import sys,os,argparse
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from time import time

def ParseArg():
    p=argparse.ArgumentParser( description = 'Recover the real fragment between PCR primers', epilog = 'Library dependency: Bio, itertools')
    p.add_argument('input1',type=str,metavar='reads1',help='forward input fastq/fasta file')
    p.add_argument('input2',type=str,metavar='reads2',help='reverse input fastq/fasta file')
    p.add_argument('primer',type=str,metavar='primer',help='file contianing primer sequences')
    p.add_argument("-v","--verbose",action='store_true',help='specify to print information for each alignment')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()




#=============================================================
# Alignment Parameters
#=============================================================


#=============================================================
# Sequence Alignment
#=============================================================
def local_align(seq1,seq2,Print=False,open_gap=-2,extend_gap=-.3):
    alns=pairwise2.align.localms(seq1,seq2,2,-2,open_gap,extend_gap)
    end_seq1=True
    try:
        top=alns[0]
        score=top[2]
        loc=top[3:5] # start-end point of the second sequence
        if Print:
            print >> sys.stderr,top[0]
            print >> sys.stderr,top[1]
        if (top[0][loc[1]:]==seq1[(loc[1]-loc[0]):]) and (loc[1]-loc[0])!=len(seq1):
            end_seq1=False
    except:
        score=0
        loc=(0,len(seq1))
                 
    return score,loc,end_seq1



'''
type1:
       forward reads:                      XXXX...XXXXNAGATCGGAAGAGCGGTTCAG
                                           ||||...||||
       reverse reads: TGTGCTGCGAGAAGGCTAGANXXXX...XXXX

type2:
       forward reads: XXXXX...XXXXXXXXXXX...XXXX
                                     ||||...||||
       reverse reads:                XXXX...XXXXXXXXXXX...XXXX

type3: 
       forward reads: XXXXXXXXXXXXX

       reverse reads:                      XXXXXXXXXXXXX

'''


def main():
    args=ParseArg()
    #---------------- read pair-end seq ----------------
    print >>sys.stderr,"reading sequence file 1..."
    fastq_iter1 = SeqIO.parse(open(args.input1),"fastq")
    print >>sys.stderr,"reading sequence file 2..."
    fastq_iter2 = SeqIO.parse(open(args.input2),"fastq")
    #---------------- read primer seq ------------------
    primers=[]
    for i in open(args.primer,'r'):
        i=i.strip()
        primers.append(i)
    primer1=primers[0]
    primer2=primers[1]
    #---------------------------------------------------
    name="2_"+os.path.basename(args.input1).split('.')[0][:-3]
    output1=open("short_"+name+".fasta",'w')
    output2=open("long_"+name+".fasta",'w')
    output3=open("evenlong_"+name+"_1.fastq",'w')
    output4=open("evenlong_"+name+"_2.fastq",'w')
    output5=open("wierd_"+name+"_1.fastq",'w')
    output6=open("wierd_"+name+"_2.fastq",'w')

    print >>sys.stderr,"start search..."
    type1=type2=type3=type4=total=0
    t0=time()
    for rec1, rec2 in itertools.izip(fastq_iter1,fastq_iter2):
        seq1=str(rec1.seq)
        seq2=str(rec2.seq)
        seq2_RC=str(rec2.seq.reverse_complement())
        score3,loc3,end_seq1=local_align(seq1,seq2_RC,args.verbose)
       
        if (score3>=1.8*(loc3[1]-loc3[0])):
            if end_seq1:
                fragment=SeqRecord(rec1.seq+rec2.seq.reverse_complement()[(loc3[1]-loc3[0]):],id=rec1.id,name=rec1.name,description=rec1.description)
                SeqIO.write(fragment,output2,"fasta")
                if args.verbose:
                    print >> sys.stderr,fragment.seq
                    print >> sys.stderr, "type2\n"
                type2+=1
            else:
                score1,loc1,end1=local_align(seq1[(loc3[1]-loc3[0]):],primer1,args.verbose,-.5)
                score2,loc2,end2=local_align(seq2[(loc3[1]-loc3[0]):],primer2,args.verbose,-.5)
                if (score1+score2>=1.5*(loc1[1]-loc1[0]+loc2[1]-loc2[0])): # use a less stringent cut for the primer match
                    fragment=rec1[:(loc3[1]-loc3[0])]
                    SeqIO.write(fragment,output1,"fasta")
                    type1+=1
                    if args.verbose:
                        print >> sys.stderr,fragment.seq
                        print >> sys.stderr, "type1\n"
                else:
                    SeqIO.write(rec1,output5,'fastq')
                    SeqIO.write(rec2,output6,'fastq')
                    type4+=1
                    #score1,loc1,end1=local_align(seq1[(loc3[1]-loc3[0]):],primer1,True)
                    #score2,loc2,end2=local_align(seq2[(loc3[1]-loc3[0]):],primer2,True)
                    if args.verbose:
                        print >> sys.stderr, "type4\n"
                        print >> sys.stderr,score1,loc1[1]-loc1[0],score2,loc2[1]-loc2[0]
        else:
            score4,loc4,end4=local_align(seq1[2*len(seq1)/3:],seq2_RC[:len(seq2)/3],args.verbose,-5)
            if score4>=1.7*max(loc4[1]-loc4[0],7):
                fragment=SeqRecord(rec1.seq+rec2.seq.reverse_complement()[(loc4[1]-loc4[0]):],id=rec1.id,name=rec1.name,description=rec1.description)
                SeqIO.write(fragment,output2,"fasta")
                if args.verbose:
                    print >> sys.stderr,score4,loc4
                    print >> sys.stderr,fragment.seq
                    print >> sys.stderr, "type2\n"
                type2+=1
            else:
                SeqIO.write(rec1,output3,'fastq')
                SeqIO.write(rec2,output4,'fastq')
                if args.verbose:
                    print >> sys.stderr,score4,loc4
                    print >> sys.stderr, "type3\n"
                type3+=1

        total+=1
        if total%1000==0:
            t1=time()
            print >> sys.stderr, " %d read pairs processed and %d type1(short), %d type2(long), %d type3(evenlong), and %d type4(wierd). Time:%.2fmin\r "%(total, type1, type2, type3,type4,(t1-t0)/60),
            t0=time()
    output1.close()
    output2.close()
    output3.close()
    output4.close()
    output5.close()
    output6.close()

if __name__ == '__main__':
    main()
