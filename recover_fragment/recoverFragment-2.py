#-------------------------------------------------------------------------------
# Name:        recoverFragment
# Purpose:
#
# Author:      Pengfei
#
# Created:     25/09/2012
# Copyright:   (c) Pengfei 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import sys,os,argparse
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def ParseArg():
    p=argparse.ArgumentParser( description = 'Recover the real fragment between PCR primers', epilog = 'Library dependency: Bio, itertools')
    p.add_argument('input1',type=str,metavar='reads1',help='forward input fastq/fasta file')
    p.add_argument('input2',type=str,metavar='reads2',help='reverse input fastq/fasta file')
    p.add_argument('primer',type=str,metavar='primer',help='file contianing primer sequences')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()


def make_matrix(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros."""
    return [[0]*sizey for i in range(sizex)]


def is_complement(a, b):
    """Return True if character a is complmentary to character b"""
    assert len(a) == len(b) == 1
    return (a.upper(), b.upper()) in [
        ("A", "T"), ("T", "A"),
        ("C", "G"), ("G", "C"),
        ("A", "U"), ("U", "A")
    ]


#=============================================================
# Alignment Parameters
#=============================================================

class ScoreParam:
    """Stores the parameters for an alignment scoring function"""
    def __init__(self, match, mismatch, gap, gap_start=0):
        self.gap_start = gap_start
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

    def matchchar(self, a,b):
        """Return the score for aligning character a with b"""
        assert len(a) == len(b) == 1
        if a==b:
            return self.match
        else:
            return self.mismatch

    def __str__(self):
        return "match = %d; mismatch = %d; gap_start = %d; gap_extend = %d" % (
                self.match, self.mismatch, self.gap_start, self.gap
        )

#=============================================================
# Sequence Alignment
#=============================================================

def local_align(x, y, score=ScoreParam(10, -8, -10)):
    """Do a local alignment between x and y with the given scoring parameters.
    We assume we are MAXIMIZING."""

    # create a zero-filled matrix
    A = make_matrix(len(x) + 1, len(y) + 1)

    best = 0
    optloc = (0,0)

    # fill in A in the right order
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):

            # the local alignment recurrance rule:
            A[i][j] = max(
               A[i][j-1] + score.gap,
               A[i-1][j] + score.gap,
               A[i-1][j-1] + score.matchchar(x[i-1], y[j-1]),
               0
            )

            # track the cell with the largest score
            if A[i][j] >= best:
                best = A[i][j]
                optloc = (i,j)
    ''' 
    print "Scoring:", str(score)
    #print "A matrix ="
    #print_matrix(x, y, A)
    print "Optimal Score =", best
    print "Max location in matrix =", optloc
    '''
    # return the opt score and the best location
    return best, optloc


'''
type1:
       forward reads:                      XXXX...XXXXNAGATCGGAAGAGCGGTTCAG
                                           ||||...||||
       reverse reads: TGTGCTGCGAGAAGGCTAGANXXXX...XXXX

type2:
       forward reads: XXXXX...XXXXXXXXXXX...XXXX
                                     ||||...||||
       reverse reads:                XXXX...XXXXXXXXXXX...XXXX
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
    name="2_"+os.path.basename(args.input1).split('.')[0][:-2]
    output1=open("short_"+name+".fasta",'w')
    output2=open("long_"+name+".fasta",'w')
    output3=open("evenlong_"+name+"_1.fastq",'w')
    output4=open("evenlong_"+name+"_2.fastq",'w')
    output5=open("wierd_"+name+"_1.fastq",'w')
    output6=open("wierd_"+name+"_2.fastq",'w')

    print >>sys.stderr,"start search..."
    type1=type2=type3=type4=total=0
    for rec1, rec2 in itertools.izip(fastq_iter1,fastq_iter2):
        seq1=str(rec1.seq)
        seq2=str(rec2.seq)
        seq2_RC=str(rec2.seq.reverse_complement())
        score1,loc1=local_align(primer1,seq1)
        score2,loc2=local_align(primer2,seq2)
        if (score1>=7*loc1[0]) or (score2>=7*loc2[0]):
            if (score1>=7*loc1[0]) and ((score2>=7*loc2[0])):
                fragment=rec1[:(loc1[1]-loc1[0])]
                SeqIO.write(fragment,output1,"fasta")
                type1+=1
            else:
                SeqIO.write(rec1,output5,'fastq')
                SeqIO.write(rec2,output6,'fastq')
                type4+=1
            #print >> sys.stderr, "type1\n"
        else:
            score,loc=local_align(seq1,seq2_RC)
            if (score>=7*min(loc)) and (loc[1]<=loc[0]) and (loc[0]>=len(seq1)-2):
                fragment=SeqRecord(rec1.seq+rec2.seq.reverse_complement()[loc[1]+1:],id=rec1.id,name=rec1.name,description=rec1.description)
                SeqIO.write(fragment,output2,"fasta")
                #print >> sys.stderr, "type2\n"
                type2+=1

            elif (score>=8*min(loc)) and (loc[1]>loc[0]) and (loc[0]>=len(seq1)-10):
                SeqIO.write(rec1[:loc[0]],output1,"fasta")
                type1+=1
                #print >> sys.stderr, "type1\n"
            else:
                SeqIO.write(rec1,output3,'fastq')
                SeqIO.write(rec2,output4,'fastq')
                #print >> sys.stderr, "type3\n"
                type3+=1
        total+=1
        if total%10000==0:
            print >> sys.stderr, " %d read pairs processed and there are %d type1(short), %d type2(long), %d type3(evenlong), and %d type4(wierd). \r "%(total, type1, type2, type3,type4)
    output1.close()
    output2.close()
    output3.close()
    output4.close()
    output5.close()
    output6.close()

if __name__ == '__main__':
    main()
