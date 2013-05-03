#!/usr/bin/env python
import sys,os,argparse
import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

'''
Change zip() to itertools.izip for python 2.?
Change range() to xrange() for python 2.?
'''

'''
http://www.cs.umd.edu/class/fall2010/cmsc423/align423.py
'''

def ParseArg():
    p=argparse.ArgumentParser( description = 'Select miRNA_mRNA pairs and remove linker sequences', epilog = 'Library dependency: Bio, itertools')
    p.add_argument('input1',type=str,metavar='reads1',help='forward input fastq/fasta file')
    p.add_argument('input2',type=str,metavar='reads2',help='reverse input fastq/fasta file')
    p.add_argument('linker',type=str,metavar='linker',help='file contianing linker sequences')
    
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

def local_align(x, y, score=ScoreParam(10, -5, -7)):
    """Do a local alignment between x and y with the given scoring parameters.
    We assume we are MAXIMIZING."""

    # create a zero-filled matrix
    A = make_matrix(len(x) + 1, len(y) + 1)

    best = 0
    optloc = (0,0)

    # fill in A in the right order
    for i in xrange(1, len(x)+1):
        for j in xrange(1, len(y)+1):

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
    print "A matrix ="
    print_matrix(x, y, A)
    print "Optimal Score =", best
    print "Max location in matrix =", optloc
    '''
    # return the opt score and the best location
    return best, optloc



def locate_linker(record,linkers,reverse=False):
    end_loc=len(record)
    real_score=0
    for linker in linkers:
        if reverse == True:
            linker=str(Seq(linker,IUPAC.unambiguous_dna).reverse_complement())
        score,loc=local_align(linker,record)
        if (end_loc > (loc[1]-loc[0]))&(score > 8*loc[0]):
            end_loc = (loc[1]-loc[0])
            real_score = score
    return real_score, end_loc


'''
Seven situation:

1. NNNXXXXNN - miRNA - AUC UGG UAA UCC GUA UAA AGU AUG UUG AUG UUC CAA UAA GCA GAU CAU GUU UUU UAA GCC GUC A - mRNA
2. NNNXXXXNN - miRNA - UAA GCA GAU CAU GUU UUU UAA GCC GUC A - mRNA
3. NNNXXXXNN - miRNA - AUC UGG UAA UCC GUA UAA AGU AUG UUG AUG UUC CAA - mRNA
4. NNNXXXXNN - miRNA - AUC UGG UAA UCC GUA UAA AGU AUG UUG AUG UUC CAA
5. NNNXXXXNN - UAA GCA GAU CAU GUU UUU UAA GCC GUC A - mRNA 
6. NNNXXXXNN - mRNA
7. NNNXXXXNN - miRNA (less likely)

'''

def Main():
    args=ParseArg()
    #---------------- read pair-end seq ----------------
    fastq_iter1 = SeqIO.parse(open(args.input1),"fastq")
    fastq_iter2 = SeqIO.parse(open(args.input2),"fastq")
    #---------------------------------------------------

    #---------------- read linker seq ------------------
    linkers=[]
    for i in open(args.linker,'r'):
        i=i.strip()
        linkers.append(i)
    #---------------------------------------------------
    name1=args.input1.split("/")[-1]
    name2=args.input2.split("/")[-1]
    outfile1 = open("Paired"+name1,"w")
    outfile2 = open("Paired"+name2,"w")
    outfile3 = open("Single_miRNA.fastq",'w')
    outfile4 = open("Single_mRNA.fastq",'w')
    

    for rec1, rec2 in itertools.izip(fastq_iter1,fastq_iter2):
        real_score1, end_loc1 = locate_linker(rec1.seq,linkers)
        real_score2, end_loc2 = locate_linker(rec2.seq,linkers, reverse=True)

        # situation 1,2,3:
        if (15<end_loc1<50) & (end_loc2>20) :
            SeqIO.write(rec1[9:end_loc1],outfile1,"fastq")
            SeqIO.write(rec2[:end_loc2],outfile2,"fastq")
        # situation 4:
        elif (15<end_loc1<50) & (0<=end_loc2<3) :
            SeqIO.write(rec1[9:end_loc1],outfile3,"fastq")
        # situation 5:
        elif (9<=end_loc1<=15) & (end_loc2>20) :
            SeqIO.write(rec2[:end_loc2],outfile4,"fastq")
    outfile1.close()
    outfile2.close()
    outfile3.close()
    outfile4.close()

if __name__ == '__main__':
    Main()
    
            
            
            
            


