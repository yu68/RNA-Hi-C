#!/usr/bin/env python
# Pengfei Yu yupf05@gmail.com

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import itertools
import os,sys, argparse
import time

def ParseArg():
    p=argparse.ArgumentParser( description = 'Remove duplicated reads which have same sequences for both forward and reverse reads. Choose the one appears first.', epilog = 'Library dependency: Bio, itertools')
    p.add_argument('input1',type=str,metavar='reads1',help='forward input fastq/fasta file')
    p.add_argument('input2',type=str,metavar='reads2',help='reverse input fastq/fasta file')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()


def Main():
    '''
    http://bytesizebio.net/index.php/2011/08/25/short-bioinformatics-hacks-merging-fastq-files/
    http://stackoverflow.com/questions/1215208/how-might-i-remove-duplicate-lines-from-a-file
    '''
    import time
    t0=time.time()
    count=0
    unique=0
    Unique_seqs=set()
    args=ParseArg()
    name1=args.input1.split('/')[-1]
    # or alternative: name1=os.path.basename(args.input1)
    name2=args.input2.split('/')[-1]
    outfile1 = open("Rm_dupPE_"+name1,"w")
    outfile2 = open("Rm_dupPE_"+name2,"w")
    fastq_iter1 = SeqIO.parse(open(args.input1),"fastq")
    fastq_iter2 = SeqIO.parse(open(args.input2),"fastq")
    for rec1, rec2 in itertools.izip(fastq_iter1, fastq_iter2):
        count=count+1
        if str((rec1+rec2).seq) not in Unique_seqs:
            SeqIO.write(rec1,outfile1,"fastq")
            SeqIO.write(rec2,outfile2,"fastq")
            Unique_seqs.add(str((rec1+rec2).seq))
            unique=unique+1
    outfile1.close()
    outfile2.close()
    time=time.time()-t0
    print "from %i pair-end read pairs, there are %i unique pairs, %i PCR duplictes removed." %(count, unique, count-unique)
    print "Time cost: %.2f s"%(time)

if __name__ == '__main__':
    Main()
