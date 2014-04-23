import sys,argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import scipy


def ParseArg():
    p=argparse.ArgumentParser( description='Draw histogram for sequence length distribution',epilog='Library dependency:Bio, matplotlib, scipy')
    group=p.add_mutually_exclusive_group()
    group.add_argument("-f","--fasta",action='store_true',help='add this option for fasta input file')
    group.add_argument("-q","--fastq",action='store_true',help='add this option for fastq input file')
    p.add_argument("-i","--input",dest='input',type=str,help='input fasta/fastq file')
    p.add_argument("-o","--output",dest='output',type=str,help='output file name, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

args=ParseArg()
file=args.input

output=args.output

if args.fastq:
    type="fastq"
elif args.fasta:
    type="fasta"

seqLen=[]
reads=SeqIO.parse(open(file),type)
for read in reads:
    seqLen.append(len(read.seq))


if max(seqLen)<=150:
    x1=0
    x2=150
elif max(seqLen)<=300:
    x1=150
    x2=300
print >>sys.stderr,"start drawing..."
plt.hist(seqLen,bins=150,normed=1,range=(x1,x2))
density=gaussian_kde(scipy.array(seqLen,float))
xs = np.linspace(x1,x2,x2-x1)
density.covariance_factor = lambda : 0.05
density._compute_covariance()
print "KDE density computated"
plt.plot(xs,density(xs),label='KDE',color='r')
plt.xlabel("Length (nt)")
plt.ylabel("Density")
plt.title("Histogram")
plt.legend()
plt.savefig(output)
print "the average flagment length is %f, the SD is %f" %(scipy.mean(seqLen),scipy.std(seqLen)) 

reads.close()

