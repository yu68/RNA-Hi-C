import sys,argparse,os
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import scipy.stats
plt.ioff()

def ParseArg():
    p=argparse.ArgumentParser(description="find intersections (overlaps) between two interaction sets",epilog="require 'random'&'numpy'&'scipy' module if set '-p'")
    p.add_argument("-a","--filea",type=str,required=True,help="file for interaction set a")
    p.add_argument("-b","--fileb",type=str,required=True,help="file for interaction set b")
    p.add_argument("-s","--start",type=int, default=7,help="start column number of the second part in each interaction (0-based),  default:7")
    p.add_argument('-n',"--nbase",type=int,default=1,help="number of overlapped nucleotides for each part of interactions to call intersections, default: 1")
    p.add_argument('-t',"--type",type=int,default=1,help="Output type: 1. print each overlap in three lines (overlap ID, interaction_a, interaction_b); 2. only output all interaction_a that overlapping with interaction_b")
    p.add_argument('-o','--output',type=str,help="specify output file")
    p.add_argument("-p",'--pvalue',action='store_true',help="calculate p-values based on 100times permutations")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()

class Bed:
    def __init__(self,x):
        self.chro=x[0]
        self.start=int(x[1])
        self.end=int(x[2])
    def overlap(self,bed2,n):
        return (self.chro==bed2.chro) and (self.start+n<bed2.end) and (bed2.start+n<self.end)

def read_interaction(File,s):
    '''
    s: start column number for second part of interaction
    '''
    a=open(File,'r')
    for l in a.read().split('\n'):
        if l.strip()=="": continue
        lsep=l.split('\t')
        bed1=Bed(lsep[0:3])
        bed2=Bed(lsep[s:(s+3)])
        yield (bed1,bed2,l)

args=ParseArg()
N=args.nbase # overlaped nucleotide number
out=open(args.output,'w')

# tabix file b
os.system("sort -k1,1 -k2,2n "+args.fileb+" > temp.txt")
os.system("bgzip temp.txt")
os.system("tabix -p bed temp.txt.gz")


m=1 # store the number of intersections in set a
num=0 # store number of interactions in set a that have been proessed
for a in read_interaction(args.filea,args.start):
    k=0
    num+=1
    print >>sys.stderr,"Prcessing interactions in set a: %d\r"%(num),
    os.system("tabix temp.txt.gz %s:%i-%i > temp2.txt"%(a[0].chro,a[0].start+N,a[0].end-N))
    for b in read_interaction("temp2.txt",args.start):        
        if a[0].overlap(b[0],N) and a[1].overlap(b[1],N):
            k=k+1
            if args.type==1:
                print >>out, "# %d-%d"%(m,k)
                print >>out,a[2]
                print >>out,b[2]
    if k>0:
        if args.type==2:
            print >>out,a[2]
        m+=1 

real_inter_n = m-1

print >> sys.stderr, "Number of interactions in set a overlapped with those in set b is: %d"%(real_inter_n)

if not args.pvalue:
    os.system("rm temp.txt.gz*")
    os.system("rm temp2.txt")
    sys.exit(0)

# for permutation here
part1_a=[]
part2_a=[]
for a in read_interaction(args.filea,args.start):
    part1_a.append(a[0])
    part2_a.append(a[1])

print >> sys.stderr, "\nStart permutation..."
perm_N = 100 
perm_positive = 0   # store permutation number larger than real
perm_num = []
for i in range(perm_N):
    print >>sys.stderr,"  Permutation -  %d\r"%(i),
    random.shuffle(part2_a)
    m=0
    for j in range(num):
        a0 = part1_a[j]
        a1 = part2_a[j]
        k = 0
        os.system("tabix temp.txt.gz %s:%i-%i > temp2.txt"%(a0.chro,a0.start+N,a0.end-N))        
        for b in read_interaction("temp2.txt",args.start):
            if a[0].overlap(b[0],N) and a[1].overlap(b[1],N):
                k=k+1
        if k>0:
            m+=1
    if m >= real_inter_n:
        perm_positive+=1
    perm_num.append(m)

mean = np.mean(perm_num)
sd = np.std(perm_num)
p_value = 1-scipy.stats.norm(mean, sd).cdf(real_inter_n)
print >> sys.stderr," Permutation P-value: %.4f, mean: %.2f, sd: %.2f"%(perm_positive*1.0/perm_N,mean,sd)
print >> sys.stderr," Permutation distribution P-value: %.6f  "%(p_value)

os.system("rm temp.txt.gz*")
os.system("rm temp2.txt")
