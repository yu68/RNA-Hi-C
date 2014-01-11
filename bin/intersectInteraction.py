import sys,argparse,os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
plt.ioff()

def ParseArg():
    p=argparse.ArgumentParser(description="find intersections (overlaps) between two interaction sets",epilog="require numpy and matplotlib if set '-c'")
    p.add_argument("-a","--filea",type=str,required=True,help="file for interaction set a")
    p.add_argument("-b","--fileb",type=str,required=True,help="file for interaction set b")
    p.add_argument("-s","--start",type=int, default=7,help="start column number of the second part in each interaction (0-based),  default:7")
    p.add_argument('-n',"--nbase",type=int,default=1,help="number of overlapped nucleotides for each part of interactions to call intersections, default: 1")
    p.add_argument('-o','--output',type=str,help="specify output file")
    p.add_argument("-c",'--compare',action='store_true',help="Use a set of different 'nbase' to call overlaps and find the best one. if nbase=-200, then choose from [0,-10,-20,...,-200]")
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
if args.compare:
    count={}
for a in read_interaction(args.filea,args.start):
    k=0
    num+=1
    print >>sys.stderr,"Prcessing interactions in set a: %d\r"%(num),
    os.system("tabix temp.txt.gz %s:%i-%i > temp2.txt"%(a[0].chro,a[0].start+N,a[1].end-N))
    for b in read_interaction("temp2.txt",args.start):        
        if a[0].overlap(b[0],N) and a[1].overlap(b[1],N):
            k=k+1
            print >>out, "# %d-%d"%(m,k)
            print >>out,a[2]
            print >>out,b[2]
    if k>0:
        m+=1 
    if not args.compare: continue
    for i in range(abs(N)/10+1):
        k=0
        for b in read_interaction("temp2.txt",args.start):
            if a[0].overlap(b[0],i*10*N/abs(N)) and a[1].overlap(b[1],i*10*N/abs(N)): 
                k=k+1
        if k>0:
            if i*10*N/abs(N) not in count:
                count[i*10*N/abs(N)]=1
            else:
                count[i*10*N/abs(N)]+=1

if args.compare:
    print
    x=[]
    y=[] 
    for i in range(abs(N)/10+1):
        x.append(i*10)
        y.append(count[i*10*N/abs(N)])
    x=np.asarray(x)
    y=np.asarray(y)
    fig, ax = plt.subplots()
    rects = ax.plot(x,y,'r.-')
    ax.set_ylabel("# of overlaps")
    namea=os.path.basename(args.filea).split("_")[0]
    nameb=os.path.basename(args.fileb).split("_")[0]
    ax.set_title(namea+" vs "+nameb)
    ax.set_xlabel("allowed distance")
    ax.set_ylim(0,1.1*max(y))
    plt.savefig('test.pdf')

os.system("rm temp.txt.gz*")
os.system("rm temp2.txt")
print >> sys.stderr, "Number of interactions in set a overlapped with those in set b is: %d"%(m-1)
