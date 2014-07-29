

""" Description:

@version:  0.0.1
@author:   Pengfei Yu
@contact:  
"""

import sys,argparse
from bx.bbi.bigwig_file import BigWigFile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import random


def ParseArg():
    p=argparse.ArgumentParser( description = 'plot average conservation score near linker connection point',epilog="library denpenddency: bx-python, matplotlib,numpy")
    p.add_argument("linkedPair",type=str,help="file for information of linked pairs, which is output of 'Stitch-seq_Aligner.py'")
    p.add_argument('-p','--phyloP_wig',type=str,default='/data2/sysbio/UCSD-sequencing/mouse.phyloP30way.bw',help='the bigWig file for phyloP scores,defualt: mouse.phyloP30way.bw')
    p.add_argument('-s','--span',type=int,default=500,help="span length from the linker connection point,default=500")
    p.add_argument('-w','--win_l',type=int,default=3,help='smooth window length for ploted average counts, (default:3, no smooth)')
    p.add_argument('-o','--output',type=str,default='test',help='name of output figure file,can be (.pdf, .eps, .png, .jpg,...)')
    if len(sys.argv)==1:
        print >> sys.stderr, p.print_help()
        sys.exit(0)
    return p.parse_args()

class Bed:
    def __init__(self,x,**kwargs):
        self.chr=x[0]
        self.start=int(x[1])
        self.stop=int(x[2])
        if len(x)>=6:
            self.type=x[3]
            self.name=x[4]
            self.subtype=x[5]
        self.center=int((self.start+self.stop)/2)
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def overlap(self,bed2,n):
        return (self.chr==bed2.chr) and (self.start+n<bed2.stop) and (bed2.start+n<self.stop)
    def str_region(self):
        return self.chr+":%d-%d"%(self.start,self.stop)

def read_linkedPair(File,s,span):
    '''
    this is not for general purpose
    s: start column number for second part of interaction
    span: span from interaction connection point
    '''
    a=open(File,'r')
    for l in a.read().split('\n'):
        if l.strip()=="": continue
        lsep=l.split('\t')
        if lsep[0]=='chrM' or lsep[s]=='chrM': continue
        if lsep[3]=="+":
            connectP=int(lsep[2])
        else:
            connectP=int(lsep[1])
        bed1=Bed([lsep[0],connectP-span,connectP+span],strand=lsep[3])
        if lsep[s+3]=="+":
            connectP=int(lsep[s+2])
        else:
            connectP=int(lsep[s+1])
        bed2=Bed([lsep[s],connectP-span,connectP+span],strand=lsep[s+3])
        yield (bed1,bed2,l)

def RandomBed(bed):
    '''
    generate a random region with the same chr and length as the given bed
    '''
    start = random.randint(0,50000000)
    length = bed.stop-bed.start
    Strand = random.choice("+-")
    return(Bed([bed.chr,start,start+length],strand=Strand))

def WigToCount(bed,bw):
    '''
    find intensity of wig within bed range, strand specific
    reso: resolution
    '''
    try:
      array=bw.summarize(bed.chr,bed.start,bed.stop,(bed.stop-bed.start)).sum_data
    except:
      print bed.str_region()
    if bed.strand=='-':
        array=array[::-1]
    return array

def smooth(x,window_len=11,window='hanning'):
        '''
        smooth the count for each interval.
        from:
        http://stackoverflow.com/questions/5515720/python-smooth-time-series-data
        '''
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=mp.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

def Main():
    args=ParseArg()
    linkedPair=args.linkedPair
    phyloP_wig=args.phyloP_wig
    w=args.win_l
    s=9
    span=args.span
    arrays1=[]
    arrays2=[]
    arrays3=[]
    arrays4=[]
    bw=BigWigFile(open(phyloP_wig))
    n=0
    for i in read_linkedPair(linkedPair,s,span):
        bed1=i[0]
        bed2=i[1]
        bed3=RandomBed(bed1)
        bed4=RandomBed(bed2)
        arrays1.append(WigToCount(bed1,bw))
        arrays2.append(WigToCount(bed2,bw))
        arrays3.append(WigToCount(bed3,bw))
        arrays4.append(WigToCount(bed3,bw))
        n+=1
        if n%100==0:
            print >>sys.stderr, "%d th linkedPair now\r"%(n),
        if n>10000:
            break
    arrays1=np.array(arrays1)
    arrays2=np.array(arrays2)
    arrays3=np.array(arrays3)
    arrays4=np.array(arrays4)
    average1=np.sum(arrays1,axis=0)/arrays1.shape[0]
    average2=np.sum(arrays2,axis=0)/arrays2.shape[0]
    average3=np.sum(arrays3,axis=0)/arrays3.shape[0]
    average4=np.sum(arrays4,axis=0)/arrays4.shape[0]
    average2=average2[::-1]
    average4=average4[::-1]
    average1=smooth(average1,w)
    average2=smooth(average2,w)
    average3=smooth(average3,w)
    average4=smooth(average4,w)
    plt.figure(figsize=(8,4))
    plt.xlim=(-span,3*span)
    plt.plot(range(-span,span),average1,label="part1",color="blue")
    plt.plot(range(span,3*span),average2,label="part2",color="red")
    plt.plot(range(-span,span),average3,"b--",label="part1_ctl",color="#8787F9")
    plt.plot(range(span,3*span),average4,"r--",label="part2_ctl",color="#F98C8C")
    plt.bar(-100,0.05,100,facecolor='blue',edgecolor='blue',bottom=-0.05,lw=0.1)
    plt.bar(-10,0.05,10,facecolor='k',edgecolor='k',bottom=-0.05,lw=0.1)
    plt.bar(2*span,0.05,100,facecolor='red',edgecolor='red',bottom=-0.05,lw=0.1)
    plt.bar(2*span,0.05,10,facecolor='k',edgecolor='k',bottom=-0.05,lw=0.1)
    plt.xticks([-span,0,span,2*span,3*span],[-span,0,"...",0,span])
    plt.bar(0,0.01,2*span,facecolor='k',edgecolor='k',bottom=-0.03,lw=0.1)
    plt.ylabel("Conservation Pylop Score")
    plt.xlabel("length (bp)")
    plt.ylim(-0.1,1.0)
    plt.tight_layout()
    plt.savefig(args.output)
    plt.close()
if __name__=="__main__":
    Main()
