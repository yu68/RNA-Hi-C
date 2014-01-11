import sys,argparse,os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
#plt.ioff()

def ParseArg():
    p=argparse.ArgumentParser(description="plot linked pairs around a given interaction. information of linked pairs are stored in file '*_fragment_paired_align.txt'",epilog="Require: matplotlib, numpy")
    p.add_argument("interaction",type=str,help="Interaction file from output of 'Select_strongInteraction_pp.py'")
    p.add_argument("linkedPair",type=str,help="file for information of linked pairs, which is output of 'Stitch-seq_Aligner.py'")
    p.add_argument('-n',type=int,default=1,help="Plot all linked pairs for 'n'th interaction in the interaction file, default=1")
    p.add_argument('-s','--start',type=int,nargs='+',default=(7,8),help='start column number of the second region in interaction file and linkedPair file, default=(7,8)')
    p.add_argument('-d','--distance',type=int,default=10,help='the plus-minus distance (unit: kbp) flanking the interaction regions to be plotted, default=10')
    p.add_argument('-o','--output',type=str,help="output plot file, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()

class Bed:
    def __init__(self,x):
        self.chro=x[0]
        self.start=int(x[1])
        self.end=int(x[2])
        if len(x)>=6:
            self.type=x[3]
            self.name=x[4]
            self.subtype=x[5]
        self.center=int((self.start+self.end)/2)
    def overlap(self,bed2,n):
        return (self.chro==bed2.chro) and (self.start+n<bed2.end) and (bed2.start+n<self.end)
    def str_region(self):
        return self.chro+":%d-%d"%(self.start,self.end)

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


def transform(x1,start1,end1,start2,end2):
    '''
    transfrom from scale 2 to scale 1
    '''
    return int(1.0*(x1-start1)/(end1-start1)*(end2-start2)+start2)

def Main():
    args=ParseArg()
    distance=args.distance*1000

    print "\nChecking if linkedPair file is tabixed..."
    if not os.path.isfile(args.linkedPair+".tbi"):
        print "  tabix-ing..."
        os.system("sort -k1,1 -k2,2n "+args.linkedPair+" > temp.txt")
        os.system("bgzip temp.txt")
        os.system("tabix -p bed temp.txt.gz")
        linkedPair='temp.txt.gz'
    else:
        linkedPair=args.linkedPair
    print "  linkedPair file is tabixed."

    # start column number for second regions
    # s1 for interaction file and s2 for linkedPair file
    (s1,s2)=args.start

    print "\nExtracting interaction information..."
    Interactions=open(args.interaction,'r')
    l=Interactions.read().split('\n')[args.n-1].split('\t')
    part1=Bed(l[0:6])
    part2=Bed(l[s1:(s1+6)])
    start1=part1.start-distance
    end1=part1.end+distance
    start2=part2.start-distance
    end2=part2.end+distance

    
    # initialize figure
    print "\n Start plot interaction: "+part1.str_region()+" <-> "+part2.str_region()
    col1="#4F81BD"
    col2="#C0504D"
    fig = plt.figure(figsize=(8,4))
    ax1 = plt.subplot(111,frameon=False,yticks=[])
    plt.tick_params(axis="y",which="both",left="off",right="off",labelleft="off")  # remove y ticks
    plt.subplots_adjust(top=0.75)
    ax2 = ax1.twiny()
    ax1.set_xlim(start1,end1)
    ax2.set_xlim(start2,end2)
    ax1.set_ylim(0,1)
    ax2.set_ylim(0,1)

    #set x ticks withour offset
    locs=ax1.get_xticks()
    ax1.set_xticklabels(map(lambda x: "%i"%x, locs),fontsize=8)
    locs=ax2.get_xticks()
    ax2.set_xticklabels(map(lambda x: "%i"%x, locs),fontsize=8)
    
    y_1=0.36
    y_2=0.64
    ax1.add_patch(matplotlib.patches.Rectangle((part1.start,y_1),part1.end-part1.start,0.04,color=col1))
    ax2.add_patch(matplotlib.patches.Rectangle((part2.start,y_2),part2.end-part2.start,0.04,color=col2))
    ax1.plot([start1,end1],[y_1+0.02,y_1+0.02],color=col1,linewidth=1,alpha=0.7)
    ax2.plot([start2,end2],[y_2+0.02,y_2+0.02],color=col2,linewidth=1,alpha=0.7)


    print "\nQuery linkedPairs within +-%dkbp of interaction"%(distance/1000)
    os.system("tabix "+linkedPair+" %s:%i-%i > temp2.txt"%(part1.chro,part1.start-distance,part1.end+distance))
    print "\nList of linked pairs plotted: "
    for b in read_interaction("temp2.txt",s2):
        if part1.overlap(b[0],-distance) and part2.overlap(b[1],-distance):
            x1_2_start=transform(b[0].start,start1,end1,start2,end2)
            x1_2_end=transform(b[0].end,start1,end1,start2,end2)
            ax2.plot([(x1_2_start+x1_2_end)/2,b[1].center],[y_1+0.02,y_2+0.02],"ko-",
                     markersize=3,markeredgewidth=0,alpha=0.3)
            print "  "+b[0].str_region()+" <-> "+b[1].str_region()
    plt.text(0.5, 1.12, part1.str_region()+" <-> "+part2.str_region(),
         horizontalalignment='center',
         fontsize=10,
         transform = ax1.transAxes)
    ax1.text(part1.center,y_1-0.03,"|".join([part1.type,part1.name,part1.subtype]),
             verticalalignment='center', horizontalalignment='center',fontsize=10,color=col1)
    ax2.text(part2.center,y_2+0.07,"|".join([part2.type,part2.name,part2.subtype]),
             verticalalignment='center', horizontalalignment='center',fontsize=10,color=col2)
    ax1.set_ylim(0.1,0.9)
    ax2.set_ylim(0.1,0.9)
    ax1.text(start1, 0.15, part1.chro,horizontalalignment='left',fontsize=10)
    ax2.text(start2, 0.85, part2.chro,horizontalalignment='left',fontsize=10)
    plt.savefig(args.output)
    
    # remove temp file
    if not os.path.isfile(args.linkedPair+".tbi"):
        os.system("rm temp.txt.gz*")
    os.system("rm temp2.txt")

if __name__=="__main__":
    Main()
