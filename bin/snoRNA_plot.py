import sys,os,argparse
import numpy as np
import matplotlib.pyplot as plt
import string
from matplotlib.colors import LogNorm
from pylab import cm

def ParseArg():
    p=argparse.ArgumentParser(description='plot read stack in snoRNA regions')
    p.add_argument("ID",type=str,help="Ensembl ID of snoRNA")
    p.add_argument("linkedPair",type=str,help="file for information of linked pairs, which is output of 'Stitch-seq_Aligner.py'")
    p.add_argument('-o','--output',type=str,help="output plot file, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()

def revcomp(seq, rev_table):
    return seq.translate(rev_table)[::-1]


def Main():
    GeneLoc='/home/yu68/bharat-interaction/new_lincRNA_data/all_RNAs-rRNA_repeat.txt'
    args=ParseArg()
    ID = args.ID
    gene=os.popen('grep "'+ID+'" '+GeneLoc).read().split('\n')[0].split('\t')
    geneStart=int(gene[1])
    geneEnd=int(gene[2])
    geneChro=gene[0]
    geneName=gene[4]
    geneStrand=gene[5]
    if not os.path.isfile(args.linkedPair):
        print "LinkedPair file is not exist, please check!!"
        sys.exit(0)
    if geneStrand=='+':
        Pairs=os.popen('grep "'+ID+'" '+args.linkedPair+" | awk '$7==$16' | sort -k3,3n -k11,11n").read().split('\n')
        print len(Pairs)
    else:
        Pairs=os.popen('grep "'+ID+'" '+args.linkedPair+" | awk '$7==$16' | sort -k2,2n -k12,12n").read().split('\n')

    fig = plt.figure(figsize=(4,4))
    ax = plt.subplot(121,frameon=False,yticks=[],xticks=[])
    ax2 = plt.subplot(122)
    ax.set_xlim(1.2*geneStart-0.2*geneEnd,1.2*geneEnd-0.2*geneStart)
    
    #set x ticks withour offset
    locs=ax.get_xticks()
    ax.set_xticklabels(map(lambda x: "%i"%x, locs),fontsize=8)
    
    ax.bar(geneStart,0.035,geneEnd-geneStart,facecolor="#5E5E5E",edgecolor='k',bottom=0.01,lw=1)
    if geneStrand=="+":
        geneName+="-->"
    else:
        geneName="<--"+geneName
    ax.text((geneStart+geneEnd)/2,0.02,geneName,horizontalalignment='center',fontsize=8)

    n_pairs = len(Pairs)
    
    space= 0.95/len(Pairs)  # space for one pair
    i=0
    col1="#4F81BD"
    col2="#C0504D"

    geneLen = geneEnd-geneStart+1
    matrix = np.zeros((geneLen,geneLen))
    matrix = matrix+0.01        
    for p in Pairs:
        if p.strip()=='': continue
        p=p.split('\t')
        assert(geneChro==p[0]==p[9])
        p1_start=int(p[1])
        p1_end=int(p[2])
        p2_start=int(p[10])
        p2_end=int(p[11])
        y_bottom=i*space+0.05
        ax.bar(p1_start,0.7*space,p1_end-p1_start,facecolor=col1,edgecolor=col1,bottom=y_bottom,lw=0.1)
        ax.bar(p2_start,0.7*space,p2_end-p2_start,facecolor=col2,edgecolor=col2,bottom=y_bottom,lw=0.1)
        rev_table=string.maketrans('ACGTacgt', 'TGCAtgca')
        if p[3]=='+':
            print p1_end-geneStart,p[4][-10:],
            x=p1_end-5-geneStart
            ax.bar(p1_end-5,0.7*space,5,facecolor='k',edgecolor='k',bottom=y_bottom,lw=0.1)
        else:
            print p1_start-geneStart,revcomp(p[4][:10],rev_table),
            x=p1_start+5-geneStart
            ax.bar(p1_start,0.7*space,5,facecolor='k',edgecolor='k',bottom=y_bottom,lw=0.1)
        if p[12]=='+':
            print p2_end-geneStart,revcomp(p[13][-10:],rev_table)
            y=p2_end-5-geneStart
            ax.bar(p2_end-5,0.7*space,5,facecolor='k',edgecolor='k',bottom=y_bottom,lw=0.1)
        else:
            print p2_start-geneStart,p[13][:10]
            y=p2_start+5-geneStart
            ax.bar(p2_start,0.7*space,5,facecolor='k',edgecolor='k',bottom=y_bottom,lw=0.1)
        x = min(max(x,0),geneLen-1)
        y = min(max(y,0),geneLen-1)
        matrix[x,y] +=1
        i=i+1
    cax2=ax2.matshow(matrix,norm=LogNorm(vmin=0.01, vmax=matrix.max()))
    fig.colorbar(cax2)
    plt.savefig(args.output)
    plt.show()

if __name__=="__main__":
    Main()
