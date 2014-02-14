import sys,os,argparse
import matplotlib.pyplot as plt
import string

def ParseArg():
    p=argparse.ArgumentParser(description='plot read stack in snoRNA regions')
    p.add_argument("ID",type=str,help="Ensembl ID of snoRNA")
    p.add_argument("partner",type=str,help="name of the snoRNA interaction partner")
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
    partner = args.partner
    gene=os.popen('grep "'+ID+'" '+GeneLoc).read().split('\n')[0].split('\t')
    geneStart=int(gene[1])
    geneEnd=int(gene[2])
    geneChro=gene[0]
    geneName=gene[4]
    geneStrand=gene[5]
    if not os.path.isfile(args.linkedPair):
        print "LinkedPair file is not exist, please check!!"
        sys.exit(0)
    if geneStrand=="+":
        Pairs=os.popen('grep "'+ID+'" '+args.linkedPair+" | awk '$16==\""+partner+"\"' | sort -k3,3n ").read().split('\n')
        Pairs+=os.popen('grep "'+ID+'" '+args.linkedPair+" | awk '$7==\""+partner+"\"' | sort -k11,11n ").read().split('\n')
    else:
        Pairs=os.popen('grep "'+ID+'" '+args.linkedPair+" | awk '$16==\""+partner+"\"' | sort -k2,2n ").read().split('\n')
        Pairs+=os.popen('grep "'+ID+'" '+args.linkedPair+" | awk '$7==\""+partner+"\"' | sort -k12,12n ").read().split('\n')
    fig = plt.figure(figsize=(4,4))
    ax = plt.subplot(111,frameon=False,yticks=[],xticks=[])
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
   
    rev_table=string.maketrans('ACGTacgt', 'TGCAtgca')
 
    for p in Pairs:
        if p.strip()=='': continue
        p=p.split('\t')
        y_bottom=i*space+0.05
        if p[15]==partner:
          p_start=int(p[1])
          p_end=int(p[2])
          col=col1
          if p[3]=='+':
            ax.bar(p_end-5,0.7*space,5,facecolor='k',edgecolor='k',bottom=y_bottom,lw=0.1)
          else:
            ax.bar(p_start,0.7*space,5,facecolor='k',edgecolor='k',bottom=y_bottom,lw=0.1)
        else:
          p_start=int(p[10])
          p_end=int(p[11])
          col=col2
          if p[12]=='+':
            ax.bar(p_end-5,0.7*space,5,facecolor='k',edgecolor='k',bottom=y_bottom,lw=0.1)
          else:
            ax.bar(p_start,0.7*space,5,facecolor='k',edgecolor='k',bottom=y_bottom,lw=0.1)
        ax.bar(p_start,0.7*space,p_end-p_start,facecolor=col,edgecolor=col,bottom=y_bottom,lw=0.1)
        i=i+1
    plt.savefig(args.output)
    plt.show()

if __name__=="__main__":
    Main()
