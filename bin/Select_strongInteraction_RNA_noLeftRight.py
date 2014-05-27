import sys,argparse
from time import time
import copy,math
from data_structure import *
from scipy.stats import hypergeom

def ParseArg():
    p=argparse.ArgumentParser(description="find strong interactions between RNAs from annotation of RNAs, intron regions are removed, no self interactions",epilog="need Scipy for hypergeometric distribution")
    p.add_argument("-i","--input",type=str,required=True,help="input file which is the output file of Stitch-seq-Aligner.py")
    p.add_argument("-M","--min_clusterS",type=int,default=5,help="minimum number of segments allowed in each cluster, default:5")
    p.add_argument("-m","--min_interaction",type=int,default=3,help="minimum number of interactions to support a strong interaction, default:3")
    p.add_argument('-p',"--p_value",type=float,default=0.05,help="the p-value based on hypergeometric distribution to call strong interactions, default: 0.05")
    p.add_argument('-o','--output',type=str,help="specify output file")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()

def SingleFragment(p1,p2):
    '''
    determine if mapped pair is single RNA fragment or not.
      1. no linker
      2. same chromosome
      3. same RNA
    '''
    if (p1.end-p1.start==90) and (p2.end-p2.start==100) and (p1.chr==p2.chr) and p1.name==p2.name:
        return True


def Main():
    t1=time()
    args=ParseArg()
    inp = open(args.input, 'r')
    min_clusterS=args.min_clusterS
    min_interaction=args.min_interaction
    p_value=args.p_value
    output=open(args.output,'w')

    #store count of RNA for part1 and part2
    part={}


    k=0
    
    print >> sys.stderr,"# Inputing data..."
    interaction = {}  # store number of interactions for different RNA

    Types = ["snoRNA","protein_coding","snRNA","lincRNA","tRNA","misc_RNA","pseudogene","miRNA","antisense","sense_intronic","non_coding","processed_transcript"]
    for line in inp.read().split('\n'):
        if line=='': continue
        line=line.strip().split('\t')
        p1=annotated_bed(line[0:8],id=k,cluster=1)
        p2=annotated_bed(line[9:],id=k,cluster=1)
        if SingleFragment(p1,p2): continue
        k+=1
        if p1.subtype=="intron" or p2.subtype=="intron": continue
        if p1.type in Types:
            p1_name = p1.chr+":"+p1.name
            if p1_name not in part:
                part[p1_name]=1
            else:
                part[p1_name]+=1  
        if p2.type in Types:
            p2_name = p2.chr+":"+p2.name
            if p1.chr == p2.chr and p1.name == p2.name: continue # count once for self-interaction
            if p2_name not in part:
                part[p2_name]=1
            else:
                part[p2_name]+=1
        if p1.type in Types and p2.type in Types:
            if p1_name == p2_name: continue 
            if p1_name>p2_name:
                temp = p1
                p1 = p2
                p2 = temp
            inter_name = p1.chr+":"+p1.name+"--"+p2.chr+":"+p2.name
            if inter_name not in interaction:
                interaction[inter_name]=[copy.deepcopy(p1),copy.deepcopy(p2)]
            else:
                interaction[inter_name][0].Update(p1.start,p1.end)
                interaction[inter_name][1].Update(p2.start,p2.end)
                interaction[inter_name][0].cluster+=1
        if k%20000==0: 
            print >> sys.stderr,"  Reading %d pairs of segments\r"%(k),
    print >> sys.stderr,"Get total %d pairs."%(k)

    print >>sys.stderr,"   number of different RNAs is %d          "%(len(part))
    
    total = k # total pairs used
    n=0
    k=0  # record number of strong interactions
    for i in interaction:
        n+=1
        count = interaction[i][0].cluster
        if count < min_interaction: continue
        p1_name = i.split("--")[0]
        p2_name = i.split("--")[1]
        P1 = interaction[i][0]
        P2 = interaction[i][1]
        P1.cluster = part[p1_name]
        P2.cluster = part[p2_name]
        if part[p1_name]<min_clusterS or part[p2_name]<min_clusterS: continue
        real_p=1-hypergeom.cdf(count,total,part[p1_name],part[p2_name])
        if real_p<=p_value:
            k=k+1
            try:
                log_p = math.log(real_p)
            except:
                log_p = -float("Inf")
            print >> output, str(P1)+'\t'+str(P2)+'\t%d\t%.4f'%(count,log_p)
        if n%100==0: print >> sys.stderr, "  Progress ( %d / %d )\r"%(n,len(interaction)),
    print >> sys.stderr,"# Find %d strong interactions. Cost time: %.2f s"%(k,time()-t1)

if __name__=="__main__":
    Main()
