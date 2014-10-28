import sys,argparse,os
from xplib import DBI
import numpy as np
import subprocess


def ParseArg():
    p=argparse.ArgumentParser(description="Statistics for snoRNA-mRNA interactions, percentage of splicing intermediates within all snoRNA-mRNA interactions",epilog="Require: matplotlib, numpy")
    p.add_argument("interaction",type=str,help="Interaction file from output of 'Select_strongInteraction_pp.py', or linked fragment pair file from output of 'Stitch-seq_Aligner.py'")
    p.add_argument("-A","--Annotation",type=str,help="Annotation bed file for locations of mRNA genes")
    p.add_argument('-s','--start',type=int,default=7,help='start column number of the second region in interaction file, default=7')
    p.add_argument('-o','--output',type=str,help="Set to output all snoRNA-mRNA interactions as splicing intermediates, not output if not set")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
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

def read_interaction(File,s):
    '''
    s: start column number for second part of interaction
    '''
    a=open(File,'r')
    for l in a.read().split('\n'):
        if l.strip()=="": continue
        lsep=l.split('\t')
        if lsep[3] in ['+','-']:  # for linked fragment
            bed1=Bed(lsep[0:3]+lsep[5:8],strand=lsep[3])
            bed2=Bed(lsep[s:(s+3)]+lsep[(s+5):(s+8)],strand=lsep[s+3])
        else:  # for strong interactions
            bed1=Bed(lsep[0:6])
            bed2=Bed(lsep[s:(s+6)])
        yield (bed1,bed2,l)


def Main():
    args=ParseArg()
    interaction = args.interaction
    s = args.start
    annotation = args.Annotation
    
    if args.output:
        out = open(args.output,'w')    

    n_snoRNAmRNA=0
    n_splicingIntermediate=0
    for line in read_interaction(interaction,s):
        bed1 = line[0]
        bed2 = line[1]
        if bed1.type=="snoRNA" and bed2.type=="protein_coding":
            snoRNA = bed1
            mRNA = bed2
            n_snoRNAmRNA+=1
        elif bed1.type=="protein_coding" and bed2.type=="snoRNA":
            snoRNA = bed2
            mRNA = bed1
            n_snoRNAmRNA+=1
        else:
            continue
        if bed1.chr!=bed2.chr: continue
        geneBed = subprocess.check_output("grep %s %s"%(mRNA.name,annotation),shell=True)
        gene_region = Bed(geneBed.split("\t")[0:3])
        if gene_region.overlap(snoRNA,1):
            n_splicingIntermediate+=1
            if args.output:
                print >>out, line[2]
    if s==7:
        print >>sys.stdout, "\nsnoRNA-mRNA-interactions:\t%d"%(n_snoRNAmRNA)
    elif s==9:
        print >>sys.stdout, "\nsnoRNA-mRNA-linkedPair:\t%d"%(n_snoRNAmRNA)
    print >>sys.stdout, "SplicingIntermediates:\t%d"%(n_splicingIntermediate)
    print >>sys.stdout, "Percentage\t%.4f%%\n"%(n_splicingIntermediate*100.0/n_snoRNAmRNA)

if __name__=="__main__":
    Main() 
