import os,sys,argparse
from xplib.Annotation import *
from xplib import TableIO
from xplib import DBI

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'This program can achieve cis-annotation for bed entries (intergenic,intron,exon,UTR', epilog='Library dependency : xplib')
    p.add_argument('-i','--input',dest="input",type=str,default="stdin",help="input bed file")
    p.add_argument('-I','--input_format',dest="input_format",action="store",default="bed",help="input file format")
    p.add_argument('-o','--output',dest="output",type=str,default="stdout",help="output file")
    p.add_argument('-a','--annotations',dest="db",type=str,help="feature annotation file, gene table file from UCSC browser")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def judge_exon(bed,gene):
    "judge if bed block in exon or intron when it is in CDS"
    if gene.exon_count==1:
        return ("Exon")
    else:
        flag=0
        for exon in gene.Exons():
            if (bed.start>=exon.start)&(bed.stop<=exon.stop):
                flag=1
        if flag==1:
            return ("Exon")
        else:
            return ("Intron")




def Main():
    global args
    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",args.output,"to write. Using stdout instead"
            out=sys.stdout

    dbi=DBI.init(args.db,"genebed")
    count={}
    count["Intergenic"]=0
    for x in TableIO.parse(args.input,args.input_format):
        flag=0
        gene=""
        for hit in dbi.query(x):
            flag=1
            if hit.align_id==gene:
                continue
            gene=hit.align_id
            #print hit
            #print hit.cds_start,hit.cds_stop
            if (hit.cds_start==hit.cds_stop):
                if hit.align_id[0:3]=="Mir":
                    loc="MiRNA"
                else:
                    loc="Non-coding"
            elif hit.strand=="+":
                if x.stop<=hit.cds_start:
                    loc="5-UTR"
                elif x.start>=hit.cds_stop:
                    loc="3-UTR"
                else:
                    loc=judge_exon(x,hit)
                        
            else:
                if x.stop<=hit.cds_start:
                    loc="3-UTR"
                elif x.start>=hit.cds_stop:
                    loc="5-UTR"
                else:
                    loc=judge_exon(x,hit)
            print >>out,"\t".join (str(f) for f in [x.chr,x.start,x.stop,x.id,x.score,x.strand,hit.align_id,loc])
            if count.has_key(loc):
                count[loc]+=1
            else:
                count[loc]=1

        if flag==0:
            print >>out, "\t".join (str(f) for f in [x.chr,x.start,x.stop,x.id,x.score,x.strand,"None","Intergenic"])
            count["Intergenic"]+=1
    
    out2=open(args.output.split(".")[0]+".cisStat","w")
    for key in sorted(count.keys()):
        print >>out2,key+"\t"+str(count[key])

if __name__=="__main__":
    Main()
            
