import argparse,sys
from xplib import TableIO
from xplib import DBI

def ParseArg():
    p=argparse.ArgumentParser(description = "RNA composition for aligned sample using ensembl annotation",epilog="Library dependency: bam2x")
    p.add_argument("-i","--input",type=str,help="input aligned file")
    p.add_argument("-f","--format",type=str,default="bam",help="input file format, default: bam")
    p.add_argument("-o","--output",type=str,default="stdout",help="output txt file")
    p.add_argument("-a","--annotation",dest="db",type=str,help="annotation bed file")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def Main():
    global args,out
    args=ParseArg()
    if args.output=="stdout":
        out=sys.stdout
    else:
        try:
            out=open(args.output,"w")
        except IOError:
            print >>sys.stderr,"can't open file ",args.output,"to write. Using stdout instead"
            out=sys.stdout

    count={}
    dbi=DBI.init(args.db,"bed")
    for bed in TableIO.parse(args.input,args.format):
        flag=0
        typ="non"
        name="."
        for hit in dbi.query(bed):
            if flag==0:
                name=hit.id.split(".")[1]
                typ=hit.id.split(".")[0]
                flag=1
        if count.has_key(typ):
            count[typ]+=1
        else:
            count[typ]=1
        print >>out, "\t".join (str(f) for f in [bed.chr,bed.start,bed.stop,bed.id,name,bed.strand,typ])

    print >>out, "\n".join ("#"+typ+"\t%d"%(count[typ]) for typ in count.keys())


if __name__=="__main__":
    Main()

