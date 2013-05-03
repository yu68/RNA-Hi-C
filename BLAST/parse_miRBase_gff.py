#-------------------------------------------------------------------------------
# Name:        parse_miRBase_gff
# Purpose:     parse gff3 file from miRBase and separate into pre-/mature miRNA bed file
#
# Author:      Pengfei
#
# Created:     09/12/2012
# Copyright:   (c) Pengfei 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os,sys,argparse

class gff(object):
    '''
    gff3 file from miRBase website:
        chr1	.	miRNA_primary_transcript	12425986	12426106	.	+	.	ID=MI0021869_1;accession_number=MI0021869;Name=mmu-mir-6341
    '''
    def __init__(self,x=None,**kwargs):
        self.name="noname"
        self.type=""
        self.id=""
        self.strand="."
        if x is not None:
            x=x.split("\t")
            self.chr="chr"+str(x[0]).strip()
            self.type=x[2].strip()
            self.start=int(x[3])
            self.end=int(x[4])
            self.strand=x[6].strip()
            info=x[8].split(";")
            try:
                self.id=info[0].split("=")[1]
                self.name=info[2].split("=")[1]
            except:
                pass
    def __str__(self):
        string=self.chr+"\t"+str(self.start)+"\t"+str(self.end)+"\t"+self.name+"\t"+self.type
        string+="\t"+str(self.strand)
        return string
    def __len__(self):
        return self.end-self.start
    def __cmp__(self,other):
	    return cmp(self.chr,other.chr) or cmp(self.start,other.start) or cmp(self.end,other.end)

def GffIterator(handle,**kwargs):
    '''
    '''
    if type(handle)==type("s"):
        try:
            handle=handle.strip()
            x=handle.split(".")
            if x[-1]=="gz":
                handle=gzip.open(handle,"r")
            else:
                handle=open(handle,"r")

        except:
            raise ValueError("Can't open file %s"%handle)
    for line in handle:
        line=line.strip()
        if len(line)==0:
            continue
        if line[0]=="#":
            continue
        yield gff(line)

def ParseArg():
    p=argparse.ArgumentParser( description = "parse gff3 file from miRBase and separate into pre-/mature miRNA bed file")
    p.add_argument("-i","--input",dest="input",type=str,help="input gff3 file")
    p.add_argument("-o","--output",dest="output",type=str,help="output_prefix,two file generated: PREFIX_mature.bed,PREFIX_pre.bed")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()


def main():
    global args
    args=ParseArg()
    if args.output=="stdout":
        out_p=sys.stdout
        out_m=sys.stdout
    else:
        try:
            out_m=open(args.output+"_mature.bed","w")
            out_p=open(args.output+"_pre.bed","w")
        except:
            pass
    for record in GffIterator(open(args.input)):
        if record.type=="miRNA_primary_transcript":
            print >>out_p,"\t".join(str(f) for f in [record.chr,record.start,record.end,record.name,record.id,record.strand])
        if record.type=="miRNA":
            print >>out_m,"\t".join(str(f) for f in [record.chr,record.start,record.end,record.name,record.id,record.strand])
    out_m.close()
    out_p.close()


if __name__ == '__main__':
    main()
