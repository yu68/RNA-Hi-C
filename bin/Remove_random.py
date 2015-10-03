import sys, argparse
from Bio import SeqIO

def ParseArg():
    p=argparse.ArgumentParser( description = 'Remove head and/or tail sequences')
    p.add_argument('fin',type=str,metavar='input',help='Input fasta file')
    p.add_argument('output',type=str,metavar='output',help='Output fasta file')
    p.add_argument('-H','--head',type=int,default=0, help="Number of nucleotides cut at head. default: 0") 
    p.add_argument('-t','--tail',type=int,default=0, help="Number of nucleotides cut at tail. default: 0")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

args = ParseArg()

pair=open(args.fin,"r")
output=open(args.output,"w")


seq_file=SeqIO.parse(pair,"fasta")

i=1
for line in seq_file:
#  print line.description
#  print line.seq
    line.seq=line.seq[args.head:] if args.tail == 0 else line.seq[args.head:-args.tail]
    i+=1
    print >>output,">"+line.description
    print >>output,line.seq

print i
