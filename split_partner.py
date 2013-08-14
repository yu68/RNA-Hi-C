import sys,os,argparse
from Bio.Blast import NCBIXML
import itertools
from Bio import SeqIO
from time import time


'''
dir(hsp):
['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand']
'''


def ParseArg():
    
    p=argparse.ArgumentParser( description = 'DESCRIPTION: Run BLAST, find linker sequences and split two parts connected by linkers', epilog='')
    p.add_argument("input",type=str,help="the input fasta file containing fragment sequences")
    p.add_argument("-e","--evalue",dest="evalue",type=float,default=0.00001,help="cutoff evalues, only choose alignment with evalue less than this cutoffs (default: 1e-5).")
    p.add_argument("--linker_db",dest="linker_db",type=str,help="BLAST database of linker sequences",default="~/Stitch-seq/blast_db/linker.fa")
    p.add_argument("--blast_path",dest="blast_path",type=str,help="path for the local blast program",default="~/Softwares/ncbi-blast-2.2.27+/bin/blastn")
    p.add_argument("-o","--output",dest="output",type=str,help="output file containing sequences of two sepatated parts")
    p.add_argument("-t","--trim",type=int,default=10,help="trim off the first this number of nt as index, default:10")
    p.add_argument("-b","--batch",type=int,default=100000,help="batch this number of fragments for BLAST at a time. default: 100000")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def blast_align(fasta,blast_path,linker_db):
    os.system(blast_path+" -task blastn -outfmt 5 -num_threads 6 -evalue 0.1 -db "+linker_db+" -query ./temp/"+fasta+" > "+name+"temp_blast_linker.xml")
    linker_records=NCBIXML.parse(open(name+"temp_blast_linker.xml"))
    os.system("rm ./temp/"+fasta)
    return (linker_records)

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch



def main():
    #initialization
    n=0 # total number of query seq
    align_no_linker=0 # aligned_reads
    os.system("mkdir temp") # create temp folder

    args=ParseArg()
    output=open(args.output,'w')
    trim_n=args.trim
    
    global name
    name=args.output.split(".")[0].split("_")[0]
  
    blast_path=args.blast_path
    linker_db=args.linker_db
    
    # E-values
    evalue=float(args.evalue)
    
    # determine input sequence file type
    types="fastq"
    if args.input.split(".")[-1] in ["fa","fasta"]:
        types="fasta"

    seq_file=SeqIO.parse(args.input,types)



    # output of different groups:
    output1=open("NoLinker_"+args.input,'w')
    output2=open("frontOnly_"+args.input,'w')
    output3=open("backOnly_"+args.input,'w')
    output4=open("Paired1_"+args.input,'w')
    output5=open("Paired2_"+args.input,'w')
    
    ###################################
    ##    start editing from here    ## 
    ###################################
    
    for i, batch in enumerate(batch_iterator(seq_file, args.batch)):
        t0=time()
        filename=name+"group_%i.fasta" % (i+1)
        handle=open("./temp/"+filename, "w")
        count=SeqIO.write(batch,handle,"fasta")
        handle.close()
        print "Wrote %i records to %s" % (count,filename)
        
        linker_records = blast_align(filename,blast_path,linker_db)
        print "BLAST aligned for %s." % (filename)
        

        print "Start to parse BLAST results for %s" %(filename)
        for linker_record, fragment in itertools.izip(linker_records,batch):
            n=n+1
            line=''
            start=len(fragment) # record the start of all linker alignment
            end=0 # record the end of all linker alignment
            pos={} # position of all alignment for one fragment
            counts={} # counts of linker alignment for one fragment
            Types="None"
            for alignment in linker_record.alignments:
                expect=evalue
                for hsp in alignment.hsps:
                    if hsp.expect < expect:
                        start=min(hsp.query_start,hsp.query_end,start)
                        end=max(hsp.query_end,hsp.query_start,end)
                        if alignment.hit_def in counts.keys():
                            counts[alignment.hit_def]+=1
                        else:
                            counts[alignment.hit_def]=1
                        pos[",".join (str(f) for f in [hsp.query_start,hsp.query_end])]=alignment.hit_def
            stat=";".join (f+':%d'%(counts[f]) for f in sorted(counts.iterkeys()))
            pos_order=";".join (str(f) for f in sorted(pos.iterkeys()))
            order=";".join (str(pos[f]) for f in sorted(pos.iterkeys()))

            if start>end:
                align_no_linker+=1
                if len(fragment)-trim_n>10:    
                    #XXXXXXXXXX11111111222222222 or XXXXXXXXXX11111111111111111 or XXXXXXXXXX2222222222222222
                    SeqIO.write(fragment[trim_n:],output1,types)
                    Types="NoLinker"
            elif (start>trim_n+10) and (end<len(fragment)-10):
                #XXXXXXXXXX111111LLLLLLL2222222
                SeqIO.write(fragment[trim_n:start],output4, types)
                SeqIO.write(fragment[end:],output5, types)
                Types="Paired"
            elif (start>trim_n+10):
                #XXXXXXXXXX1111111LLLLLLL
                SeqIO.write(fragment[trim_n:start],output2, types)
                Types="FrontOnly"
            elif (end<len(fragment)-10):
                #XXXXXXXXXXLLLLLL22222222
                SeqIO.write(fragment[end:],output3,types)
                Types="BackOnly"

            print >> output, "\t".join (str(f) for f in [fragment.id,stat,pos_order,Types,order])
        t1=time()
        print "After %s, got %i sequences, %i align to linkers. This batch takes %.2f min."%(filename,n,n-align_no_linker,(t1-t0)/60)

    output.close()


if __name__=="__main__":
    main()
