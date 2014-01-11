#!/usr/bin/env python
# Pengfei Yu yupf05@gmail.com

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
    
    p=argparse.ArgumentParser( description = 'DESCRIPTION: Run BLAST, find linker sequences and split two parts connected by linkers', epilog='Library dependency: Bio, itertools')
    p.add_argument("input",type=str,help="the input fasta file containing fragment sequences of type1 and type2")
    p.add_argument("type3_1",type=str,help="read_1 for evenlong (type3) fastq file")
    p.add_argument("type3_2",type=str,help="read_2 for evenlong (type3) fastq file")
    p.add_argument("-e","--evalue",dest="evalue",type=float,default=0.00001,help="cutoff evalues, only choose alignment with evalue less than this cutoffs (default: 1e-5).")
    p.add_argument("--linker_db",dest="linker_db",type=str,help="BLAST database of linker sequences",default="~/Stitch-seq/blast_db/linker.fa")
    p.add_argument("--blast_path",dest="blast_path",type=str,help="path for the local blast program",default="~/Softwares/ncbi-blast-2.2.27+/bin/blastn")
    p.add_argument("-o","--output",dest="output",type=str,help="output file containing sequences of two sepatated parts")
    p.add_argument("-t","--trim",type=int,default=10,help="trim off the first this number of nt as index, default:10")
    p.add_argument("-b","--batch",type=int,default=200000,help="batch this number of fragments for BLAST at a time. default: 100000")
    p.add_argument("-l","--length",type=int,default=15,help="shortest length to be considered for each part of the pair, default: 15")
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
            
    
    
def batch_iterator2(iterator1,iterator2,batch_size):
    """ batch iterator for paired-end files"""        
    entry1 = True #Make sure we loop once
    while entry1 :
        batch1 = []
        batch2 = []
        while len(batch1) < batch_size :
            try :
                entry1 = iterator1.next()
                entry2 = iterator2.next()
            except StopIteration :
                entry1 = None
            if entry1 is None :
                #End of file
                break
            batch1.append(entry1)
            batch2.append(entry2)
        if batch1 :
            yield batch1, batch2



def main():
    #initialization
    n=0 # total number of query seq
    align_no_linker=0 # aligned_reads
    os.system("mkdir temp") # create temp folder

    args=ParseArg()
    output=open(args.output,'w')
    trim_n=args.trim
    min_l=args.length
    
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
    
    # determine input type3 sequence files type
    types2="fastq"
    if args.type3_1.split(".")[-1] in ["fa","fasta"]:
        types2="fasta"
    
    seq_type3_1=SeqIO.parse(args.type3_1,types2)
    seq_type3_2=SeqIO.parse(args.type3_2,types2)

    n_p=0
    n_b=0
    n_f=0
    
    # output of different groups:
    output1=open("NoLinker_"+args.input,'w')
    output2=open("frontOnly_"+args.input,'w')
    output3=open("backOnly_"+args.input,'w')
    output4=open("Paired1_"+args.input,'w')
    output5=open("Paired2_"+args.input,'w')
    
    ###################################
    ##    start editing from here    ## 
    ###################################
    '''
    THis loop is for recovered fragment (type1&2)
    '''
    print "split partners for recovered fragments"
    for i, batch in enumerate(batch_iterator(seq_file, args.batch)):
        t0=time()
        filename=name+"group_%i.fasta" % (i+1)
        handle=open("./temp/"+filename, "w")
        count=SeqIO.write(batch,handle,"fasta")
        handle.close()
        #print "Wrote %i records to %s" % (count,filename)
        
        linker_records = blast_align(filename,blast_path,linker_db)
        #print "BLAST aligned for %s." % (filename)
        

        #print "Start to parse BLAST results for %s" %(filename)
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
                if len(fragment)-trim_n>min_l:    
                    #XXXXXXXXXX11111111222222222 or XXXXXXXXXX11111111111111111 or XXXXXXXXXX2222222222222222
                    SeqIO.write(fragment[trim_n:],output1,types)
                    Types="NoLinker"
            elif (start>trim_n+min_l) and (end<len(fragment)-min_l):
                #XXXXXXXXXX111111LLLLLLL2222222
                SeqIO.write(fragment[trim_n:start],output4, types)
                SeqIO.write(fragment[end:].reverse_complement(fragment.id),output5, types)
                n_p+=1
                Types="Paired"
            elif (start>trim_n+min_l):
                #XXXXXXXXXX1111111LLLLLLL
                SeqIO.write(fragment[trim_n:start],output2, types)
                n_f+=1
                Types="FrontOnly"
            elif (end<len(fragment)-min_l):
                #XXXXXXXXXXLLLLLL22222222
                SeqIO.write(fragment[end:].reverse_complement(fragment.id),output3,types)
                Types="BackOnly"
                n_b+=1   
            print >> output, "\t".join (str(f) for f in [fragment.id,stat,pos_order,Types,order])
        t1=time()
        print >>sys.stderr,"Type1-2, with_linker: (%i/%i). P:%i B:%i F:%i. Time: %.2f min.\r"%(n-align_no_linker,n,n_p,n_b,n_f,(t1-t0)/60),
    
    '''
    This loop is for evenlong(type3) reads to find "Paired"/"BackOnly"/"FrontOnly"    
    '''
    i=0    
    align_no_linker=0
    n=0
    print "\nsplit partners for type3 paired reads"
    for batch1, batch2 in batch_iterator2(seq_type3_1, seq_type3_2, args.batch):
        t0=time()
        filename1="type3_group_%i_1.fasta" % (i+1)
        filename2="type3_group_%i_2.fasta" % (i+1)
        handle1=open("./temp/"+filename1, "w")
        count1=SeqIO.write(batch1,handle1,"fasta")
        handle1.close()
        handle2=open("./temp/"+filename2, "w")
        count2=SeqIO.write(batch2,handle2,"fasta")
        handle2.close()
        #if count1==count2:
            #print "Wrote %i records to %s and %s" % (count1,filename1,filename2)
        
        linker_records1 = blast_align(filename1,blast_path,linker_db)
        linker_records2 = blast_align(filename2,blast_path,linker_db)
        #print "BLAST aligned for %s.and %s" % (filename1, filename2)
        i=i+1
        j=0
        for linker_record1, linker_record2 in itertools.izip(linker_records1,linker_records2):
            n=n+1
            read1=batch1[j]
            read2=batch2[j]
            j+=1
            if read1.id!=read2.id:
                print "ERROR!! ID not match for paired type3 reads"
                sys.exit(0)
            start1=len(read1)
            for alignment in linker_record1.alignments:
                expect=evalue
                for hsp in alignment.hsps:
                    if hsp.expect < expect:
                        start1=min(hsp.query_start,hsp.query_end,start1)
            start2=len(read2)
            for alignment in linker_record1.alignments:
                expect=evalue
                for hsp in alignment.hsps:
                    if hsp.expect < expect:
                        start2=min(hsp.query_start,hsp.query_end,start2)
          
            if (start1==len(read1) and start2==len(read2)):
                align_no_linker+=1
            if (start1>trim_n+min_l) and (start2>min_l):
                SeqIO.write(read1[trim_n:start1],output4, types)
                SeqIO.write(read2[:start2],output5, types)
                n_p+=1
                Types="Paired"
            elif start1>trim_n+min_l:
                SeqIO.write(read1[trim_n:start1],output2, types)
                n_f+=1
                Types="FrontOnly"
            elif start2>min_l:
                n_b+=1
                SeqIO.write(read2[:start2],output3,types)
                Types="BackOnly"
        t1=time()
        print >>sys.stderr,"Type3, with_linker: (%i/%i). P:%i B:%i F:%i. Time: %.2f min.\r"%(n-align_no_linker,n,n_p,n_b,n_f,(t1-t0)/60),
                

            

    output.close()


if __name__=="__main__":
    main()
