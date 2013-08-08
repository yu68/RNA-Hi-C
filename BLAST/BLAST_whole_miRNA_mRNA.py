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
    
    p=argparse.ArgumentParser( description = 'DESCRIPTION: Run BLAST and select positive results from BLAST XML output for both miRNA and mRNA database', epilog='')
    p.add_argument("input",type=str,help="the input fastq/fasta file containing reads sequences")
    p.add_argument("-e","--evalue",dest="evalue",nargs='+',default=0,help="cutoff evalues, only choose alignment with evalue less than this cutoffs (need two, first for miRNA and second for mRNA, default: 1e-5, 1e-15).")
    p.add_argument("--mRNA_db",dest="mRNA_db",type=str,help="BLAST database of mRNA sequences",default="~/Stitch-seq/blast_db/mouse_RefSeq_mRNA.fasta")
    p.add_argument("--miRNA_db",dest="miRNA_db",type=str,help="BLAST database of miRNA sequences",default="~/Stitch-seq/blast_db/mouse_miRNA_mature.fa")
    p.add_argument("--blast_path",dest="blast_path",type=str,help="path for the local blast program",default="~/Softwares/ncbi-blast-2.2.27+/bin/blastn")
    p.add_argument("-o","--output",dest="output",type=str,help="output prefix. The program will generate three output file: prefix_all.txt - output all matches for reads; prefix_both.txt - output only reads with both miRNA and mRNA matches; prefix_stringent.txt - with requirement in both but also the one miRNA match need to be before mRNA matches in the read.")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def blast_align(fasta,blast_path,miRNA_db,mRNA_db):
    os.system(blast_path+" -task blastn -outfmt 5 -num_threads 6 -evalue 1e-3 -db "+miRNA_db+" -query "+fasta+" > "+args.output+"temp_blast_miRNA.xml")
    os.system(blast_path+" -task blastn -outfmt 5 -num_threads 6 -evalue 1e-5 -db "+mRNA_db+" -query "+fasta+" > "+args.output+"temp_blast_mRNA.xml")
    os.system("rm "+fasta)
    miRNA_records=NCBIXML.parse(open(args.output+"temp_blast_miRNA.xml"))
    mRNA_records=NCBIXML.parse(open(args.output+"temp_blast_mRNA.xml"))
    return (miRNA_records,mRNA_records)

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
    align_mi=0
    align_m=0    

    global args
    args=ParseArg()
    out_all=open(args.output+"_all.txt",'w')
    out_both=open(args.output+"_both.txt",'w')
    out_string=open(args.output+"_stringent.txt",'w')
    
    blast_path=args.blast_path
    miRNA_db=args.miRNA_db
    mRNA_db=args.mRNA_db    
    
    # E-values
    if args.evalue==0:
        evalue_mi=1e-5
        evalue_m=1e-15
    else:
        evalue_mi=float(args.evalue[0])
        evalue_m=float(args.evalue[1])
    
    # determine input sequence file type
    types="fastq"
    if args.input.split(".")[-1] in ["fa","fasta"]:
        types="fasta"

    seq_file=SeqIO.parse(args.input,types)
    
    ###################################
    ##    start editing from here    ## 
    ###################################
    
    for i, batch in enumerate(batch_iterator(seq_file, 100000)):
        t0=time()
        filename=args.output+"group_%i.fasta" % (i+1)
        handle=open(filename, "w")
        count=SeqIO.write(batch,handle,"fasta")
        handle.close()
        print "Wrote %i records to %s" % (count,filename)
        
        [miRNA_records,mRNA_records] = blast_align(filename,blast_path,miRNA_db,mRNA_db)
        print "BLAST aligned for %s ." % (filename)
        
        
        print "Start to parse BLAST results for %s" %(filename)
        for mi_record,m_record in itertools.izip(miRNA_records,mRNA_records):
            temp_output=''
            mi_indic=0 # whether there are miRNA alignment
            m_indic=0  # whether there are mRNA alignment
            mi_end=150  #shortest miRNA aligned end in query sequence
            n=n+1
            if (mi_record.query!=m_record.query):
                print >>sys.stderr,"The two query seqs from miRNA and mRNA results are not matched!"
                break
            query_name=mi_record.query.split('\t')[0]
            temp_output=query_name+'\n'
            for alignment in mi_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < evalue_mi:
                        mi_indic=1
                        line="\t".join (str(f) for f in [hsp.query_start,hsp.query_end,alignment.title,hsp.sbjct,hsp.sbjct_start,hsp.sbjct_end,hsp.expect,hsp.score])
                        print >>out_all,query_name+'\t'+line
                        temp_output=temp_output+line+'\n'
                        if mi_end>max(hsp.query_start,hsp.query_end):
                            mi_end=max(hsp.query_start,hsp.query_end)
        
            if mi_indic==0:
                mi_end=0
            temp_output2=temp_output
            for alignment in m_record.alignments:
                for hsp in alignment.hsps:
                    if (hsp.expect < evalue_m):
                        m_indic=1
                        line="\t".join (str(f) for f in [hsp.query_start,hsp.query_end,alignment.title,hsp.sbjct,hsp.sbjct_start,hsp.sbjct_end,hsp.expect,hsp.score])
                        print >>out_all,query_name+'\t'+line
                        temp_output=temp_output+line+'\n'
                        if (min(hsp.query_start,hsp.query_end)>mi_end):
                            temp_output2=temp_output2+line+'\n'
            if mi_indic+m_indic>=2:
                out_both.write(temp_output)
                if len(temp_output2.split("\n"))>=4:
                    out_string.write(temp_output2)
            if mi_indic==1:
                align_mi+=1
            if m_indic==1:
                align_m+=1
        print "After %s, got %i sequences, %i align to miRNA and %i align to mRNA." % (filename,n,align_mi,align_m)
        t1=time()
        print "Processing %s takes %f min. \n" %(filename,(t1-t0)/60)

    out_all.close()
    out_both.close()
    out_string.close()


if __name__=="__main__":
    main()
