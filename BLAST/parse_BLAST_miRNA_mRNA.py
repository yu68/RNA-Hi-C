import sys,os,argparse
from Bio.Blast import NCBIXML
import itertools

'''
dir(hsp):
['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand']
'''


def ParseArg():
    
    p=argparse.ArgumentParser( description = 'DESCRIPTION: select positive results from BLAST XML output for both miRNA and mRNA database', epilog='')
    p.add_argument("mi_xml",type=str,help="the input BLAST result xml file for miRNA db")
    p.add_argument("m_xml",type=str,help="the input BLAST result xml file for mRNA db")
    p.add_argument("-e","--evalue",dest="evalue",nargs='+',default=0,help="cutoff evalues, only choose alignment with evalue less than this cutoffs (need two, first for miRNA and second for mRNA, default: 1e-5, 1e-10).")
    p.add_argument("-o","--output",dest="output",type=str,help="output file contains sequenceswith both miRNA and mRNA alignment")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def main():
    #initialization
    n=0 # total number of query seq
    align_mi=0
    align_m=0    


    args=ParseArg()
    miRNA_result=open(args.mi_xml)
    mRNA_result=open(args.m_xml)
    miRNA_records=NCBIXML.parse(miRNA_result)
    mRNA_records=NCBIXML.parse(mRNA_result)
    output=open(args.output,'w')
    
    
    # E-values
    if args.evalue==0:
        evalue_mi=1e-5
        evalue_m=1e-15
    else:
        evalue_mi=float(args.evalue[0])
        evalue_m=float(args.evalue[1])
    
    for mi_record,m_record in itertools.izip(miRNA_records,mRNA_records):
        temp_output=''
        mi_indic=0 # whether there are miRNA alignment
        m_indic=0  # whether there are mRNA alignment
        mi_end=150  #shortest miRNA aligned end in query sequence
        n=n+1
        if (mi_record.query!=m_record.query):
            print >>sys.stderr,"The two query seqs from miRNA and mRNA results are not matched!"
            break
        temp_output=mi_record.query+'\n'
        for alignment in mi_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < evalue_mi:
                    mi_indic=1
                    line="\t".join (str(f) for f in [hsp.query_start,hsp.query_end,alignment.title,hsp.sbjct,hsp.sbjct_start,hsp.sbjct_end,hsp.expect,hsp.score])
                    temp_output=temp_output+line+'\n'
                    if mi_end>max(hsp.query_start,hsp.query_end):
                        mi_end=max(hsp.query_start,hsp.query_end)
        
        if mi_indic==0:
            mi_end=0

        for alignment in m_record.alignments:
            for hsp in alignment.hsps:
                if (hsp.expect < evalue_m) and (min(hsp.query_start,hsp.query_end)>mi_end):
                    m_indic=1
                    line="\t".join (str(f) for f in [hsp.query_start,hsp.query_end,alignment.title,hsp.sbjct,hsp.sbjct_start,hsp.sbjct_end,hsp.expect,hsp.score])
                    temp_output=temp_output+line+'\n'
        if mi_indic+m_indic>=2:
            output.write(temp_output)
        if mi_indic==1:
            align_mi+=1
        if m_indic==1:
            align_m+=1
    print n,align_mi,align_m

if __name__=="__main__":
    main()
