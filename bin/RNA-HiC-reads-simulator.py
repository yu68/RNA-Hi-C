import sys, os, argparse
import string, copy
from xplib import TableIO
from xplib.Annotation import Bed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
from Bio.SeqRecord import SeqRecord
import random,time
from scipy import stats
from xplib import DBI
from Annotation import *

def ParseArg():
    p=argparse.ArgumentParser( description = 'generate simulated sequences', epilog = 'Library dependency: Bio, itertools')
    p.add_argument('linker',type=str,metavar='linker',help='file contianing linker sequences')
    p.add_argument('barcode',type=str,metavar='barcode',help='file contianing barcode sequences')
    p.add_argument('-n','--number',type=int,dest='num',help='number of simulated sequences,default:10000',default=10000)
    p.add_argument('-p','--parameter',nargs='+',type=float,help='parameters for percentage of LinkerOnly, Nolinker, RNA1-linker, linker-RNA2, RNA1-linker-RNA2,default:[.1,.3,.1,.3,.2]',default=[.1,.3,.1,.3,.2])
    p.add_argument('-l','--length',type=int,dest='len',help='length of generated reads',default=100)
    p.add_argument('-g',"--genomeFa",type=str,default="/home/yu68/Software/bowtie-0.12.7/indexes/mm9.fa",help="genomic sequence,need to be fadix-ed")
    p.add_argument('-s','--samtool_path',dest='spath', type=str,metavar='samtool_path',help="path for the samtool program",default='samtools')
    p.add_argument('-a','--annotation',type=str,help='If specified, include the RNA type annotation for each aligned pair, need to give bed annotation RNA file')
    p.add_argument("-A","--annotationGenebed",dest="db_detail",type=str,help="annotation bed12 file for lincRNA and mRNA with intron and exon")
    p.add_argument("-o","--output",type=str,help="output file for information of each simulated read pairs, including fragment_len, classes, linker_num, and information for RNA1 and RNA2")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def selectType(length):
    ''' select RNA types based on a length-related distribution  '''
    Types = ["miRNA","protein_coding","lincRNA","snoRNA","snRNA","tRNA"]
    if length<50:
        pk = [.2,.2,.1,.2,.2,.1]
    else:
        pk = [.05,.4,.2,.2,.1,.05]
    xk = range(6)
    custm = stats.rv_discrete(name='custm', values=(xk, pk))
    Type = Types[custm.rvs(size=1)]
    return Type

def randRegion(length,RNAs):
    '''  get random sequences coordinate given desire seq length and a database of all RNAs'''
    Type = selectType(length)
    RNA = random.choice(RNAs[Type])
    start = random.randrange(RNA.start,max(RNA.stop-length,RNA.start)+1)
    return Bed([RNA.chr,start,start+length],strand=RNA.strand),Type

rev_table=string.maketrans('ACGTacgt', 'TGCAtgca')
def revcomp(seq, rev_table):
    return seq.translate(rev_table)[::-1]

def fetchSeq(chro,start,end,strand,fasta,s_path):
    ''' s_path is the path of samtools  '''
    region = '%s:%d-%d'%(chro,start,end)
    seq="".join(os.popen(s_path+" faidx "+fasta+" "+region).read().split('\n')[1:])
    if strand=="-":
        seq = revcomp(seq,rev_table)
    return seq      

def generatePairs(fragment,read_len):
    ''' generate pair-end read pairs based on fragment sequences and read length  '''
    P7 = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGAT" 
    P5 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCT" 
    fragment_RC = str(Seq(fragment,IUPAC.unambiguous_dna).reverse_complement())
    if len(fragment)>read_len:
        read1 = fragment[:read_len]
        read2 = fragment_RC[:read_len]
    elif len(fragment)+2*len(P7)>=read_len:
        read1 = (fragment+P7*2)[:read_len]
        read2 = (fragment_RC+P5*2)[:read_len]
    else:
        N_num = read_len - len(fragment) - 2*len(P7)
        read1 = fragment+P7*2+'N'*N_num
        read2 = fragment+P5*2+'N'*N_num
    return read1,read2 
    

def Main():
    args=ParseArg()
    fastq1=open("simulated_"+str(args.num)+"_read_R1.fastq","w")
    fastq2=open("simulated_"+str(args.num)+"_read_R2.fastq","w")


    RNA=TableIO.parse(args.annotation,'bed')

    # create a dictionary for all RNAs pools except rRNA
    RNAs = {}
    for b in RNA:
        if b.id.startswith('rRNA'): continue
        if b.chr.startswith('chrM') or b.chr.startswith('chrNT'): continue
        Type = b.id.split(".")[0]
        if Type in RNAs:
            RNAs[Type].append(b)
        else:
            RNAs[Type]=[b]
    
    #---------------- read linker seq ------------------
    linkers=[]
    for i in open(args.linker,'r'):
        i=i.strip()
        linkers.append(i)
    #---------------------------------------------------

    #---------------- read barcode ---------------------
    barcodes=[]
    for i in open(args.barcode,'r'):
        i=i.strip()
        barcodes.append(i)
    #---------------------------------------------------
    
 
    # sample different classes: LinkerOnly, Nolinker, RNA1-linker, linker-RNA2, RNA1-linker-RNA2
    xk = range(5)
    pk = args.parameter
    custm = stats.rv_discrete(name='custm', values=(xk, pk))
    Class_index = custm.rvs(size=args.num)
    # specify output
    out = open(args.output,'w')


    # initiate the annotation database
    if args.db_detail:
        print >> sys.stderr, " # Index for the annotation database"
        dbi1=DBI.init(args.annotation,"bed")
        dbi2=DBI.init(args.db_detail,"bed")
        dbi3=DBI.init("/home/yu68/bharat-interaction/new_lincRNA_data/mouse.repeat.txt","bed")

    print >> sys.stderr, " # Start to simulate reads"
    t0 = time.time()
    for i in range(0,args.num):
        pair_id = "read_"+str(i)
        # barcode
        randSeq = "".join([random.choice("ACGT") for x in range(6)])
        barcode = randSeq[0:4]+barcodes[0]+randSeq[4:6]

        index = Class_index[i]  # index for different classes of fragments
        # Sample RNA1 and RNA2
        RNA1_len = random.randrange(15,150)
        b,Type = randRegion(RNA1_len,RNAs)
        RNA1_seq = fetchSeq(b.chr,b.start,b.stop,b.strand,args.genomeFa,args.spath)
        if args.db_detail:
            [name1,typ1,subtype1]=annotation(b,dbi1,dbi2,dbi3)
            RNA1_str = "\t".join(str(f) for f in [b.chr,b.start,b.stop,b.strand,name1,typ1,subtype1])
        else:
            RNA1_str = "\t".join(str(f) for f in [b.chr,b.start,b.stop,b.strand,Type])
        RNA2_len = random.randrange(15,150)
        b,Type = randRegion(RNA2_len,RNAs)
        RNA2_seq = fetchSeq(b.chr,b.start,b.stop,b.strand,args.genomeFa,args.spath)
        if args.db_detail:
            [name2,typ2,subtype2]=annotation(b,dbi1,dbi2,dbi3)
            RNA2_str = "\t".join(str(f) for f in [b.chr,b.start,b.stop,b.strand,name2,typ2,subtype2])
        else:
            RNA2_str = "\t".join(str(f) for f in [b.chr,b.start,b.stop,b.strand,Type])
       
        # fragment is the recovered cDNA fragment       
        if index == 1:  # single RNA or RNA1-RNA2
            if random.choice([0,1])==0:  # single RNAs
                fragment = barcode+RNA1_seq+RNA2_seq  
                print >> out, pair_id+"\t%d\tRNA1-RNA2\t0"%(len(fragment))+"\t"+RNA1_str+'\t'+RNA2_str
            else:  # RNA1-RNA2
                fragment = barcode+RNA1_seq
                print >> out, pair_id+"\t%d\tsingleRNA\t0"%(len(fragment))+"\t"+RNA1_str
        else:
            linker_n = random.choice([1,2])  # number of linkers in fragment
            linker = "".join([linkers[0]]*linker_n)
            if index == 0:
                fragment = barcode+linker
                print >> out, pair_id+"\t%d\tlinkerOnly\t%d"%(len(fragment),linker_n)
            elif index == 2:
                fragment = barcode+RNA1_seq+linker
                print >> out, pair_id+"\t%d\tRNA1-linker\t%d"%(len(fragment),linker_n)+"\t"+RNA1_str
            elif index == 3:
                fragment = barcode+linker+RNA2_seq
                print >> out, pair_id+"\t%d\tlinker-RNA2\t%d"%(len(fragment),linker_n)+"\t"+RNA2_str
            elif index == 4:
                fragment = barcode+RNA1_seq+linker+RNA2_seq
                print >> out, pair_id+"\t%d\tRNA1-linker-RNA2\t%d"%(len(fragment),linker_n)+"\t"+RNA1_str+"\t"+RNA2_str

        read1,read2 = generatePairs(fragment,args.len)
        score=[]
        for j in range(0, args.len):
            score.append(random.randrange(10,40))
        record1 = SeqRecord(Seq(read1,generic_dna),id=pair_id)
        record1.letter_annotations["phred_quality"] = score
        record2 = SeqRecord(Seq(read2,generic_dna),id=pair_id)
        record2.letter_annotations["phred_quality"] = score
        SeqIO.write(record1,fastq1,"fastq")
        SeqIO.write(record2,fastq2,"fastq")
     
        if i%100==0:
            print >>sys.stderr, "generate pairs %d\r"%(i),
    fastq1.close()
    fastq2.close()
    out.close()
    print time.time()-t0

if __name__ == '__main__':
    Main()

