#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Name:        Stitch-seq_Aligner
# Purpose:
#
# Author:      Pengfei
#
# Created:     21/08/2012
# Modified:    03/09/2013
# Copyright:   (c) Pengfei 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import sys, os, argparse
import pysam
import itertools,string
from Bio import SeqIO
from xplib import TableIO
from xplib import DBI
from xplib.Annotation import Bed
from Annotation import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def ParseArg():
    p=argparse.ArgumentParser( description = 'Align miRNA-mRNA pairs for Stitch-seq. print the alignable miRNA-mRNA pairs with coordinates', epilog = 'Library dependency: Bio, pysam, itertools')
    p.add_argument('input1',type=str,metavar='part1_reads',help='paired part1 fasta/fastq file')
    p.add_argument('input2',type=str,metavar='part2_reads',help='paired part2 fasta/fastq file')
    p.add_argument('bowtie_path',type=str,metavar='bowtie_path',help="path for the bowtie program")
    p.add_argument('blat_path', type=str, metavar='blat_path', help="path for the blat program (for miRNA only)")
    p.add_argument('-f','--fastq_to_fasta_path', dest='f2fpath', default='fastq_to_fasta', type=str, metavar='fastq_to_fasta_path', help="path for the fastq_to_fasta program (for miRNA only)")
    p.add_argument('-b','--bowtie2',action="store_true",help="set to use bowtie2 (--sensitive-local) for alignment")
    p.add_argument('-u','--unique',action="store_true",help="set to only allow unique alignment")
    p.add_argument('-s','--samtool_path',dest='spath', type=str,metavar='samtool_path',help="path for the samtool program",default='samtools')
    p.add_argument('--ref',type=str,metavar='rr',default="mm9",help="Reference seq for both parts (or the first part if --ref2 is specified), references will be used in the order specified",nargs="*")

    # notice that rna.fa file should be fetched from Ensembl with the following attributes specified in the following order:
    # Ensembl Transcript ID
    # cDNA sequences (this is the actual sequence and is not in header)
    # Transcript start (bp)
    # Transcript end (bp)
    # Strand
    # 5' UTR start
    # 5' UTR end
    # 3' UTR start
    # 3' UTR end
    # Associated gene name
    # Transcript biotype

    p.add_argument('--reftype',type=str,metavar='rt',default="genome",help="Reference types: (miRNA, genome, transcript, other)",nargs="*")
    p.add_argument('--ref2',type=str,metavar='rr',default="none",help="Reference seq for both parts (or the first part if --ref2 is specified), references will be used in the order specified",nargs="*")
    p.add_argument('--ref2type',type=str,metavar='rt',default="none",help="Reference types: (miRNA, genome, transcript, other)",nargs="*")
    p.add_argument('-a','--annotation',type=str,help='If specified, include the RNA type annotation for each aligned pair, need to give bed annotation RNA file')
    p.add_argument("-A","--annotationGenebed",dest="db_detail",type=str,help="annotation bed12 file for lincRNA and mRNA with intron and exon")
    p.add_argument("-nostr", "--ignore_strand", dest="nostr", action = "store_true", help="Reads mapped onto the wrong strand will be considered as not mapped by default. Set this flag to ignore strand information.")


    if len(sys.argv)==1:
        #print (p.print_help(),file=sys.stderr)
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

rev_table=string.maketrans('ACGTacgtN', 'TGCAtgcaN')
def revcomp(seq, rev_table):
    return seq.translate(rev_table)[::-1]

def blat_align(b_path, read, ref, fastqToFasta):
    # this is used for miRNA mapping only
    # b_path: blat path;
    # fastqToFasta: fastq_to_fasta path;
    # 1. use blat to map reads to miRNA references
    # 2. filter out the wrong ones:
    #   a) with wrong strand or
    #   b) not completely mapped to miRNA (have overlangs longer than 2 on both sides in total)
    # 3. annotate the correct entries with miRNA
    # 4. because blat does not generate unmapped reads, generate unmap file for the next ref

    # readfile = read.split("/")[-1].split(".")[0] + ".sam"
    tmp = read.split("/")
    filenametmp = tmp[-1].split(".")
    if not read.split(".")[-1] in ["fa","fasta"]:   # fastq needs to be converted
        fastafile = "/".join(tmp[:-1]) + "/" + ".".join(filenametmp[:-1]) + ".fasta"
        os.system(fastqToFasta + " -i " + read + " -o " + fastafile);
    else:
        fastafile = read

    output = "./" + "/".join(tmp[:-1]) + "/" + ".".join(filenametmp[:-1]) + ".blatresult"

    os.system(b_path+ " -stepSize=2 -minScore=15 -tileSize=6 "+ref+" "+read+" "+output)

    return output

def get_readdict(readfilename):
    seqdict = dict()

    fseq = open(readfilename, 'r')
    seqname = ''
    oldseq = ''

    for line in fseq:
        if line.startswith('>') or line.startswith('@'):
            if oldseq and seqname:
                seqdict[seqname] = oldseq
            seqname = line.strip('>@').split(" ")[0]
            oldseq = ''
        elif line.startswith('+'):
            # fastq
            if oldseq and seqname:
                seqdict[seqname] = oldseq
            seqname = ''
            oldseq = ''
        elif seqname:
            oldseq += line.strip()

    if oldseq and seqname:
        seqdict[seqname] = oldseq

    fseq.close()
    return seqdict

def blat_annotation(outputfilename, typename, readfilename, unmapfilename, anno = True, requireUnique = False, posstrand = True, strandenforced = False, mismatchthr = 2, results_dict = dict()):
    # results_dict is the dictionary for all previous results, can be none
    # key is the read name, value is the annotation string
    # this function put all qualified match in outputfilename
    # then put all unqualified match together with unmappedfile
    # return value is the new annotation dictionary
    newdict = dict()
    seqdict = get_readdict(readfilename)
    multiset = set()

    fres = open(outputfilename, 'r')
    # first populate the dictionary
    for x in xrange(5):
        fres.readline()

    for line in fres:
        tokens = line.split()
        strand = (tokens[8] == '+')
        if (not strandenforced) or strand == posstrand:
            # correct strand or strand not enforced
            if strand != posstrand:
                # incorrect strand, but not enforced
                # add the 21th column saying that strand is wrong
                strandcol = 'NonProperStrand'
            else:
                strandcol = 'ProperStrand'
            
            match = int(tokens[0])
            mismatch = int(tokens[1])
            readname = tokens[9]
            readlen = int(tokens[10])
            readstart = int(tokens[11])
            readend = int(tokens[12])

            mismatch += (readlen - readend) + readstart

            if mismatch <= mismatchthr:
                # this one is a match
                if readname in seqdict:
                    # first time happen
                    if anno:
                        curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seqdict[readname], typename, typename, tokens[13], '.']
                        if not strandenforced:
                            curr_anno_arr.append(strandcol)
                        curr_anno = '\t'.join(curr_anno_arr)
                    else:
                        curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seqdict[readname], typename]
                        if not strandenforced:
                            curr_anno_arr.append(strandcol)
                        curr_anno = '\t'.join(curr_anno_arr)
                    newdict[readname] = curr_anno
                    del seqdict[readname]
                elif (readname in newdict) or (readname in multiset):
                    # this is a multimatch
                    # and not first time
                    multiset.add(readname)
                    if requireUnique:
                        # then this multi-match will be discarded
                        del newdict[readname]
                    else:
                        if type(newdict[readname]) is str:
                            newdict[readname] = [newdict[readname]]
                        seq = newdict[readname][0].split()[4]
                        if anno:
                            curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seq, typename, typename, tokens[13], '.', strandcol]
                            newdict[readname].append('\t'.join(curr_anno_arr))
                        else:
                            curr_anno_arr = [tokens[13], tokens[15], tokens[16], tokens[8], seq, typename, strandcol]
                            newdict[readname].append('\t'.join(curr_anno_arr))

    fres.close()

    # merge annotation
    newanno = dict(results_dict.items() + newdict.items())

    # write unmapped reads
    funmap = open(unmapfilename, 'w')
    for key in seqdict:
        funmap.write(">" + key + os.linesep)
        funmap.write(seqdict[key] + os.linesep)

    funmap.close()

    return newanno

class ensemblSeq:

    def getSubType(self, start, end):
        # notice that this is relative to the RNA sequence
        # therefore, just compare directly with self.utrlength5 or self.utrlength3 should be enough
       
        # return format: [type, name, subtype]
        # first find the midpoint of the region
        midpoint = (start + end) / 2
        if midpoint < self.utrlength5:
            return "utr5"
        elif midpoint >= self.length - self.utrlength3:
            return "utr3"
        else:
            if self.biotype == "lincRNA" or self.biotype == "protein_coding":
                return "exon"
            elif "intron" in self.biotype:
                return "intron"
            else:
                return "."
    
    def getAnnotation(self, start, end):
        return [self.biotype, self.genename, self.getSubType(start, end)]

    def __init__(self, name, seq):
        # notice that the ">" in name is chopped
        self.longname = name.strip()
        tokens = name.split("|")
        self.transID = tokens[0]
        self.strand = (int(tokens[3]) > 0)
        self.tss = int(tokens[1]) if self.strand else int(tokens[2])
        self.tes = int(tokens[2]) if self.strand else int(tokens[1])
        self.genename = tokens[8]
        if not self.genename:
            self.genename = self.transID
        self.biotype = tokens[9]

        utrstarts5 = (tokens[4] if self.strand else tokens[5]).split(";")
        utrends5 = (tokens[5] if self.strand else tokens[4]).split(";")
        utrstarts3 = (tokens[6] if self.strand else tokens[7]).split(";")
        utrends3 = (tokens[7] if self.strand else tokens[6]).split(";")

        self.utrlength5 = 0
        self.utrlength3 = 0
        if tokens[4]:
            for i in xrange(len(utrstarts5)):
                self.utrlength5 += abs(int(utrends5[i]) - int(utrstarts5[i]))

        if tokens[6]:
            for i in xrange(len(utrstarts3)):
                self.utrlength3 += abs(int(utrends3[i]) - int(utrstarts3[i]))

        self.seq = seq.strip()
        self.length = len(self.seq)


def otherlib_annotation(outputbam, anno, readfilename, unmapfilename, annotationfile = 'misc', requireUnique = False, posstrand = True, strandenforced = False, results_dict = dict()):

    # outputbam is a pysam file
    # annotationfile can be:
    # "misc", "guess" from RNA name
    # or a known annotation file (for example, rna.fa from NCBI)

    newdict = dict()
    funmap = open(unmapfilename, 'w')

    if anno and annotationfile != 'misc':
        refdic = dict()
        # this will be a dictionary of ensembl entries
        fref = open(annotationfile, 'r')
        refname = ''
        refseq = ''
        refkey = ''

        for line in fref:
            if line.startswith(">"):
                if refkey:
                    # there is an old ref there
                    if refkey in refdic:
                        print sys.stderr, refkey
                    refdic[refkey] = ensemblSeq(refname, refseq)
                    refseq = ''
                    refkey = ''
                    refname = ''
                refkey = line.strip(">").split("|")[0]
                refname = line.strip(">").strip()
            else:
                refseq += line.strip()

        if refkey:
            if refkey in refdic:
                print sys.stderr, refkey
            refdic[refkey] = ensemblSeq(refname, refseq)
        fref.close()

    for record in outputbam:
        name = record.qname.split(" ")[0]
        try:
            record.opt('XS')
            unique = False
        except:
            unique = True
        if ((not requireUnique) or unique) and ((not strandenforced) or record.is_reverse != posstrand) and not record.is_unmapped:
            # this is a suitable entry
            strand = "+"
            if record.is_reverse:
                strand = "-"
            if record.is_reverse == posstrand:
                # strand is wrong
                strandcol = 'NonProperStrand'
            else:
                strandcol = 'ProperStrand'
            if anno:
                if annotationfile == 'misc':
                    rnaname = outputbam.getrname(record.tid).split(" ")[0]
                    rnaid = rnaname
                    if rnaname.startswith("NONMMUT"):
                        # This is a lncRNA from NONCODE
                        rnatype = "ncRNA"
                        rnasubtype = "."
                    else:
                        # tRNA
                        rnatype = "tRNA"
                        rnasubtype = "."
                else:
                    rnaensembl = refdic[outputbam.getrname(record.tid).split("|")[0]]
                    rnaid = rnaensembl.transID
                    [rnatype, rnaname, rnasubtype] = rnaensembl.getAnnotation(record.aend - record.alen + 1, record.aend)

                curr_anno_arr = (str(f) for f in [rnaid, record.aend - record.alen + 1, record.aend, strand, record.seq, annotationfile.split("/")[-1], rnatype, rnaname, rnasubtype, strandcol])
                newdict[name] = '\t'.join(curr_anno_arr)

            else:
                curr_anno_arr = (str(f) for f in [rnaid, record.aend - record.alen + 1, record.aend, strand, record.seq, annotationfile.split("/")[-1], strandcol])
                newdict[name] = '\t'.join(curr_anno_arr)

        elif record.is_unmapped:
            # write unmapped reads
            seq = record.seq
            if record.is_reverse:
                seq = revcomp(record.seq, rev_table)
            unmap_rec = SeqRecord(Seq(seq, IUPAC.unambiguous_dna), id = name, description='')
            SeqIO.write(unmap_rec, funmap, "fasta")

    funmap.close()
    
    newanno = dict(results_dict.items() + newdict.items())
    return newanno

def bowtie_align(b_path,read,ref,s_path,bowtie2):
    # b_path: bowtie path;
    # s_path: samtools path;
    # bowtie2: logic, true/false
    
    sam=read.split("/")[-1].split(".")[0]+".sam"
    if read.split(".")[-1] in ["fa","fasta"]:   # allow fasta and fastq for read
        foption=" -f"
    else:
        foption=""


    if ref.split(".")[-1] in ["fa","fasta"]:
        base=ref.split("/")[-1].split(".")[0]
        os.system("rm "+read.split("/")[-1]+".log")
        os.system(b_path+"-build "+ref+" "+base+" >> "+read+".log 2>&1")
        if not bowtie2:
            os.system(b_path+ foption+" -a --best --strata -n 1 -l 15 -e 200 -p 9 -S "+base+" "+read+" "+sam+" >> "+read.split("/")[-1]+".log 2>&1")
        else:
            os.system(b_path+ " -x "+base+foption+" -U "+read+" --sensitive-local -p 8 --reorder -t -S "+sam+" >> "+read.split("/")[-1]+".log 2>&1")
    else:
        os.system("rm "+read.split("/")[-1]+".log")
        if not bowtie2:
            os.system(b_path+ foption+" -a --best --strata -n 1 -l 15 -e 200 -p 9 -S "+ref+" "+read+" "+sam+" >> "+read.split("/")[-1]+".log 2>&1")
        else:
            os.system(b_path+ " -x "+ref+foption+" -U "+read+" --sensitive-local -p 8 --reorder -t -S "+sam+" >> "+read.split("/")[-1]+".log 2>&1")
    bam=read.split("/")[-1].split(".")[0]+".bam"
    os.system(s_path+ " view -Sb -o "+bam +" "+sam)
    os.system("rm "+sam)
    pysam.sort("-n",bam,"temp")
    align=pysam.Samfile("temp.bam","rb")
    os.system("rm temp.bam")
    os.system(s_path+ " sort "+bam+ " "+"sort_"+read.split("/")[-1].split(".")[0])
    os.system("rm "+bam)
    return align

def Included(record,RequireUnique):
    # record is a pysam read
    # non-unique alignment in Bowtie2 has 'XS' tag: https://www.biostars.org/p/18965/ 
    if RequireUnique:
        try:
            record.opt('XS')
            unique=False
        except:
            unique=True
    else:
        unique=True # not consider unique
    return (not record.is_unmapped)&unique

def genome_annotation(outputbam, annotationfile, detail, readfilename, unmapfilename, strandenforced = False, posstrand = True, requireUnique = False, results_dict = dict()):
    # annotationfile is annotation file
    # detail is db_detail file

    if annotationfile:
        dbi1=DBI.init(annotationfile,"bed")
        dbi2=DBI.init(detail,"bed")
        dbi3=DBI.init("/home/yu68/bharat-interaction/new_lincRNA_data/mouse.repeat.txt","bed")
    
    newdict = dict()
    funmap = open(unmapfilename, 'w')

    for record in outputbam:
        # print >> sys.stderr, record.qname
        IsMapped = False

        if Included(record, requireUnique):
            strand = ("+" if posstrand else "-")
            if record.is_reverse:
                strand = ("-" if posstrand else "+")
            if annotationfile:
                bed=Bed([outputbam.getrname(record.tid), record.pos, record.aend,'.',0.0,strand])
                [typ, name, subtype, strandcol] = annotation(bed,dbi1,dbi2,dbi3)
                if (not strandenforced) or strandcol == 'ProperStrand':
                    curr_anno_arr = (str(f) for f in [outputbam.getrname(record.tid), record.pos, record.aend, strand, record.seq, 'genome', typ, name, subtype, strandcol])
                    newdict[record.qname] = '\t'.join(curr_anno_arr)
                    IsMapped = True
            else:
                curr_anno_arr = (str(f) for f in [outputbam.getrname(record.tid), record.aend - record.alen + 1, record.aend, strand, record.seq, 'genome', strandcol])
                newdict[record.qname] = '\t'.join(curr_anno_arr)
                IsMapped = True

        if not IsMapped:
            # output all pairs that cannot be mapped on both sides as unmaped pairs into two fasta file
            seq = record.seq
            if record.is_reverse:
                seq = revcomp(record.seq, rev_table)
            unmap_rec = SeqRecord(Seq(seq, IUPAC.unambiguous_dna), id = record.qname, description='')
            SeqIO.write(unmap_rec, funmap, "fasta")
    
    funmap.close()
    
    newanno = dict(results_dict.items() + newdict.items())
    return newanno


def Main():
    args=ParseArg()

    reflist = [args.ref, args.ref]
    reftypelist = [args.reftype, args.reftype]
    readfilelist = [args.input1, args.input2]
    annodictlist = [dict(), dict()]
    strandslist = [True, False]
    strandenforced = True

    if args.nostr:
        strandenforced = False

    if (not type(args.ref2) is str) or args.ref2 != 'none':
        reflist[1] = args.ref2
        reftypelist[1] = args.ref2type

    for i in xrange(len(reflist)):
        reflistentry = reflist[i]
        reftypelistentry = reftypelist[i]
        rawreadfile = readfilelist[i]
        strand = strandslist[i]
        # mapping for every single ends
        # also put the annotations to the file
        # then merge all the annotations from every single step
        # then try to match all the annotations

        readfile = rawreadfile
        print >> sys.stderr, 'Mapping ' + str(i) + ' ' + rawreadfile + ' ...'

        for reffile, reftype in itertools.izip(reflistentry, reftypelistentry):
        
            # unmapped read file
            print >> sys.stderr, reftype
            tmp = readfile.split(".")
            unmap_read = ".".join(tmp[:-1])+"_unmap."+tmp[-1]

            if reftype.lower() == "mirna":
                outputfile = blat_align(args.blat_path, readfile, reffile, args.f2fpath)
                annodictlist[i] = blat_annotation(outputfile, reftype, readfile, unmap_read, args.annotation, args.unique, strand, strandenforced, 2, annodictlist[i])
            else:
                outputbam = bowtie_align(args.bowtie_path, readfile, reffile, args.spath, args.bowtie2)
                if reftype.lower() == "genome":
                    annodictlist[i] = genome_annotation(outputbam, args.annotation, args.db_detail, readfile, unmap_read, strandenforced, strand, args.unique, annodictlist[i])
                else:
                    if reftype.lower() == 'other':
                        annofile = 'misc'
                    else:
                        annofile = reffile.split('/')[-1]
                    annodictlist[i] = otherlib_annotation(outputbam, args.annotation, readfile, unmap_read, annofile, args.unique, strand, strandenforced, annodictlist[i])

            readfile = unmap_read

    for entry, value in annodictlist[0].iteritems():
        if entry in annodictlist[1]:
            # this is a pair
            if type(value) is str:
                if type(annodictlist[1][entry]) is str:
                    print '\t'.join(str(f) for f in [value, entry, annodictlist[1][entry],'OneToOne'])
                else:
                    for item2 in annodictlist[1][entry]:
                        print '\t'.join(str(f) for f in [value, entry, item2, 'OneToMany'])
            else:
                for item1 in value:
                    if type(annodictlist[1][entry]) is str:
                        print '\t'.join(str(f) for f in [item1, entry, annodictlist[1][entry], 'ManyToOne'])
                    else:
                        for item2 in annodictlist[1][entry]:
                            print '\t'.join(str(f) for f in [item1, entry, item2, 'ManyToMany'])
                    
if __name__ == '__main__':
    Main()




