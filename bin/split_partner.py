#!/usr/bin/env python
# Pengfei Yu yupf05@gmail.com

import sys,os,argparse
from Bio.Blast import NCBIXML
import itertools
from Bio import SeqIO
from time import time
#import pp
import threading
from Queue import Queue

# lock=thread.allocate_lock()
# lock1=thread.allocate_lock()
# lock2=thread.allocate_lock()
# lock3=thread.allocate_lock()
# lock4=thread.allocate_lock()
# lock5=thread.allocate_lock()

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
    p.add_argument("-d","--evalue_end",dest="evalue_end",type=float,default=0.001,help="cutoff evalues at end, only choose alignment at the end with evalue less than this cutoffs (default: 1e-3).")
    p.add_argument("--linker_db",dest="linker_db",type=str,help="BLAST database of linker sequences",default="/home/yu68/Stitch-seq/blast_db/linker.fa")
    p.add_argument("--blast_path",dest="blast_path",type=str,help="path for the local blast program",default="/home/yu68/Softwares/ncbi-blast-2.2.27+/bin/blastn")
    p.add_argument("-o","--output",dest="output",type=str,help="output file containing sequences of two sepatated parts")
    p.add_argument("-t","--trim",type=int,default=10,help="trim off the first this number of nt as index, default:10")
    p.add_argument("-b","--batch",type=int,default=200000,help="batch this number of fragments for BLAST at a time. default: 200000")
    p.add_argument("-r","--release",action='store_true',help="set to allow released criterion for Paired fragment in Type 3, include those ones with no linker in two reads")
    p.add_argument("-l","--length",type=int,default=15,help="shortest length to be considered for each part of the pair, default: 15")
    p.add_argument("-k","--linker_trim_length",type=int,default=10,help="linker length to be truncated for Paired fragment in Type 3, default: 10")
    p.add_argument("-p","--parallel",type=int,default=5,help="Number of threads")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def blast_align(fasta,blast_path,linker_db):
    fasta_name=fasta.split(".")[0]
    os.system(blast_path+" -task blastn -outfmt 5 -num_threads 6 -evalue 0.1 -db "+linker_db+" -query ./temp/"+fasta+" > ./temp/"+fasta_name+"_blast_linker.xml")
    linker_records=NCBIXML.parse(open("./temp/"+fasta_name+"_blast_linker.xml"))
#    os.system("rm ./temp/"+fasta)
#    os.system("rm ./temp/"+fasta_name+"_blast_linker.xml")
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

def RecoverType12(filename,files):
#    print filename
    # fasta_name=filename.split(".")[0]
    # os.system(blast_path+" -task blastn -outfmt 5 -num_threads 6 -evalue 0.1 -db "+linker_db+" -query ./temp/"+filename+" > ./temp/"+fasta_name+"_blast_linker.xml")
    # linker_records=NCBIXML.parse(open("./temp/"+fasta_name+"_blast_linker.xml"))
    # os.system("rm ./temp/"+fasta)
    linker_records = blast_align(filename,blast_path,linker_db)
        #print "BLAST aligned for %s." % (filename)
        

        #print "Start to parse BLAST results for %s" %(filename)
    batch_temp=SeqIO.parse(open('./temp/'+filename,"rU"),types)

    for linker_record, fragment in itertools.izip(linker_records,batch_temp):
    #    print fragment
    #    print type(fragment)
    #    print len(fragment)
        files.Add("n")
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
                    start=min(hsp.query_start - 1, hsp.query_end - 1, start)
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
            files.Add("no_linker")
            if len(fragment)-trim_n>min_l:    
                #XXXXXXXXXX11111111222222222 or XXXXXXXXXX11111111111111111 or XXXXXXXXXX2222222222222222
                files.Output_Type(1,types,fragment[trim_n:])
                # g_args.lock1.acquire()
                # SeqIO.write(fragment[trim_n:],g_args.output1,types)
                # g_args.lock1.release()
                Types="NoLinker"
        elif (start>trim_n+min_l) and (end<len(fragment)-min_l):
            #XXXXXXXXXX111111LLLLLLL2222222
            files.lock_evenlong.acquire()
            files.Output_Type(4,types,fragment[trim_n:start])
            # g_args.lock4.acquire()
            # SeqIO.write(fragment[trim_n:start],g_args.output4, types)
            # g_args.lock4.release()
            files.Output_Type(5,types,fragment[end:].reverse_complement(fragment.id))
            files.lock_evenlong.release()
            # g_args.lock5.acquire()
            # SeqIO.write(fragment[end:].reverse_complement(fragment.id),g_args.output5, types)
            # g_args.lock5.release()
   #         n_p+=1
            Types="Paired"
        elif (start>trim_n+min_l):
            #XXXXXXXXXX1111111LLLLLLL
            files.Output_Type(2,types,fragment[trim_n:start])
            # g_args.lock2.acquire()
            # SeqIO.write(fragment[trim_n:start],g_args.output2, types)
            # g_args.lock2.release()
    #        n_f+=1
            Types="FrontOnly"
        elif (end<len(fragment)-min_l):
            #XXXXXXXXXXLLLLLL22222222
            files.Output_Type(3,types,fragment[end:].reverse_complement(fragment.id))
            # g_args.lock3.acquire()
            # SeqIO.write(fragment[end:].reverse_complement(fragment.id),g_args.output3,types)
            # g_args.lock3.release()
            Types="BackOnly"
    #        n_b+=1   

        files.Output_detail([fragment.id,stat,pos_order,Types,order])
        # g_args.lock.acquire()
        # print >> g_args.output, "\t".join (str(f) for f in [fragment.id,stat,pos_order,Types,order])
        # g_args.lock.release()


def RecoverType3(filename1,filename2,files):
#    print filename1
#    print filename2
    linker_records1 = blast_align(filename1,blast_path,linker_db)
    linker_records2 = blast_align(filename2,blast_path,linker_db)

    #print "BLAST aligned for %s.and %s" % (filename1, filename2)
    j=0

    batch_temp1=SeqIO.parse(open('./temp/'+filename1,"rU"),"fasta")
    batch_temp2=SeqIO.parse(open('./temp/'+filename2,"rU"),"fasta")
 #   print batch_temp1
    for linker_record1, linker_record2 in itertools.izip(linker_records1,linker_records2):
        files.Add("n")
        read1=batch_temp1.next()
        read2=batch_temp2.next()
 #       print read1,read2

        j+=1
        if read1.id!=read2.id:
            print "ERROR!! ID not match for paired type3 reads"
            sys.exit(0)
        start1=len(read1)
        for alignment in linker_record1.alignments:
            expect=evalue
            for hsp in alignment.hsps:
                if max(hsp.query_start, hsp.query_end) >= len(read1):
                    expect = evalueEnd
                else:
                    expect = evalue
                if hsp.expect < expect:
                    start1=min(hsp.query_start - 1, hsp.query_end - 1, start1)
        start2=len(read2)
        for alignment in linker_record2.alignments:
            expect=evalue
            for hsp in alignment.hsps:
                if max(hsp.query_start, hsp.query_end) >= len(read2):
                    expect = evalueEnd
                else:
                    expect = evalue
                if hsp.expect < expect:
                    start2=min(hsp.query_start - 1, hsp.query_end - 1, start2)
      
        if (start1==len(read1) and start2==len(read2)):
            files.Add("no_linker")
            if not release: continue
            start1 = len(read1) - linkerTrimLen
            start2 = len(read2) - linkerTrimLen
        if (start1>trim_n+min_l) and (start2>min_l):
            files.lock_evenlong.acquire()
            files.Output_Type(4,types,read1[trim_n:start1])
#            SeqIO.write(read1[trim_n:start1],output4, types)
            files.Output_Type(5,types,read2[:start2])
            files.lock_evenlong.release()
#            SeqIO.write(read2[:start2],output5, types)
 #           n_p+=1
            Types="Paired"
        elif start1>trim_n+min_l:
            files.Output_Type(2,types,read1[trim_n:start1])
#            SeqIO.write(read1[trim_n:start1],output2, types)
 #           n_f+=1
            Types="FrontOnly"
        elif start2>min_l:
#            n_b+=1
            files.Output_Type(3,types,read2[:start2])
#            SeqIO.write(read2[:start2],output3,types)
            Types="BackOnly"

#    print >>sys.stderr,"Type3, with_linker: (%i/%i). P:%i B:%i F:%i. Time: %.2f min.\r"%(n-align_no_linker,n,n_p,n_b,n_f,(t1-t0)/60),


class OutputFiles():
    def __init__(self,Input,output,**kwargs):
        self.CountType={2:0,3:0,4:0,"n":0,"no_linker":0}
        self.lock=threading.Lock()
        self.lock1=threading.Lock()
        self.lock2=threading.Lock()
        self.lock3=threading.Lock()
        self.lock_evenlong=threading.Lock()
        self.lock_n=threading.Lock()
        self.lock_noLinker=threading.Lock()
        # self.trim_n=Args.trim
        # self.min_l=Args.length
        # self.blast_path=Args.blast_path
        # self.linker_db=Args.linker_db
        # self.evalue=float(Args.evalue)
        self.output=open(output,'w')
        self.output1=open("NoLinker_"+Input,'w')
        self.output2=open("frontOnly_"+Input,'w')
        self.output3=open("backOnly_"+Input,'w')
        self.output4=open("Paired1_"+Input,'w')
        self.output5=open("Paired2_"+Input,'w')
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])

    def Output_Type(self,f_n,type,read):
        output_dic={1:self.output1,2:self.output2,3:self.output3,4:self.output4,5:self.output5}
        lock_dic={1:self.lock1,2:self.lock2,3:self.lock3}
        if f_n in [1,2,3]:
            lock_dic[f_n].acquire()
            SeqIO.write(read,output_dic[f_n],type)
            if f_n in self.CountType:
                self.CountType[f_n]+=1
            lock_dic[f_n].release()
        elif f_n in [4,5]:
            SeqIO.write(read,output_dic[f_n],type)
            if f_n in self.CountType:
                self.CountType[f_n]+=1

    def Add(self,t):
        lock_dic={"n":self.lock_n,"no_linker":self.lock_noLinker}
        if t in self.CountType:
            lock_dic[t].acquire()
            self.CountType[t]+=1
            lock_dic[t].release()
        else:
            print >>sys.stderr,"ERROR!"


    def Output_detail(self,l):
        self.lock.acquire()
        print >> self.output, "\t".join (str(f) for f in l)
        self.lock.release()

def do_stuff(q,obj):
    while not q.empty():
#        try:
        item=q.get()
        RecoverType12(item,obj)
#       do_work(obj,item)
        q.task_done()
 #       except:
 #           print >>sys.stderr,"Error!"
 #           os._exit(0)

def do_stuff2(q,obj):
    while not q.empty():
#        try:
        item=q.get()
        RecoverType3(item[0],item[1],obj)
#       do_work(obj,item)
        q.task_done()
 #       except:
 #           print >>sys.stderr,"error!"
 #           os._exit(0)


def main():
    #initialization
#    n=0 # total number of query seq
#    align_no_linker=0 # aligned_reads
    os.system("mkdir temp") # create temp folder
    
    args=ParseArg()
#    output=open(args.output,'w')
    global trim_n,min_l
    trim_n=args.trim
    min_l=args.length
    
    name=args.output.split(".")[0].split("_")[0]
    global blast_path,linker_db
    blast_path=args.blast_path
    linker_db=args.linker_db

    # create blast database if not yet
    if not os.path.isfile(linker_db+".nhr"):
        print >> sys.stderr, "  ## blast database not detected, making blast db for seq"
        blastdir,_ = os.path.split(blast_path)
        os.system(blastdir+"/makeblastdb -in "+linker_db+" -dbtype 'nucl' -title "+os.path.splitext(linker_db)[0])

    
    # E-values
    global evalue
    evalue=float(args.evalue)

    # End E-values
    global evalueEnd
    evalueEnd = float(args.evalue_end)
    
    # Linker trim length
    global linkerTrimLen
    linkerTrimLen = int(args.linker_trim_length)

    #release
    global release
    release=args.release
  
    # determine input sequence file type
    global types
    types="fastq"
    if args.input.split(".")[-1] in ["fa","fasta"]:
        types="fasta"

    seq_file=SeqIO.parse(args.input,types)
    
    # determine input type3 sequence files type
    global types2
    types2="fastq"
    if args.type3_1.split(".")[-1] in ["fa","fasta"]:
        types2="fasta"
    
    seq_type3_1=SeqIO.parse(args.type3_1,types2)
    seq_type3_2=SeqIO.parse(args.type3_2,types2)

    Files=OutputFiles(args.input,args.output)
    num_thread=args.parallel

    # n_p=0
    # n_b=0
    # n_f=0

    # ncpus=args.parallel
    # ppservers = ()
    # job_server = pp.Server(ncpus, ppservers=ppservers)
  #  print job_server

    
    ###################################
    ##    start editing from here    ## 
    ###################################
    '''
    THis loop is for recovered fragment (type1&2)
    '''
    print "split partners for recovered fragments"
    t0=time()
    num_thread=args.parallel
    q=Queue(maxsize=0)
    for i, batch in enumerate(batch_iterator(seq_file, args.batch)):
        filename=name+"group_%i.fasta" % (i+1)
        handle=open("./temp/"+filename, "w")
        count=SeqIO.write(batch,handle,"fasta")
        handle.close()
        #print "Wrote %i records to %s" % (count,filename)
        q.put(filename)
#        RecoverType12(filename,batch,types,G_Arg)
  #      jobs.append(job_server.submit(RecoverType12,(filename,batch,types,G_Arg),(NCBIXML,SeqIO,),("itertools",),globals=globals()))
  
#    for j in jobs:
#        j()
    for i in range(num_thread):
        worker=threading.Thread(target=do_stuff,args=(q,Files,))
        worker.setDaemon(True)
        worker.start()

    q.join()
    
    t1=time()
    print >>sys.stderr,"Type1-2, with_linker: (%i/%i). P:%i B:%i F:%i. Time: %.2f min."%(Files.CountType["n"]-Files.CountType["no_linker"],Files.CountType["n"],Files.CountType[4],Files.CountType[3],Files.CountType[2],(t1-t0)/60) 

    '''
    This loop is for evenlong(type3) reads to find "Paired"/"BackOnly"/"FrontOnly"    
    '''
    i=0
    for counttype in Files.CountType.keys():
        Files.CountType[counttype]=0
 #   align_no_linker=0
 #   n=0
    print "\nsplit partners for type3 paired reads"
    t0=time()
    q2=Queue(maxsize=0)
    for batch1, batch2 in batch_iterator2(seq_type3_1, seq_type3_2, args.batch):
        filename1="type3_group_%i_1.fasta" % (i+1)
        filename2="type3_group_%i_2.fasta" % (i+1)
        handle1=open("./temp/"+filename1, "w")
        count1=SeqIO.write(batch1,handle1,"fasta")
        handle1.close()
        handle2=open("./temp/"+filename2, "w")
        count2=SeqIO.write(batch2,handle2,"fasta")
        handle2.close()
        i=i+1
        q2.put([filename1,filename2])

        #if count1==count2:
            #print "Wrote %i records to %s and %s" % (count1,filename1,filename2)
    for i in range(num_thread):
        worker=threading.Thread(target=do_stuff2,args=(q2,Files,))
        worker.setDaemon(True)
        worker.start()

    q2.join() 
    t1=time()
    print >>sys.stderr,"Type3, with_linker: (%i/%i). P:%i B:%i F:%i. Time: %.2f min."%(Files.CountType["n"]-Files.CountType["no_linker"],Files.CountType["n"],Files.CountType[4],Files.CountType[3],Files.CountType[2],(t1-t0)/60) 

    os.system("rm -r ./temp")


if __name__=="__main__":
    main()
