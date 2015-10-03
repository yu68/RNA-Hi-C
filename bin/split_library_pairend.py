#!/usr/bin/env python
# Pengfei Yu yupf05@gmail.com

import sys,os,argparse
from Bio import SeqIO
from Bio.SeqIO import QualityIO
import argparse
import threading
#from Queue import Queue
#import yappi
#import pstats

#def enable_thread_profiling():
#    '''Monkey-patch Thread.run to enable global profiling.
        
#    Each thread creates a local profiler; statistics are pooled
#    to the global stats object on run completion.'''
#    threading.Thread.stats = None
#    thread_run = threading.Thread.run
    
#    def profile_run(self):
#        self._prof = cProfile.Profile()
#        self._prof.enable()
#        thread_run(self)
#        self._prof.disable()
        
#        if threading.Thread.stats is None:
#            threading.Thread.stats = pstats.Stats(self._prof)
#        else:
#            threading.Thread.stats.add(self._prof)
            
#    threading.Thread.run = profile_run
    
def get_thread_stats():
    stats = getattr(threading.Thread, 'stats', None)
    if stats is None:
        raise ValueError, ('Thread profiling was not enabled,'
                'or no threads finished running.')
    return stats

def ParseArg():
    '''Parse the argument'''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -q example.F1.fastq example.R1.fastq -b barcode.txt', epilog = 'Library dependency: Bio')
    group=p.add_mutually_exclusive_group()
    group.add_argument("-f","--fasta",action='store_true',help='add this option for fasta input file')
    group.add_argument("-q","--fastq",action='store_true',help='add this option for fastq input file')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('input1',type=str,help='input fastq/fasta file 1 for pairend data (contain barcodes)')
    p.add_argument('input2',type=str,help='input fastq/fasta file 2 for pairend data')
    p.add_argument('-b','--barcode',dest='barcode',type=str,help='barcode file')
    p.add_argument('-r','--range',dest='range',nargs='+',default=0,help='set range for barcode location within reads,default is full read')
    p.add_argument('-t','--trim',action='store_true',help='trim sequence before and within barcode')
    p.add_argument('-m','--max_score',dest='max_score',type=int,default=2, help="max(mismatch+indel) allowed for barcode match, otherwise move reads into 'unassigned' file. default: 2") 
    p.add_argument('-p','--parallel',type=int,default=5,help="number of thread")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()



def fuzzy_substring(needle, haystack):
    """Calculates the fuzzy match of needle in haystack,
    using a modified version of the Levenshtein distance
    algorithm.
    The function is modified from the levenshtein function
    in the bktree module by Adam Hupp
    http://ginstrom.com/scribbles/2007/12/01/fuzzy-substring-
    matching-with-levenshtein-distance-in-python/"""
    m, n = len(needle), len(haystack)

    # base cases
    if m == 1:
        return not needle in haystack
    if not n:
        return m

    row1 = [0] * (n+1)
    minS=m
    for i in range(0,m):
        row2 = [i+1]
        for j in range(0,n):
            cost = ( needle[i] != haystack[j] )

            row2.append(   min(row1[j+1]+1, # deletion
                               row2[j]+1, #insertion
                               row1[j]+cost) #substitution
                          )
            if i == m-1:
                if row2[j+1] <= minS:
                    minS=row2[j+1]
                    end=j+1
        row1 = row2
    return minS, end
'''
TEST:
print (fuzzy_substring('ACTC', 'C_ ATCG'))
print (fuzzy_substring('ACTC', 'C_ ACTGG'))
print (fuzzy_substring("ACTAAC", "ACTAACTAGCCATGCAATGGCTAG"))
'''

def Assign_barcode(record,files):
    miScore=barcode_len
    end = 0
    title1 = record[0]
    seq1 = record[1]
    qual1 = record[2]
    title2 = record[3]
    seq2 = record[4]
    qual2 = record[5]
    if Range == 0:
        read_seq=seq1
    elif len(Range) == 2:
        read_seq=seq1[(int(Range[0])-1):int(Range[1])]
    for i in barcodes:
        score,j=fuzzy_substring(i,read_seq)
        if score<miScore:
            barcode=i
            end=j
            miScore=score
#    record2=records2.next()
    if miScore>max_score:
        files.Output_Type(title1, seq1, qual1, title2, seq2, qual2, 'unassign')
  #      SeqIO.write(record,output1['unassign'],type)
  #      SeqIO.write(record2,output2['unassign'],type)
    else:
        files.Output_Type(title1, seq1, qual1, title2, seq2, qual2, barcode, end)
#          SeqIO.write(record[end:],output1[barcode],type)
#          SeqIO.write(record,output1[barcode],type)
#      SeqIO.write(record2,output2[barcode],type)


def do_stuff(io):
    while True:
        #if q.empty():
        #    continue
        try:
            item=io.Read_Type()
            Assign_barcode(item,io)
#       do_work(obj,item)
        except:
            # end of file
            break

def push_stuff(q, num_thread, rec1, rec2):
    try:
        while True:
            record=rec1.next()
            record2=rec2.next()
            q.put((record,record2))
    except:
        for i in range(num_thread):
            q.put("kill")

class MainIO():
    def __init__(self, record1, record2, name1, name2):
        self.rec1 = record1
        self.rec2 = record2
        self.output1={}
        self.output2={}
        self.writeLocker={}
        self.readLocker = threading.Lock()
        for b in barcodes:
            self.output1[b]=open(b+name1,'w')
            self.output2[b]=open(b+name2,'w')
            self.writeLocker[b]=threading.Lock()
            #self.locker2[b]=threading.Lock()

    def Read_Type(self):
        self.readLocker.acquire()
        try:
            (title1, seq1, qual1) = self.rec1.next()
            (title2, seq2, qual2) = self.rec2.next()
        finally:
            self.readLocker.release()
        return (title1, seq1, qual1, title2, seq2, qual2)

    def Output_Type(self,title1,seq1,qual1,title2,seq2,qual2,barcode,end=0):
        self.writeLocker[barcode].acquire()
        if end > 0 and Trim:
            self.output1[barcode].write(("@%s" + os.linesep + "%s" + os.linesep 
                + "+" + os.linesep + "%s" + os.linesep) 
                % (title1, seq1[end:], qual1[end:]))
        else:
            self.output1[barcode].write(("@%s" + os.linesep + "%s" + os.linesep 
                + "+" + os.linesep + "%s" + os.linesep)
                % (title1, seq1, qual1))
        self.output2[barcode].write(("@%s" + os.linesep + "%s" + os.linesep 
            + "+" + os.linesep + "%s" + os.linesep)
            % (title2, seq2, qual2))
        self.writeLocker[barcode].release()

    
def Main():
    args=ParseArg()

#    enable_thread_profiling()
#    yappi.set_clock_type("wall")
#    yappi.start()

#    global type
#    type = "fastq"
#    if args.fastq:
#       type="fastq"
#    elif args.fasta:
#       type="fasta"

    global Range
    Range=args.range

    global Trim
    Trim=args.trim

    global max_score
    max_score=args.max_score

    name1=args.input1.split('/')[-1]
    name2=args.input2.split('/')[-1]
    #----------- read barcode ----------
    global barcodes
    global barcode_len
    barcodes=[]
    for i in open(args.barcode,'r'):
        i=i.strip()
        barcodes.append(i)
        barcode_len=len(i)

    barcodes.append('unassign')


    #-----------------------------------

    records=QualityIO.FastqGeneralIterator(open(args.input1,"rU"))
    records2=QualityIO.FastqGeneralIterator(open(args.input2,"rU"))
    
    Files=MainIO(records, records2, name1, name2)

    print "start to assign sequence to different barcodes..."    
    print "----------"
    num_thread=args.parallel
#    q=Queue(maxsize=10000)

#    feeder = threading.Thread(target = push_stuff, args = (q, num_thread, records, records2))
#    feeder.start()

    workers = []
    for i in range(num_thread):
        worker=threading.Thread(target=do_stuff,args=(Files,))
        workers.append(worker)
        #worker.setDaemon(True)
        worker.start()

#    print >>sys.stderr,"Finish reading records"


#    q.join()

#    feeder.join()
    for i in range(num_thread):
        workers[i].join()

    #get_thread_stats().print_stats()
#    yappi.get_func_stats().print_all()
#    yappi.get_thread_stats().print_all()

    # for barcode in barcodes:
    #     output1[barcode].close()
    #     output2[barcode].close()

        
if __name__=="__main__":
    Main()

