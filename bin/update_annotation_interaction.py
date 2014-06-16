import sys,os
from xplib import DBI
from Annotation import annotation

if len(sys.argv)==1 or sys.argv[1] in ["-h","--help"]:
    print >> sys.stderr, "\nUsage:"
    print >> sys.stderr, "  python %s {interaction_file} > updated_interaction_file \n"%(os.path.basename(__file__))
    sys.exit(0)

class Bed:
    def __init__(self,x,**kwargs):
        self.chr=x[0]
        self.start=int(x[1])
        self.stop=int(x[2])
        if len(x)>=6:
            self.type=x[3]
            self.name=x[4]
            self.subtype=x[5]
        self.center=int((self.start+self.stop)/2)
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def overlap(self,bed2,n):
        return (self.chr==bed2.chr) and (self.start+n<bed2.stop) and (bed2.start+n<self.stop)
    def str_region(self):
        return self.chr+":%d-%d"%(self.start,self.stop)


def read_interaction(File,s):
    '''
    s: start column number for second part of interaction
    '''
    a=open(File,'r')
    for l in a.read().split('\n'):
        if l.strip()=="": continue
        lsep=l.split('\t')
        if lsep[3] in ['+','-']:
            bed1=Bed(lsep[0:3],strand=lsep[3])
            bed2=Bed(lsep[s:(s+3)],strand=lsep[s+3])
        else:
            bed1=Bed(lsep[0:3])
            bed2=Bed(lsep[s:(s+3)])
        yield (bed1,bed2,lsep)

# annotation files
db="/home/yu68/bharat-interaction/new_lincRNA_data/all_RNAs-rRNA_repeat.txt"
db_detail="/home/yu68/bharat-interaction/new_lincRNA_data/Ensembl_mm9.genebed"
db_repeat="/home/yu68/bharat-interaction/new_lincRNA_data/mouse.repeat.txt"
print >>sys.stderr, "Indexing annotation files..."
ref_allRNA=DBI.init(db,"bed") # the DBI init file for bed6 file of all kinds of RNA
ref_detail=DBI.init(db_detail,"bed") # the DBI init file for bed12 file of lincRNA and mRNA with intron, exon, UTR
ref_repeat=DBI.init(db_repeat,"bed")

print >>sys.stderr, "Start to update..."
for l in read_interaction(sys.argv[1],7):
    l[2][3:6] = annotation(l[0],ref_allRNA,ref_detail,ref_repeat)
    l[2][10:13] = annotation(l[1],ref_allRNA,ref_detail,ref_repeat)
    print "\t".join(l[2])

    
