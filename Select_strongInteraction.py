import sys, argparse
from time import time
from operator import attrgetter
from scipy.stats import hypergeom

class annotated_bed():
    def __init__(self,x=None,**kwargs):
        if type(x)==type("str"):
            x=x.split("\t")
        self.chr=str(x[0]).strip()
        self.start=int(x[1])
        self.end=int(x[2])
        try:
            self.seq=str(x[3]).strip()
            self.type=str(x[4]).strip()
            self.name=str(x[5]).strip()
            self.subtype=str(x[6]).strip()
        except:
            pass
        for key in kwargs.keys():
            setattr(self,key,kwargs[key])
    def __lt__(self, other): # for the purpose of ordering clusters
        if self.chr==other.chr:
            return ((self.start<other.start)&(self.end<other.end))
        else:
            return (self.chr<other.chr)
    def overlap(self,other): # for the purpose of finding overlap between regions
        return ((self.chr==other.chr)&(self.end>other.start)&(self.start<other.end))
    def Cluster(self,c,c_info):
        self.cluster=c
        self.c_info=c_info
    def Update(self,S,E): # used for update or expand cluster locations
        self.start=min(self.start,S)
        self.end=max(self.end,E)
    def __str__(self):
        return "\t".join(str(f) for f in [self.chr,self.start,self.end,self.type,self.name,self.subtype,self.cluster])
    # self.cluster become number of regions within cluster for cluster_pool object


################################################
## implement of weighted quick union ##########

def get_root(lis,item):
    """
    return the root of an item
    """
    root=item
    while lis[root] != root:
        root = lis[root]
    return root    
def connected(lis,item1,item2):
    """
    Checks if two items are connected
    or not.
    Connected means both items have the
    same root
    """
    root1=get_root(lis,item1)
    root2= get_root(lis,item2)
    return True if root1==root2 else False
 
def quick_union(lis,tree_sz,item1,item2):
    """
    Join two nodes(join the root of smaller tree to root of bigger tree)
    """
    root1=get_root(lis,item1)
    root2= get_root(lis,item2)
    
    if root1 != root2:
        if tree_sz[root1]==tree_sz[root2]:       #if size of both tree is same
            lis[root2]=root1
            tree_sz[root1]+=tree_sz[root2]
            tree_sz[root2]=0                      # set the size of tree was merged to 0
 
        elif tree_sz[root1] < tree_sz[root2]:   
            lis[root1]=root2
            tree_sz[root2]+=tree_sz[root1]
            tree_sz[root1]=0
 
        elif tree_sz[root1] > tree_sz[root2]:
            lis[root2]=root1
            tree_sz[root1]+=tree_sz[root2]
            tree_sz[root2]=0    
    #else:
       # print "Already connected"
    return [lis,tree_sz]

###########################################################


##########################################################
#cluster the regions together if they are connected (overlapped) using Union-find approach
def cluster_regions(part,min_clusterS):
    N=len(part)
    List=range(N)
    tree_sz=[1]*N
    cluster_loc=part
    for i in range(N):
        for j in range(i+1,N):
            if part[i].overlap(part[j]):
                [List,tree_sz]=quick_union(List,tree_sz,i,j)
                #update cluster location (expand)
                cluster_loc[get_root(List,j)].Update(min(part[i].start,part[j].start),max(part[i].end,part[j].end))
            if part[i]<part[j]:
                break
    
    c_pool=[] # cluster pool
    for i in range(N):
        c=get_root(List,i)
        c_info="%d:%d"%(cluster_loc[c].start,cluster_loc[c].end)
        part[i].Cluster(get_root(List,i),c_info)
        c_pool.append(c)

    cluster_pool={}
    for c in set(c_pool):
        if c_pool.count(c)>=min_clusterS:
            cluster_pool[c]=cluster_loc[c]
            cluster_pool[c].cluster=c_pool.count(c)
    return (cluster_pool)


#parameters

def ParseArg():
    p=argparse.ArgumentParser(description="find strong interactions from paired genomic location data",epilog="need Scipy for hypergeometric distribution")
    p.add_argument("-i","--input",type=str,required=True,help="input file which is the output file of Stitch-seq-Aligner.py")
    p.add_argument("-M","--min_clusterS",type=int,default=5,help="minimum number of segments allowed in each cluster")
    p.add_argument("-m","--min_interaction",type=int,default=3,help="minimum number of interactions to support a strong interaction")
    p.add_argument('-p',"--p_value",type=float,default=0.05,help="the p-value based on hypergeometric distribution to call strong interactions")
    p.add_argument('-o','--output',type=str,help="specify output file")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()


def Main():
    t1=time()

    args=ParseArg()
    inp = open(args.input, 'r')
    min_clusterS=args.min_clusterS
    min_interaction=args.min_interaction
    p_value=args.p_value
    output=open(args.output,'w')

    #store genomic location of part1 and part2
    part1=[]
    part2=[]


    k=0
    for line in inp.read().split('\n'):
        if line=='': continue
        line=line.strip().split('\t')
        k=k+1
        part1.append(annotated_bed(line[0:7],id=k))
        part2.append(annotated_bed(line[8:],id=k))
    
    if len(part1)!=len(part2):
        print >> sys.stderr, "number of regions in two part not match!!"
        sys.exit(0)

    # sort in genomic order, easy for clustering
    part1=sorted(part1, key=attrgetter('start'))
    part1=sorted(part1, key=attrgetter('chr'))
    part2=sorted(part2, key=attrgetter('start'))
    part2=sorted(part2, key=attrgetter('chr'))


    print >>sys.stderr,"# Generating clusters for two parts...\n"
    cluster_pool1=cluster_regions(part1,min_clusterS)
    cluster_pool2=cluster_regions(part2,min_clusterS)

    # sort back to pair two parts
    part1=sorted(part1, key=attrgetter('id'))
    part2=sorted(part2, key=attrgetter('id'))

    print >>sys.stderr,"# cluster number for part1 is %d"%(len(cluster_pool1))
    print >>sys.stderr,"# cluster number for part1 is %d"%(len(cluster_pool2))

    c_interaction=[]
    for i in range(len(part1)):
        region1=str(part1[i])
        region2=str(part2[i])
        print region1+'\t'+region2
        c_interaction.append("%d--%d"%(part1[i].cluster,part2[i].cluster))

    print >> sys.stderr,"# finding strong interactions from clusters..."
    for i in cluster_pool1:
        for j in cluster_pool2:
            interaction="%d--%d"%(i,j)
            count=c_interaction.count(interaction)
            if count<min_interaction: continue
            count1=cluster_pool1[i].cluster
            count2=cluster_pool2[j].cluster
            real_p=1-hypergeom.cdf(count,len(part1),count1,count2)
            if real_p<=p_value:
                print >> output,str(cluster_pool1[i])+'\t'+str(cluster_pool2[j])+'\t%d\t%.5f'%(count,real_p)
    

    print >> sys.stderr,"# Cost time: %.2f s"%(time()-t1)


if __name__=="__main__":
    Main()
