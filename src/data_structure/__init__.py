

class annotated_bed():
    """
    To store, compare, cluster for the genomic regions with RNA annotation information. Utlized in the program :ref:`Select_stronginteraction_pp.py<Step6>`
    """
    
    def __init__(self,x=None,**kwargs):
        """
        Initiation example:
        
        >>> str="chr13\t40975747\t40975770\tATTAAG...TGA\tprotein_coding\tgcnt2\tintron"
        >>> a=annotated_bed(str)
        or
        >>> a=annotated_bed(chr="chr13",start=40975747,end=40975770,type="protein_coding",)
        """
        if x is not None:
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
    def overlap(self,other): # for the purpose of finding overlap between regions
        """
        Find overlap between regions

        :param other: another :class:`annotated_bed` object
        :returns: boolean

        """
        return ((self.chr==other.chr)&(self.end>other.start)&(self.start<other.end))

    def __lt__(self, other): # for the purpose of ordering clusters
        """
        Compare two objects self and other when they are not **overlapped**
        
        :param other: another :class:`annotated_bed` object
        :returns: boolean -- "None" if overlapped.

        Example:
        
        >>> a=annotated_bed(chr="chr13",start=40975747,end=40975770)
        >>> b=annotated_bed(chr="chr13",start=10003212,end=10005400)
        >>> print a>b
        False

        """
        if self.overlap(other): return "None"
        if self.chr==other.chr:
            return ((self.start<other.start)&(self.end<other.end))
        else:
            return (self.chr<other.chr)

    def Cluster(self,c):
        """
        Store cluster information of self object
        
        :param c: cluster index
        
        Example:
        
        >>> a=annotated_bed(chr="chr13",start=40975747,end=40975770)
        >>> a.Cluster(3)
        >>> print a.cluster
        3

        .. Note::
        
           a.cluster will be the count information when a become a cluster object in :ref:`Select_stronginteraction_pp.py<Step6>`

        """
        self.cluster=c
        #self.c_info=c_info
    def Update(self,S,E): # used for update or expand cluster locations
        """
        Update the upper and lower bound of the cluster after adding segments using Union-Find.
        
        :param S: start loc of the newly added genomic segment
        :param E: end loc of the newly added genomic segment
        
        Example:
        
        >>> a=annotated_bed(chr="chr13",start=40975747,end=40975770)
        >>> a.Update(40975700,40975800)
        >>> print a.start, a.end
        40975700 40975800
        """

        self.start=min(self.start,S)
        self.end=max(self.end,E)
    def __str__(self):
        """
        Use print function to output the cluster information (chr, start, end, type, name, subtype,cluster)
        
        Example:
        
        >>> str="chr13\t40975747\t40975770\tATTAAG...TGA\tprotein_coding\tgcnt2\tintron"
        >>> a=annotated_bed(str)
        >>> a.Cluster(3)
        >>> a.Update(40975700,40975800)
        >>> print a
        "chr13\t40975700\t40975800\tprotein_coding\tgcnt2\tintron\t3"
        """
        return "\t".join(str(f) for f in [self.chr,self.start,self.end,self.type,self.name,self.subtype,self.cluster])
    # self.cluster become number of regions within cluster for cluster_pool object
