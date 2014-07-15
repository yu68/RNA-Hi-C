
import sys,os,argparse
import string
import json
import numpy as np
from RNAstructure import RNAstructure
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()


def ParseArg():
    p=argparse.ArgumentParser(description="plot RNA structure with distribution of digested end, refine structure with loc of digested end")
    p.add_argument("ID",type=str,help="Ensembl gene ID of RNA")
    p.add_argument("linkedPair",type=str,help="file for information of linked pairs, which is output of 'Stitch-seq_Aligner.py'")
    p.add_argument('-g',"--genomeFa",type=str,default="/home/yu68/Software/bowtie-0.12.7/indexes/mm9.fa",help="genomic sequence,need to be fadix-ed")
    p.add_argument("-R","--RNAstructureExe",type=str,default="/home/yu68/Software/RNAstructure/exe/",help="folder of RNAstrucutre suite excutable")
    p.add_argument("-a","--acceptDot",type=str,help="accepted structure in dot format, for comparing of accuracy, no comparison if not set")
    p.add_argument("-o","--output",help="output distribution of digested sites with dot structures, can be format of eps, pdf, png,...")
    p.add_argument('-s','--samtool_path',dest='spath', type=str,metavar='samtool_path',help="path for the samtool program",default='samtools')
    p.add_argument("-v",'--varna',type=str,help="path for the VARNA visualization for RNA", default="/home/yu68/Software/VARNA/VARNAv3-9.jar")
    p.add_argument("-c","--colorMapStyle",type=str,help='style of color map, choose from: "red", "blue", "green", "heat", "energy", and "bw",default:"heat"',default="heat")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

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

def relativeLoc(geneStart,geneEnd,geneStrand,loc):
    ''' find relative location of 'loc' on the gene, 0 indexed  '''
    if geneStrand=="+":
        rela_loc = loc-geneStart
    else:
        rela_loc = geneEnd-loc
    return rela_loc

def Count_DS_SS(geneSeq,geneDot,count,ax=None,title=""):
    ''' count of digested sites in non base pairing and base pairing nucleotide with box plot '''
    count_SS=[]  # store single strand (non base pairing) count
    count_DS=[]  # store double strand (base pairing) count
    for i in range(len(geneSeq)):
        if geneDot[i]==".":
            count_SS.append(count[i])
        else:
            count_DS.append(count[i])
    if ax!=None:
        ax.boxplot([count_SS,count_DS],0,"",patch_artist=True,widths=0.5)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_title(title)
        ax.set_ylabel("Digested site count")
        ax.set_xticklabels(["SS","DS"])
    print >> sys.stderr, np.mean(count_SS), np.mean(count_DS)
    return count_SS,count_DS

def DrawShadowDot(geneDot,ax,ymin,ymax,col="g"):
    for i in range(len(geneDot)):
        if geneDot[i]==".":
            rect = matplotlib.patches.Rectangle((i+0.5,ymin),1,ymax-ymin,color=col,alpha=0.3)
            ax.add_patch(rect)

def generateJSON(geneDot,count,name,out_file_name):
    '''
    generate JSON file for uploading into RNA 2D browser for visualization
    '''
    RNA = {}
    count=[int(x) for x in count]
    RNA["ideograms"]=[{"id":name,"length":len(geneDot),"color":"grey"}]
    RNA["plottracks"]=[{"name":"coverage","color":"blue","values":[{"chr":name,"values":count}]}]
    RNA["structs"]=[{"chr":name,"struct":geneDot}]
    out = open(out_file_name,'w')
    print >>out, json.dumps(RNA, indent=2)


def Main():
    GeneLoc='/home/yu68/bharat-interaction/new_lincRNA_data/all_RNAs-rRNA_repeat.txt'
    args=ParseArg()
    ID = args.ID
    gene=os.popen('grep "'+ID+'" '+GeneLoc).read().split('\n')[0].split('\t')
    geneStart=int(gene[1])
    geneEnd=int(gene[2])
    geneChro=gene[0]
    geneName=gene[4]
    geneStrand=gene[5]
    if not os.path.isfile(args.linkedPair):
        print "LinkedPair file is not exist, please check!!"
        sys.exit(0)
    if geneStrand=='+':
        Pairs=os.popen('grep "'+ID+'" '+args.linkedPair+" | awk '$7==$16' | sort -k3,3n -k11,11n").read().split('\n')
    else:
        Pairs=os.popen('grep "'+ID+'" '+args.linkedPair+" | awk '$7==$16' | sort -k2,2n -k12,12n").read().split('\n')

    geneSeq = fetchSeq(geneChro,geneStart,geneEnd,geneStrand,args.genomeFa,args.spath)
    print >> sys.stderr, geneName, len(geneSeq)
    print >> sys.stderr, "\nsequence: "
    print >> sys.stderr, "\t", geneSeq

    # predict strucutre use RNAstructure
    RNA_prog = RNAstructure(exe_path=args.RNAstructureExe)
    seq, geneDot = RNA_prog.Fold(geneSeq,ct_name="temp_predict.ct")
    assert(seq==geneSeq)
    
    # count for digested point in each position
    count = np.zeros(len(geneSeq))
    for p in Pairs:
        if p.strip()=="": continue
        p=p.split('\t')
        assert(geneChro==p[0]==p[9])
        p1_start=int(p[1])
        p1_end=int(p[2])
        p2_start=int(p[10])
        p2_end=int(p[11])
        if p[3]=="+" and p[12]=="-":
            try:
                count[relativeLoc(geneStart,geneEnd,geneStrand,p1_start)]+=1
                count[relativeLoc(geneStart,geneEnd,geneStrand,p1_end)]+=1
            except: pass
            try:
                count[relativeLoc(geneStart,geneEnd,geneStrand,p2_end)]+=1
                count[relativeLoc(geneStart,geneEnd,geneStrand,p2_start)]+=1
            except: pass
        elif p[3]=="-" and p[12]=="+":
            try:
                count[relativeLoc(geneStart,geneEnd,geneStrand,p1_end)]+=1
                count[relativeLoc(geneStart,geneEnd,geneStrand,p1_start)]+=1
            except: pass
            try:
                count[relativeLoc(geneStart,geneEnd,geneStrand,p2_start)]+=1
                count[relativeLoc(geneStart,geneEnd,geneStrand,p1_end)]+=1
            except: pass
    print >> sys.stderr, ",".join(["%.2f"%(math.log(f+1)) for f in count])

    assert(len(geneSeq)==len(geneDot))
    

    sso=open("temp_sso.txt",'w')
    for i in range(len(geneSeq)):
        if count[i]!=0:
            print >>sso, "%d\t%.5f"%(i+1,-math.log(count[i]+1)/2)
    sso.close()

    # refine structure using single strand constraint
    seq, geneDot_refine = RNA_prog.Fold(geneSeq,ct_name="temp_refine.ct",sso_file="temp_sso.txt")
    assert(seq==geneSeq)


    # compare with accepted structure
    if args.acceptDot!=None:
        os.system("dos2unix %s"%(args.acceptDot))
        os.system("%s/dot2ct %s temp_accept.ct"%(args.RNAstructureExe,args.acceptDot))
        sen_predict, PPV_predict, Number_predict = RNA_prog.scorer("temp_predict.ct","temp_accept.ct")
        sen_refine, PPV_refine, Number_refine = RNA_prog.scorer("temp_refine.ct","temp_accept.ct")
        print >> sys.stderr, "Prediction: %.2f, %.2f, %d"%(sen_predict, PPV_predict, Number_predict)
        print >> sys.stderr, "Refine: %.2f, %.2f, %d"%(sen_refine, PPV_refine, Number_refine)
        # predicted structure again
        seq, geneDot = RNA_prog.Fold(geneSeq,ct_name="temp_predict.ct",Num=Number_predict)
        assert(seq==geneSeq)
        # refine structure using single strand constraint again
        seq, geneDot_refine = RNA_prog.Fold(geneSeq,ct_name="temp_refine.ct",sso_file="temp_sso.txt",Num=Number_refine)
        assert(seq==geneSeq)

    print >> sys.stderr, "structure: "
    print >> sys.stderr, "\t", geneDot

   
 
    print >> sys.stderr, "refined structure: "
    print >> sys.stderr, "\t", geneDot_refine

    # output JSON file for visualization
    generateJSON(geneDot,count,geneName,"structure_%s-%s.json"%(ID,geneName))
    generateJSON(geneDot_refine,count,geneName,"structure_%s-%s_refine.json"%(ID,geneName))



    # plot digested site distribution
    fig = plt.figure(figsize=(10,4))
    ax1 = plt.subplot2grid((1,8),(0,0),colspan=5,frameon=False)  # plot the digested site distribution
    ax2 = plt.subplot2grid((1,8),(0,5))  # predicted counts of digested site in predicted structure
    ax3 = plt.subplot2grid((1,8),(0,6))  # predicted counts of digested site in refined structure
    ax4 = plt.subplot2grid((1,8),(0,7))  # predicted counts of digested site in accepted structure


    # count digested site in ds ss for each structure
    print >> sys.stderr, "Predicted Structure count:"
    Count_DS_SS(geneSeq,geneDot,count,ax2,"Predicted")

    print >> sys.stderr, "Refined structure count:"
    Count_DS_SS(geneSeq,geneDot_refine,count,ax3,"Refined")
    

    ax1.bar(range(len(count)),[math.log(f+1) for f in count],edgecolor = "none",facecolor="#545454")
    ax1.set_ylabel("log(count+1)")
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.set_xlim((1,len(geneSeq)))
    DrawShadowDot(geneDot,ax1,1.5,2.0,col='g')
    DrawShadowDot(geneDot_refine,ax1,2.5,3.0,col='y')
    if args.acceptDot!=None:
        ax1.text(0.01,0.90,"Prediction: SEN- %.3f; PPV- %.3f"%(sen_predict, PPV_predict),horizontalalignment="left",verticalalignment="bottom",transform=ax1.transAxes,color='g',fontsize=10)
        ax1.text(0.01,0.99,"Refined: SEN- %.3f; PPV- %.3f"%(sen_refine, PPV_refine),horizontalalignment="left",verticalalignment="bottom",transform=ax1.transAxes,color='y',fontsize=10)



    # draw structure with color map
    cmd = ["java -cp"]
    cmd.append(args.varna)
    cmd.append("fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN")
    cmd.append('"%s"'%(geneSeq))
    cmd.append("-structureDBN")
    cmd.append('"%s"'%(geneDot))
    cmd.append("-colorMap")
    cmd.append('"%s"'%(";".join([str(math.log(f+1,2)) for f in count])))
    cmd.append("-colorMapStyle %s"%(args.colorMapStyle))
    cmd.append("-o structure_%s-%s_cm.eps"%(ID,geneName))
    cmd = " ".join(cmd)
    os.system(cmd)

    # draw refined structure with color map
    cmd = ["java -cp"]
    cmd.append(args.varna)
    cmd.append("fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN")
    cmd.append('"%s"'%(geneSeq))
    cmd.append("-structureDBN")
    cmd.append('"%s"'%(geneDot_refine))
    cmd.append("-colorMap")
    cmd.append('"%s"'%(";".join([str(math.log(f+1,2)) for f in count])))
    cmd.append("-colorMapStyle %s"%(args.colorMapStyle))
    cmd.append("-o structure_%s-%s_cm_refine.eps"%(ID,geneName))
    cmd = " ".join(cmd)
    os.system(cmd)

    if args.acceptDot==None:
        plt.savefig(args.output)
        sys.exit(0)
    # draw accept structure with color map
    accept = open(args.acceptDot).read().split('\n')
    geneSeq = accept[1]
    geneDot_accept = accept[2]

    print >> sys.stderr, "Accept structure count: "
    Count_DS_SS(geneSeq,geneDot_accept,count,ax4,"Accepted")

    ax1.set_xticks(range(len(count)))
    ax1.set_xticklabels(list(geneDot_accept))
    plt.savefig(args.output)

    cmd = ["java -cp"]
    cmd.append(args.varna)
    cmd.append("fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN")
    cmd.append('"%s"'%(geneSeq))
    cmd.append("-structureDBN")
    cmd.append('"%s"'%(geneDot_accept))
    cmd.append("-colorMap")
    cmd.append('"%s"'%(";".join([str(math.log(f+1,2)) for f in count])))
    cmd.append("-colorMapStyle %s"%(args.colorMapStyle))
    cmd.append("-o structure_%s-%s_cm_accept.eps"%(ID,geneName))
    cmd = " ".join(cmd)
    os.system(cmd)

    os.system("rm temp*")

if __name__=="__main__":
    Main()
