import argparse
import sys
import string,random
from pysam import Fastafile
from RNAstructure import RNAstructure
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from altschulEriksonDinuclShuffle import dinuclShuffle

if len(sys.argv)<2:
    print "Usage:"
    print "  python freeEnergy.py fragment_file.txt figure.pdf"
    print "\nOutput:"
    print "  txt file with two groups of free energies"
    print "  pdf file for the plot of distribution of free energies, the name the pdf file could be assigned in the command line"
    sys.exit(0)

mappedPair=open(sys.argv[1],'r')
output_fig = sys.argv[2]
output_txt = sys.argv[2].split(".")[0]+".txt"


rev_table=string.maketrans('ACGTacgt', 'TGCAtgca')
def revcomp(seq, rev_table):
    return seq.translate(rev_table)[::-1]

def shuffleseq(seq1):
    return "".join([random.choice("ACGT") for i in range(len(seq1))])
    #temp=list(seq1)
    #random.shuffle(temp)
    #return "".join(temp)

RNA_prog = RNAstructure(exe_path="/home/yu68/Software/RNAstructure/exe/")
energies=[]
random_energies=[]

output = open(output_txt,'w')

print >> output, "Paire_ID\tType_1\tName_1\tType_2\tName_2\treal_FE\trandom_FE\tLength"

i=0
for l in mappedPair.read().split('\n'):
    if l.strip()=="": continue
    l=l.split('\t')
    seq1=l[4]
    seq2=l[13]
    strand1=l[3]
    strand2=l[12]
    if strand1=='-':
        seq1=revcomp(seq1,rev_table)
    if strand2=='+':
        seq2=revcomp(seq2,rev_table)
    seq1=seq1.replace("N","T")
    seq2=seq2.replace("N","T")
    energy=RNA_prog.DuplexFold(seq1,seq2)
    printline = l[8]+"\t"+l[5]+"\t"+l[6]+"\t"+l[14]+"\t"+l[15]+"\t%.4f"%(energy)
    energies.append(energy)
    try:
        r_seq1=shuffleseq(seq1.replace("T","U"))
        r_seq2=shuffleseq(seq2.replace("T","U"))
    except:
        print seq1,seq2
    #r_seq1=shuffleseq(seq1)
    #r_seq2=shuffleseq(seq2)
    energy=RNA_prog.DuplexFold(r_seq1,r_seq2)
    printline += '\t%.4f\t%.1f'%(energy,(len(seq1)+len(seq2))/2.0)
    random_energies.append(energy)
    i+=1
    print >>output, printline
    if i%100==0:
        print >>sys.stderr, "processing %d %.2f %.2f\r"%(i,sum(energies)/float(i),sum(random_energies)/float(i)),
    if i==60000:
        break

energies=np.array(energies)
random_energies=np.array(random_energies)
print np.mean(energies),np.mean(random_energies),"         "
print np.std(energies),np.std(random_energies),"         "

density = stats.gaussian_kde(energies)
density.covariance_factor = lambda : .10
density._compute_covariance()
density_r = stats.gaussian_kde(random_energies)
density_r.covariance_factor = lambda : .10
density_r._compute_covariance()
xs = np.linspace(min(energies),0,500)
plt.figure(figsize=(6,4))
plt.plot(xs,density(xs),label="ChimericRNA")
plt.plot(xs,density_r(xs),label="Random")
plt.xlim(-60,0)
_,pvalue1=stats.wilcoxon(energies,random_energies)
_,pvalue2=stats.ttest_rel(energies,random_energies)
plt.legend()
print "wilcoxon signed-rank test:",pvalue1,"Paired t-test:",pvalue2
plt.savefig(output_fig)
plt.close()    
