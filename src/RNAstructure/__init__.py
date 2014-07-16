# -------------
# Pengfei Yu
#

"""Interface class for `RNAstructure <http://rna.urmc.rochester.edu/RNAstructure.html>`_ executable programs. """


import os
import re
import subprocess
import tempfile
import sys

class RNAstructure(object):
    """
    Interface class for `RNAstructure <http://rna.urmc.rochester.edu/RNAstructure.html>`_ executable programs.
    """
    def __init__(self, exe_path=None):
        """
        Initiation of object
        
        :param exe_path: the folder path of the RNAstructure executables
        
        Example:

        >>> from RNAstructure import RNAstructure
        >>> RNA_prog = RNAstructure(exe_path="/home/yu68/Software/RNAstructure/exe/")
        """
        if exe_path==None:
            self.exe_path=''
        else:
            self.exe_path = exe_path
        self.datapath=self.exe_path+"/../data_tables/"
        os.environ['DATAPATH'] = self.datapath
    
    def DuplexFold(self,seq1=None,seq2=None,dna=False):
        '''
        Use "DuplexFold" program to calculate the minimum folding between two input sequences

        :param seq1,seq2: two DNA/RNA sequences as string, or existing fasta file name
        :param dna: boolean input. Specify then DNA parameters are to be used
        :returns: minimum binding energy, (unit: kCal/Mol)

        Example:

        >>> from RNAstructure import RNAstructure
        >>> RNA_prog = RNAstructure(exe_path="/home/yu68/Software/RNAstructure/exe/")
        >>> seq1 = "TAGACTGATCAGTAAGTCGGTA"
        >>> seq2 = "GACTAGCTTAGGTAGGATAGTCAGTA"
        >>> energy=RNA_prog.DuplexFold(seq1,seq2)
        >>> print energy
        '''
        cmd = [os.path.join(self.exe_path,"DuplexFold")]
        #sequences
        seq1_file=None
        if seq1.endswith("fasta"):
            cmd.append(seq1)
        elif seq1 is not None:
            seq1_file = tempfile.NamedTemporaryFile(mode='w')
            seq1_file.write(">seq1\n")
            seq1_file.write(seq1)
            seq1_file.flush()
            cmd.append(seq1_file.name)
        seq2_file=None
        if seq2.endswith("fasta"):
            cmd.append(seq2)
        elif seq2 is not None:
            seq2_file = tempfile.NamedTemporaryFile(mode='w')
            seq2_file.write(">seq2\n")
            seq2_file.write(seq2)
            seq2_file.flush()
            cmd.append(seq2_file.name)
        # output ct file
        ct_file = tempfile.NamedTemporaryFile(mode='w')
        cmd.append(ct_file.name)
        # DNA or RNA
        if dna==True:
            cmd.append('-d')
        # Excuting program
        cmd = " ".join(cmd)
        os.system(cmd+"> /dev/null 2>/dev/null")
        if seq1_file is not None:
            seq1_file.close()
        if seq2_file is not None:
            seq2_file.close()
        #parsing result
        out = open(ct_file.name).read()
        decoded = re.search(r'ENERGY = (?P<prob>\S+)', out)
        ct_file.close()
        try:
            return float(decoded.groupdict()['prob'])
        except:
            return 0.0

    def Fold(self,seq=None,ct_name=None,sso_file=None,Num=1):
        '''
        Use "Fold" program to predict the secondary structure and output dot format.

        :param seq: one DNA/RNA sequence as string, or existing fasta file name
        :param ct_name: specify to output a ct file with this name, otherwise store in temp, default: None
        :param sso_file: give a single strand offset file, format see http://rna.urmc.rochester.edu/Text/File_Formats.html#Offset
        :param Num: choose Num th predicted structure
        :returns: dot format of RNA secondary structure and RNA sequence.

        Example:

        >>> from RNAstructure import RNAstructure
        >>> RNA_prog = RNAstructure(exe_path="/home/yu68/Software/RNAstructure/exe/")
        >>> seq = "AUAUAAUUAAAAAAUGCAACUACAAGUUCCGUGUUUCUGACUGUUAGUUAUUGAGUUAUU"
        >>> sequence,dot=RNA_prog.Fold(seq)
        >>> assert(seq==sequence)
        '''
        cmd = [os.path.join(self.exe_path,"Fold")]
        cmd2 = [os.path.join(self.exe_path,"ct2dot")]
        #sequences
        seq_file=None
        if seq.endswith("fasta"):
            cmd.append(seq)
        elif seq is not None:
            seq_file = tempfile.NamedTemporaryFile(mode='w')
            seq_file.write(">seq\n")
            seq_file.write(seq)
            seq_file.flush()
            cmd.append(seq_file.name)
        # output ct file
        if ct_name==None:
            ct_file = tempfile.NamedTemporaryFile(mode='w')
            cmd.append(ct_file.name)
            cmd2.append(ct_file.name)
        else:
            cmd.append(ct_name)
            cmd2.append(ct_name)
        # output dot file
        dot_file = tempfile.NamedTemporaryFile(mode='w')
        #cmd.append('-mfe') # only the best one
        # single strand offset
        if sso_file!=None:
            cmd.append("-sso")
            cmd.append(sso_file)     

        cmd2.append('%d'%(Num))  # Num th structure
        cmd2.append(dot_file.name)
        
        # only the best one
        #cmd.append("-m 1")       
 
        # Excuting program
        cmd = " ".join(cmd)
        cmd2 = " ".join(cmd2)
        os.environ['DATAPATH'] = self.datapath
        os.system(cmd+"> /dev/null 2>/dev/null")
        os.system(cmd2+"> /dev/null 2>/dev/null")
        if seq_file is not None:
            seq_file.close()
        #parsing result
        out = open(dot_file.name).read().split('\n')
        sequence = out[1].strip()
        dot = out[2].strip()
        return sequence,dot
 
    def scorer(self,ct_name1,ct_name2):
        '''
        Use 'scorer' pogram to compare a predicted secondary structure to an accepted structure. It calculates two quality metrics, sensitivity and PPV
        
        :param ct_name1: The name of a CT file containing predicted structure data.
        :param ct_name2: The name of a CT file containing accepted structure data, can only store one structure.
        :return: sensitivity, PPV, number of the best predicted structure.
        
        Example:

        >>> ct_name1 = "temp_prediction.ct"
        >>> ct_name2 = "temp_accept.ct"
        >>> from RNAstructure import RNAstructure
        >>> RNA_prog = RNAstructure(exe_path="/home/yu68/Software/RNAstructure/exe/")
        >>> sensitivity, PPV, Number = RNA_prog.scorer(ct_name1,ct_name2)
        '''
        cmd = [os.path.join(self.exe_path,"scorer")]
        cmd.append(ct_name1)
        cmd.append(ct_name2)
        output_file = tempfile.NamedTemporaryFile(mode='w')
        cmd.append(output_file.name)
        cmd=" ".join(cmd)
        os.system(cmd+"> /dev/null 2>/dev/null")
        out=open(output_file.name).read().split("\n")
        n=0
        maximum=0.0
        sensitivity = 0.0
        PPV = 0.0
        number = 1
        for l in out:
            if l.startswith("Score"):
                n+=1
                continue
            if l.startswith("Sensitivity"):
                o = re.search(r'Sensitivity: (?P<portion>.*) = (?P<percent>.*)%',l)
                s = float(o.groupdict()['percent'])/100
            if l.startswith("PPV"):
                o = re.search(r'PPV: (?P<portion>.*) = (?P<percent>.*)%',l)
                p = float(o.groupdict()['percent'])/100
                if s+p > maximum:
                    sensitivity = s
                    PPV = p
                    maximum = s+p
                    number=n
        return sensitivity, PPV, number


# dot to block
def dot2block(dot_string,name="Default"):
    '''
    convert dot format of RNA secondary structure into several linked blocks

    :param dot_string: the dot format of RNA secondary structure
    :param name: name of the RNA
    :return: A list of all stems, each stem is a dictionary with 'source' and 'target'

    Example:
    
    >>> from RNAstructure import dot2block
    >>> stems = dot2block("(((((...)))...(((...)))..))","RNA_X")
    >>> print stems
    [{'source': {'start': 2, 'chr': 'test', 'end': 4}, 'target': {'start': 8, 'chr': 'test', 'end': 10}}, {'source': {'start': 14, 'chr': 'test', 'end': 16}, 'target': {'start': 20, 'chr': 'test', 'end': 22}}, {'source': {'start': 0, 'chr': 'test', 'end': 1}, 'target': {'start': 25, 'chr': 'test', 'end': 26}}]    
    
    '''
    posRegister=[]
    pairRegister=[]
    stems=[]
    def dumpPairRegister(pairRegister,stems,name):
        #print "source:%d-%d; target:%d-%d"%(pairRegister[-1][0],pairRegister[0][0],pairRegister[-1][1],pairRegister[0][1])
        stems.append({"source":{"chr":name,"start":pairRegister[-1][0],"end":pairRegister[0][0]},"target":{"chr":name,"start":pairRegister[0][1],"end":pairRegister[-1][1]}})
        pairRegister=[]
        return pairRegister,stems

    for i in range(len(dot_string)):
        char = dot_string[i]
        if char==".":  continue
        if char=="(":
            posRegister.append(i)
        if char==")":
            first = posRegister.pop()
            pair=[first,i]
            if (len(pairRegister)>0):
                if (pairRegister[-1][0]-first>1) or (i-pairRegister[-1][1]>1):
                    pairRegister,stems=dumpPairRegister(pairRegister,stems,name);
            pairRegister.append(pair)
    
    if (len(pairRegister)>0):
        pairRegister,stems=dumpPairRegister(pairRegister,stems,name);

    return stems 
