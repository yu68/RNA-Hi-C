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

    def Fold(self,seq=None,ct_name=None):
        '''
        Use "Fold" program to predict the secondary structure and output dot format.

        :param seq: one DNA/RNA sequence as string, or existing fasta file name
        :param ct_name: specify to output a ct file with this name, otherwise store in temp, default: None
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
        cmd2.append('1')
        cmd2.append(dot_file.name)
        # Excuting program
        cmd = " ".join(cmd)
        cmd2 = " ".join(cmd2)
        print cmd
        print cmd2
        os.environ['DATAPATH'] = self.datapath
        print os.environ['DATAPATH']
        os.system('echo $DATAPATH')
        os.system(cmd+"")
        os.system(cmd2+"> /dev/null 2>/dev/null")
        if seq_file is not None:
            seq_file.close()
        #parsing result
        out = open(dot_file.name).read().split('\n')
        sequence = out[1].strip()
        dot = out[2].strip()
        return sequence,dot
