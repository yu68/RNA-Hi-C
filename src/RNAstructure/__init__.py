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
        >>> RNA_prog = RNAstructure(exe_path="~/Software/RNAstructure/exe/")
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

        :param seq1,seq2: two DNA/RNA sequences as string
        :param dna: boolean input. Specify then DNA parameters are to be used
        :returns: minimum binding energy, (unit: kCal/Mol)

        Example:

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
