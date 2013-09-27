=================
Analysis pipeline
=================

Overview
========

The next generation DNA sequencing based technology utilize RNA proximity ligation to transfrom RNA-RNA interactions into chimeric DNAs. Through sequencing and mapping these chimeric DNAs, it is able to achieve high-throughput mapping of nearly entire interaction networks. RNA linkers were introduced to mark the junction of the ligation and help to split the chimeric RNAs into two interacting RNAs.
This bioinformatic pipeline is trying to obtain the strong interactions from raw fastq sequencing data. The major steps are:

* :ref:`Step1`
* :ref:`Step2`
* :ref:`Step3`
* :ref:`Step4`
* :ref:`Step5`
* :ref:`Step6`

Other functions:

1. :ref:`RNA_types`
2. :ref:`find_linker`

Pipeline
========

.. _step1:

Step 1: Remove PCR duplicates.
------------------------------

Starting from the raw pair-end sequencing data, PCR duplicates should be removed as the first step if both the 10nt random indexes and the remaining sequences are exactly the same for two pairs. It is achieved by ``remove_dup_PE.py`` ::

  usage: remove_dup_PE.py [-h] reads1 reads2

  Remove duplicated reads which have same sequences for both forward and reverse
  reads. Choose the one appears first.

  positional arguments:
    reads1      forward input fastq/fasta file
    reads2      reverse input fastq/fasta file

  optional arguments:
    -h, --help  show this help message and exit

  Library dependency: Bio, itertools

The program will generate two fastq/fasta files after removind PCR duplicates and report how many read pairs has been removed. The output are prefixed with 'Rm_dupPE'

.. note::

  One pair is considered as a PCR duplicate only when the sequences of both two ends (including the 10nt random index) are the exactly same as any of other pairs.

.. _step2:

Step 2: Split library based on barcode.txt.
-------------------------------------------
After removing PCR duplicates, the libraries from different samples are separated based on 4nt barcodes in the middle of random indexes ("RRRBBBBRRR"; R: random, B: barcode). It is implemented by ``split_library_pairend.py`` ::

  usage: split_library_pairend.py [-h] [-f | -q] [-v] [-b BARCODE]
                                  [-r RANGE [RANGE ...]] [-t] [-m MAX_SCORE]
                                  input1 input2

  Example: split_library_pairend.py -q Rm_dupPE_example.F1.fastq 
           Rm_dupPE_example.R1.fastq -b barcode.txt

  positional arguments:
    input1                input fastq/fasta file 1 for pairend data (contain
                          barcodes)
    input2                input fastq/fasta file 2 for pairend data

  optional arguments:
    -h, --help            show this help message and exit
    -f, --fasta           add this option for fasta input file
    -q, --fastq           add this option for fastq input file
    -v, --version         show program's version number and exit
    -b BARCODE, --barcode BARCODE
                          barcode file
    -r RANGE [RANGE ...], --range RANGE [RANGE ...]
                          set range for barcode location within reads,default is
                          full read
    -t, --trim            trim sequence of 10nt index
    -m MAX_SCORE, --max_score MAX_SCORE
                          max(mismatch+indel) allowed for barcode match,
                          otherwise move reads into 'unassigned' file
                          default: 2.

  Library dependency: Bio

Here is a example for barcode.txt ::
  
  ACCT
  CCGG
  GGCG

The output of this script are several pairs of fastq/fasta files prefixed with the 4nt barcode sequences, together with another pair of fastq/fasta files prefixed with 'unassigned'.

For example, if the input fastq/fasta files are ``Rm_dupPE_example.F1.fastq`` and ``Rm_dupPE_example.R1.fastq``, and the barcode file is the same as above, then the output files are ::

  ACCT_Rm_dupPE_example.F1.fastq
  ACCT_Rm_dupPE_example.R1.fastq
  CCGG_Rm_dupPE_example.F1.fastq
  CCGG_Rm_dupPE_example.R1.fastq
  GGCG_Rm_dupPE_example.F1.fastq
  GGCG_Rm_dupPE_example.R1.fastq
  unassigned_Rm_dupPE_example.F1.fastq
  unassigned_Rm_dupPE_example.F1.fastq

.. _step3:

Step 3: Recover fragments for each library.
-------------------------------------------

.. _step4:

Step 4: Split partners and classify different types of fragments.
-----------------------------------------------------------------

.. _step5:

Step 5: Align both parts of "Paired" fragment to the genome.
------------------------------------------------------------

.. _step6:

Step 6: Determine strong interactions.
--------------------------------------


Other functions
===============

.. _RNA_types:

Determine the RNA types of different parts within fragments.
------------------------------------------------------------

.. _find_linker:

Find linker sequences within the library.
-----------------------------------------
