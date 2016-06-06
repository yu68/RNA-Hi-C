RNA Hi-C
========

## program for RNA-RNA interactome and RNA-chromatin interactome ##

### Documentation ###
[HTML](http://pengfeiyu.net/RNA-HiC-tools/)

### pipeline (start from pair-end fastq file): ###
1. Remove duplicates  `rm_dup\remove_dup_PE.py`
2. Split libraries based on barcode.txt `split_library_pairend.py`
3. Recover fragment for each library `recoverFragment`
4. Split partners and classify different types of fragments. `split_partner.py`
5. Algn to the genome for both parts of "paired" fragments. `Stitch-seq_Aligner.py`
6. Determine the RNA types of different parts within fragment. `RNA_composition.py`
7. (Alternative) Find liner sequences within the library. `find_linker_new.py`
8. Determine strong interactions from output of step 5. `Select_strongInteraction_pp.py` for parallel computing; `Select_strongInteraction.py` regular.

we can also do bed annotation on different cis-features. `bed_annotation.py`




### Seven situation: ###

1. NNNXXXXNN - miRNA - AUC UGG UAA UCC GUA UAA AGU AUG UUG AUG UUC CAA UAA GCA GAU CAU GUU UUU UAA GCC GUC A - mRNA
2. NNNXXXXNN - miRNA - UAA GCA GAU CAU GUU UUU UAA GCC GUC A - mRNA
3. NNNXXXXNN - miRNA - AUC UGG UAA UCC GUA UAA AGU AUG UUG AUG UUC CAA - mRNA
4. NNNXXXXNN - miRNA - AUC UGG UAA UCC GUA UAA AGU AUG UUG AUG UUC CAA
5. NNNXXXXNN - UAA GCA GAU CAU GUU UUU UAA GCC GUC A - mRNA
6. NNNXXXXNN - mRNA
7. NNNXXXXNN - miRNA (less likely)

linkers:
  * UAA GCA GAU CAU GUU UUU UAA GCC GUC A
  * AUC UGG UAA UCC GUA UAA AGU AUG UUG AUG UUC CAA


###  recover fragment   ###

```
* type1:

                                                        (reverse primer)
       forward reads:                      XXXX...XXXXNAGATCGGAAGAGCGGTTCAG
                                           ||||...||||
       reverse reads: TGTGCTGCGAGAAGGCTAGANXXXX...XXXX
                       (forward primer)
```
```
* type2:

       forward reads: XXXXX...XXXXXXXXXXX...XXXX
                                     ||||...||||
       reverse reads:                XXXX...XXXXXXXXXXX...XXXX
```


#### Library Dependency #####
1. python libraries [python 2.x]:   
..* [Biopython](http://biopython.org/wiki/Main_Page)  
..* [Pysam](https://code.google.com/p/pysam/)  
..* [BAM2X](http://bam2xwiki.appspot.com/Welcome)  
..* numpy, scipy  
..* [parallel python](http://www.parallelpython.com/) (only for `Select_strongInteraction_pp.py`)  
..* [PyCogent](http://pycogent.org/) (for annotation of RNA types) [see Notes]  

2. [Boost.Python](http://www.boost.org/doc/libs/1_54_0/libs/python/doc/index.html)  

