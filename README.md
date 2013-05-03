stitch-seq
==========

## program fror stitch-seq ##

### pipeline (start from pair-end fastq file): ###
1. Remove duplicates
2. Remove phiX (using bowtie)
3. Split libraries based on barcode.txt
4. Recover fragment for each library
5. Do BLAST find linked miRNA and mRNA
6. Find liner sequences within the library


we can also do bed annotation on different cis-features.




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


* type1:
```
                                                        (reverse primer)
       forward reads:                      XXXX...XXXXNAGATCGGAAGAGCGGTTCAG
                                           ||||...||||
       reverse reads: TGTGCTGCGAGAAGGCTAGANXXXX...XXXX
                       (forward primer)
```
* type2:
```
       forward reads: XXXXX...XXXXXXXXXXX...XXXX
                                     ||||...||||
       reverse reads:                XXXX...XXXXXXXXXXX...XXXX
```

