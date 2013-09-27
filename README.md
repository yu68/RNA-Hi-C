stitch-seq
==========

## program fror stitch-seq ##

### pipeline (start from pair-end fastq file): ###
1. Remove duplicates
2. Split libraries based on barcode.txt
3. Recover fragment for each library
4. Split partners and classify different types of fragments.
5. Algin to the genome for both parts of "paired" fragments.
6. Determine the RNA types of different parts within fragment.
7. (Alternative) Find liner sequences within the library.


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


#### Notes: ####
1. the Annotation feature need the development version of PyCogent [install instruction](http://pycogent.org/install.html#to-use-the-development-version-of-pycogent). Since we need the getTranscriptByStableId function which is described [here](https://github.com/pycogent/pycogent/issues/21).
