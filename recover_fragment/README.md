Recover fragment
================

* recover\_fragment-2.py detect _first three_ types
* recover\_fragment-2.py detect all _four_ types

```
* type1 [short]:

                                                        (reverse primer)
       forward reads:                      XXXX...XXXXNAGATCGGAAGAGCGGTTCAG
                                           ||||...||||
       reverse reads: TGTGCTGCGAGAAGGCTAGANXXXX...XXXX
                       (forward primer)
```
```
* type2 [long]:

       forward reads: XXXXX...XXXXXXXXXXX...XXXX
                                     ||||...||||
       reverse reads:                XXXX...XXXXXXXXXXX...XXXX
```
```
* type3 [even_long]:

       forward reads: XXXXX...XXXX
                                     
       reverse reads:                  XXXX...XXXXX
```
```
* type4 [wierd]:
                                                           (reverse primer)
       forward reads:                      XXXX...XXXXNAGATCGGAAGAGCGGTTCAG

       reverse reads: XXXXX...XXXX (no forward primer)

       # or the other way
```
