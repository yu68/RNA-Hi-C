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

Starting from the raw pair-end sequencing data, PCR duplicates should be removed as the first step if both the 10nt random indexes and the remaining sequences are exactly the same for two pairs.

.. _step2:

Step 2: Split library based on barcode.txt.
-------------------------------------------

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
