.. RNA-Stitch-seq-tools documentation master file, created by
   sphinx-quickstart on Fri Sep 27 11:30:10 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RNA Hi-C tools's documentation!
==========================================

**Contents:**

.. toctree::
   :maxdepth: 3

   RNA-Hi-C-tools
   
   Analysis_pipeline

   Visualization.rst

   Other_api.rst

   Data_Resources

.. note::

   RNA Hi-C tools benifits a lot from BAM2X, a convenient python interface for most common NGS datatypes. `Try BAM2X now <http://www.bam2x.net/>`_!


Updates
-------

2014-11-7:
  * Changed the mapping mechanism of RNA pairs: RNA pairs are now decoupled first, and mapped to several references (order specified from the parameters) including miRNA, genome, transcripts and other sequences, in current workflow, the order is miRNA->transcripts->other RNA sequences->genome. After mapping the two fragments from the same RNA pair are linked back and pairs with both fragments mapped will be reported. ":ref:`Stitch-seq_Aligner.py<Step5>`"
  * Mapping parameters has been changed so that if one or both fragment is mapped to multiple sequencing in one reference, all match from the best stratum will be reported with a "multimap" tag attached in the result.
  * If the RNA fragment is mapped to the reverse-complement strand of a sequence in the reference (e.g. mapped to the reverse-complement strand of an RNA sequence), this mapping will be discarded. A parameter is added to allow such mapping results be retained and marked out in the final result.

2014-10-27:
  * Add new script to detect potential splicing intermediates from snoRNA-mRNA interactions ":ref:`snoRNA_mRNA_statistics.py<snoRNA_mRNA>`"

2014-7-15:
  * Update ":ref:`RNA_structure_prediction.py<Structure>`" function to allow output of JSON files for predicted structure and refined structure (predicted structure after providing single-strand offset information). The JSON output can be uploaded into `RNA2D-browser <http://circos.zhu.land/>`_ (developed by `Xiaopeng Zhu <https://github.com/nimezhu>`_ ) to show the Circos view of secondary structure and digested location distribution.
  * Add an API function to convert dot format of RNA secondary structure into several linked blocks. see ":ref:`dot2block<dot2block>`"

2014-6-27:
  * new strong interaction list added based on whole RNA annotation using a FDR cutoff, and using ES-indirect (dual crosslinking) sample as control. See: :ref:`resources<SIFDR>`, update in ":ref:`Select_strongInteraction_pp.py<Step6>`" as well.

2014-5-15:
  * Add result :ref:`resources<Resource>` for identified strong interactions in mouse E14 cells and MEF cells.
  * New function to generate heatmap for intra-RNA interactions: :ref:`Plot_interaction_heatmap.py<VisualizationHeatmap>`.

2014-05-14:
  * Add new function to find overlap between two interaction sets based on their RNA annotations, see: :ref:`intersectInteraction_genePair.R<intersectiongene>`.
  * Allow input of two genomic regions to visualize local interactions using ``-r`` option in ":ref:`Plot_interaction.py<plotInteraction>`" function

2014-05-11:
  * Add new function to show enrichment of different types of interactions: :ref:`Interaction_type_enrichment.R<VisualizationEnrich>`.

Version 0.3.2 (2014-05-07):
  * change the name into RNA-Hi-C

2014-05-06:
  * In ":ref:`Select_strongInteraction_pp.py<Step6>`" function, now annotations are updated after doing clustering and for strong interaction. The indexing of annotation files may take some time.
  * New ":ref:`RNA_structure_prediction.py<Structure>`" function to refine RNA structure prediction based on empirical offset of free energies for single strand nucleotide.

New features in 0.3.1 (2014-05-02):
  * Add "--release" option in ":ref:`split_partner.py<Step4>`" function. Allow a Type3 read-pair considered to be a "Paired" chimeric fragment even linker does not show up.
  * Fix bugs in ":ref:`Select_strongInteraction_pp.py<Step6>`" function when the number of mapped pairs is low and some chromosomes don't have any mapped read in part1 or part2.
  * Add bowtie 2 option and Unique-align option in ":ref:`Stitch-seq_Aligner.py<Step5>`" function. 
  * Different colors for different types of interactions in the :ref:`visualization of interactome<VisualizationGlobal>`. 
  * New API for folding energies of two RNA molecules, see ":ref:`RNAstructure<rnafold>`". 
  * Allow permutation-based strategies to calculate the p-value for the overlap between two independent interaction sets in ":ref:`intersectInteraction.py<intersection>`" function

 
New features in 0.2.2:
  * ":ref:`Plot_interaction.py<plotInteraction>`" function to plot local RNA-RNA interactions. 
  * ":ref:`intersectInteraction.py<intersection>`" function to call overlap between two independent interaction sets.  

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

