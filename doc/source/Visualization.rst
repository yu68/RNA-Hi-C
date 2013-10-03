====================================
Visualization of RNA-RNA interactome
====================================

Prerequirement
--------------

This program is powered by `RCircos <http://cran.r-project.org/web/packages/RCircos/index.html>`_.

Required R packages (our program will check for the presence of these packages and install/load them automatically if not present):
 
 * argparse, RCircos, biovizBase, rtracklayer 

The program also require a python script "bam2tab.py" (already in /bin/ folder) to call coverage from `BAM2X <https://github.com/nimezhu/bam2x/blob/master/scripts/bed2tab.py>`_

Run the program to generate visualization
-----------------------------------------
.. index:: Plot_Circos.R

We will use the script "Plot_Circos.R" for this purpose. ::
  
  usage: Plot_Circos.R [-h] [-g GENOME] [-b BIN] [-o OUTPUT]
                     interaction part1 part2

  positional arguments:
    interaction           the interaction file,[required]
    part1                 aligned BAM file for part1,[required]
    part2                 aligned BAM file for part2,[required]

  optional arguments:
    -h, --help            show this help message and exit
    -g GENOME, --genome GENOME
                          genome information, choice: mm9/mm10/hg19 et.al.,
                          [default: mm9]
    -b BIN, --bin BIN     window size for the bins for coverage calling, [default: 100000.0]
    -o OUTPUT, --output OUTPUT
                          output pdf file name, [default: Interactome_view.pdf]

.. note::
  
  part1, part2 BAM files are the ones generated from :ref:`Step5:Stitch-seq_Aligner.py<Step5>` of the pipeline; Interaction txt file is the output of :ref:`Step6:Select_strongInteraction_pp.py<Step6>`.


Example of result graph
-----------------------

Example code: ::
  
  Rscript Plot_Circos.R GGCG_interaction_clusters.txt 
  sort_Paired1_fragment_GGCG.bam sort_Paired2_fragment_GGCG.bam 
  -b 100000 -o Interactome_GGCG.pdf

Result figure:

.. image:: Interactome_GGCG.JPG

Explanation:

.. raw:: html

 <ul>
 <li>The <font color="#763a7a"> purple </font>track right inside chromatin cytoband ideogram is the coverage of part1 (the first genomic regions connected with linker sequences)  of this sample.</li> 
 <li>The <font color="#0288ad"> light blue </font>track next is the coverage of part2 (the other genomic regions connected with linkers). </li>
 <li>The <font color="red">inner </font>links are the strong interactions after removing rRNA. colors represent the confidence of the interaction (the ones with lower p-values are stronger) </li>
 </ul>
