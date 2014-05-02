.. RNA-Stitch-seq-tools documentation master file, created by
   sphinx-quickstart on Fri Sep 27 11:30:10 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to RNA-Stitch-seq-tools's documentation!
================================================

**Contents:**

.. toctree::
   :maxdepth: 3

   Stitch-seq-tools
   
   Analysis_pipeline

   Visualization.rst

   Other_api.rst


.. note::

   RNA-Stitch-seq-tools benifits a lot from BAM2X, a convenient python interface for most common NGS datatypes. `Try BAM2X now <http://bam2xwiki.appspot.com/>`_!


Updates
-------

New features in 0.3.1:
  * Add "--release" option in ":ref:`split_partner.py<Step4>`" function. Allow a Type3 read-pair considered to be a "Paired" chimeric fragment even linker does not show up.
  * Fix bugs on ":ref:`Select_strongInteraction_pp.py<Step6>`" function when the number of mapped pairs is low and some chromosomes don't have any mapped read in part1 or part2.
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

