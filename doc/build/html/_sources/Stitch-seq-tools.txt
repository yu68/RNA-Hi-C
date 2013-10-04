============================================
RNA-Stitch-seq-tools |version| documentation
============================================

Overview
========

**RNA-Stitch-seq-tools** is a set of bioinformatic tools for analysis of a novel DNA sequencing based technology to detect RNA-RNA interactome and RNA-chromatin interactome (RNA-chromatin interactome is coming soon).

Below is a illustration for the experimental design of this new technology

.. image:: exp.jpg
  :align: center

.. seealso:: 

  Offline documentation.

  Download a copy of RNA-stitch-seq-tools documentation:

   * `PDF </pdf/latest/RNA-Stitch-seq-tools.pdf>`_

Installation
============

step 1: Install the dependent prerequisites:
--------------------------------------------

1. Python libraries [for python 2.x]:

  * `Biopython <http://biopython.org/wiki/Main_Page>`_
  * `Pysam <https://code.google.com/p/pysam/>`_
  * `BAM2X <http://bam2xwiki.appspot.com/Welcome>`_
  * `Numpy <http://www.numpy.org/>`_, `Scipy <http://www.scipy.org/scipylib/index.html>`_
  * `Parallel python <http://www.parallelpython.com/>`_ (Only for ``Select_strongInteraction_pp.py``)
  * `PyCogent <http://pycogent.org/>`_ (for annotation of RNA types) [see note] 

2. The `Boost.Python <http://www.boost.org/doc/libs/1_54_0/libs/python/doc/index.html>`_ C++ library

3. Other softwares needed:

  * `Bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_ (not Bowtie 2)
  * `samtools <http://samtools.sourceforge.net/>`_
  * `NCBI blast+ <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ (use blastn)

.. note::

  the Annotation feature need the development version of PyCogent (`install instruction <http://pycogent.org/install.html#to-use-the-development-version-of-pycogent>`_). Since we need the getTranscriptByStableId function which is described `here <https://github.com/pycogent/pycogent/issues/21>`_.


Step 2: Download the package
----------------------------

Clone the package from GitHub::

  git clone http://github.com/yu68/stitch-seq.git


Step 3: Add library source to your python path
----------------------------------------------

Add these lines into your ~/.bash_profile or ~/.profile ::

  Location="/path/of/RNA-Stitch-seq-tools" # change accordingly
  export PYTHONPATH="$Location/src:$PYTHONPATH"
  export PATH="$PATH:$Location/bin"


Support
=======

For issues related to the use of RNA-Stitch-seq-tools, or if you want to **report a bug or request a feature**, please contact Pengfei Yu <p3yu at ucsd dot edu>

