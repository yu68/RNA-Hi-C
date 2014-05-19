.. _Resource:

==========================================================
Resources of strong interactions from two mouse cell types
==========================================================

Description of different samples
================================

E14_WP_1
--------

:Cell line: ESC E14
:Barcode: ACCT
:Experimental Details: Actively growing E14 cells were UV irradiated (254 nm) at 200mJ/cm 
  to crosslink proteins to interacting RNAs. After cell lysis, we trim down RNAs into 
  1000-2000 nt using RNase I and remove DNA by TURBO DNase. To recover RNAs bound to 
  RNA-binding proteins, we biotin-labeled them with EZ-Link Iodoacetyl-PEG2-Biotin from 
  Pierce. RNA-protein complexes were next immobilized on Streptavidin-coated beads. The 
  beads are then saturated with free biotin, preventing it from interfering with following 
  ligation with biotin-tagged linker. A biotin-tagged RNA linker was ligated to the 5’-end 
  of immobilized RNAs. Proximity ligation was then carried out under diluted conditions 
  while the RNA-protein complexes are still bound on bead. After RNA purification by 
  Proteinase K and phenol-chloroform extraction, we specifically removed the unligated 
  biotin by first anneal a complementary DNA oligo to the biotin-tagged linker by using 
  the annealing protocol: 70oC for 5 min, 25oC for 20 min, T7 Exo for 30 min. Exonuclease 
  T7 was added to remove terminal unligated biotin at the double-stranded RNAlinker-DNAoligo 
  hybrid. T7 Exonuclease acts in the 5' to 3' direction, catalyzing the removal of 5' 
  mononucleotides from duplex DNA and RNA/DNA hybrids in the 5’ to 3’ direction. The resulted 
  RNAs were fragmented again into ~200 nt using RNase III RNA Fragmentation Module from NEB 
  (1ul of RNase III in 6 min at 37C). The RNAs were purified by column and ligated with 
  sequencing adapter, then reverse-transcribed and PCR for library construction. We applied 
  an rRNA removal step after constructing cDNA by using an rRNA removal protocol based on 
  the Duplex-Specific thermostable nuclease (DSN) enzyme using the protocol recommended by 
  `Illumina <http://supportres.illumina.com/documents/myillumina/7836bd3e-3358-4834-b2f7-80f80acb4e3f/dsn_normalization_sampleprep_application_note_15014673_c.pdf>`_. 
  The constructed cDNAs were quality-checked by Bioanalyzer. The cDNAs were next 
  subjected paired-end sequencing on HiSeq-2500 platform.
:Linker: 
  *mL5*: 5' - rCrUrA rG/iBiodT/rA rGrCrC rCrArU rGrCrA rArUrG rCrGrA rGrGrA - 3'

E14_WP_2
--------

:Cell line: ESC E14
:Barcode: GGCG
:Experimental Details: Same as E14_WP_1 but this time rRNA removal was performed right after 
  Proteinase K and phenol-chloroform treatment using the GeneRead rRNA Depletion Kit by 
  Qiagen. Furthermore, the annealing of RNA linker and complementary DNA oligo was changed 
  into: denature for 90 s at 90°C, and then anneal at -0.1°C/s to 25 °C and then incubate 
  for 25 min at 25 °C. Since after rRNA depletion the amount of RNA remained was less than 
  that obtained from E14 #1, we reduced the duration of RNA fragmentation by RNase III from 
  6 min to 3 min. However, this reduction in RNase III treatment led to large fragments than 
  desirable. 
:Linker:
  *mL5*: 5' - rCrUrA rG/iBiodT/rA rGrCrC rCrArU rGrCrA rArUrG rCrGrA rGrGrA - 3'

E14_WP_3
--------

:Cell line: ESC E14
:Barcode: AATG
:Experimental Details: To detect interactions between RNAs that are not bound to the same protein 
  but to interacting proteins, we used formaldehyde in conjunction with a second crosslinker, EGS. 
  The combination of formaldehyde and EGS crosslinks both RNA-protein and protein-protein 
  interactions thereby maximize the detection of RNA-RNA interactions that are facilitated by 
  interacting proteins. Actively grown E14 cells was crosslinked with 1.5 mM of freshly prepared 
  EthylGlycol bis(SuccinimidylSuccinate) (EGS)) for 45 minutes at room temperature and then 1% of 
  formaldehyde for 10 minutes also at room temperature. Since crosslinking by formaldehyde makes 
  the cells very rigid and less amenable to be broken down lysis buffer. Therefore, we utilized 
  sonication to fragment the protein-bound RNA into ~1000 nt size range. The remaining steps were 
  performed similarly to E14_WP_2. 
  
  Another main difference between this sample and other samples is that we didn't remove the nuclei, 
  thus effectively including RNA-RNA interactions inside the nucleus into the sample. In other 
  samples, only the cytoplasm was enriched.
:Linker:
  *mL5*: 5' - rCrUrA rG/iBiodT/rA rGrCrC rCrArU rGrCrA rArUrG rCrGrA rGrGrA - 3'

MEF_WP_1
--------

:Cell line: MEF
:Barcode: GGCG
:Experimental Details: We irradiated actively grown 1E8 MEF cells (254 nm). This time, the RNAs 
  were fragmented into 300nt size range. RNase III fragmentation was also modified accordingly 
  to adjust for smaller amount of RNAs: instead of adding 1ul of RNase III, we added only 0.5uL 
  of RNase III and then incubated at 37C for 5 min. The subsequent steps were performed using 
  the same procedure as E14_WP_2.
:Linker:
  *mL5*: 5' - rCrUrA rG/iBiodT/rA rGrCrC rCrArU rGrCrA rArUrG rCrGrA rGrGrA - 3'


Resources of Strong Interactions
================================

From merged data of E14_WP_1 and E14_WP_2:
------------------------------------------
`Download <http://systemsbio.ucsd.edu/RNA-Hi-C/Data/ACCT_GGCG_interaction_clusters.xlsx>`_

From E14_WP_3 dual crosslinking:
--------------------------------
`Download <http://systemsbio.ucsd.edu/RNA-Hi-C/Data/AATG_interaction_clusters.xlsx>`_

From MEF_WP_1 sample:
---------------------
`Download <http://systemsbio.ucsd.edu/RNA-Hi-C/Data/GGCG_MEF_interaction_clusters.xlsx>`_

