
# Tool namespace — makes 'tools' a Python package.
#
# Upstream (QC → alignment → assembly):
#   qc.py            FastQC + MultiQC command builders
#   star.py          STAR index + align command builders
#   trinity.py       Trinity assembly command builder
#
# Downstream (counting → DE → annotation → domains):
#   featurecounts.py featureCounts (Subread) command builder
#   pydeseq2.py      PyDESeq2 differential expression helpers
#   blast.py         BLASTp + BLASTx command builders
#   hmmer.py         HMMER hmmscan (Pfam) command builder
#   prosite.py       EMBOSS patmatmotifs command builder
