from Bio import SeqIO
from Bio.AlignIO import read
from Bio.Phylo.TreeConstruction import *
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

# Running Clustalw to align the sequences
clustalw_cline = ClustalwCommandline("clustalw2", infile="C:\\Users\\ricar\Downloads\\1-200613_S7_L001_R1_001.fastq\\1-200613_S7_L001_R1_001.fastq")
stdout, stderr = clustalw_cline()

# Reading the aligned sequences
aligned_sequences = AlignIO.read("C:\\Users\\ricar\Downloads\\1-200613_S7_L001_R1_001.fastq\\1-200613_S7_L001_R1_001", "clustal")

#registros = (read("C:\\Users\\ricar\Downloads\\1-200613_S7_L001_R1_001.fastq\\1-200613_S7_L001_R1_001.fastq","fastq"))

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aligned_sequences)
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)
Phylo.draw(tree)
Phylo.write(tree, "C:\\Users\\ricar\Downloads\\1-200613_S7_L001_R1_001.fastq\\tree.xml", "phyloxml")
