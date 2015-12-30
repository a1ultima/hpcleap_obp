###########
# Imports #
###########
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
import os
from Bio import AlignIO
import subprocess

#############
# Functions #
#############

def gen_multiple_sequence_alignment(geneid):

	""" Run a multiple sequence alignment (MSA) on all orthologous sequences within a FASTA file, named after <geneid>

	Args: 
		geneid = 'AGAP002191'
	Exceptions:
		None
	Retuns:
		None 
	Usage:
		gen_multiple_sequence_alignment('AGAP002191')

	"""

	cline = ClustalwCommandline(	"clustalw", 
									infile= "../data/agam_orthologues_filtered/"+geneid+".fasta",\
									outfile="../data/alignments/"+geneid+".aln" )

	stdout, stderr = cline()

def get_obp_ids(filepath):

	""" Returns a list of VB Agam OBP Ids to use as search queries for FASTA sequences of orthologous OBPs

	Args: 
		filepath =  "../data/vb_agam_obp_geneIDs.txt"
	Exceptions:
		None
	Retuns:
		None 
	Usage:
		vb = get_obp_ids("../data/vb_agam_obp_geneIDs.txt")

	"""

	fi = open("../data/vb_agam_obp_geneIDs.txt")
	vb = []
	while True:
		line = fi.readline().rstrip()
		if line == "":
			break
		vb.append(line)
	fi.close()
	return vb

def alignment_clustal_to_fasta(alignment_file_path):

	""" Converts an alignment file *.aln from "clustal" format to "fasta"

	Args: 
		alignment_file_path = "../data/alignments/AGAP004433.aln"
	Exceptions:
		None
	Retuns:
		None 
	Usage:
		alignment_clustal_to_fasta("../data/alignments/AGAP004433.aln")

	"""

	alignment = AlignIO.read(alignment_file_path, "clustal")
	AlignIO.write([alignment], alignment_file_path.replace(".aln","_fasta.aln"), "fasta")

def gen_phylogeny(alignment_fasta_path):

	""" Generate a phylogeny using FastTree (http://www.microbesonline.org/fasttree/), taking a fasta format alignment as input

	Args: 
		None 
	Exceptions:
		None
	Retuns:
		None 
	Usage:
		gen_phylogeny("../data/alignments/AGAP004433_fasta.aln")

	"""
	try:
		subprocess.call('FastTree '+alignment_fasta_path+' > '+alignment_fasta_path.replace("alignments","phylogenies").replace("_fasta.aln",".tre"), shell=True)
	except OSError:
		raise OSError('It seems FastTree is not installed... get it here: http://www.microbesonline.org/fasttree/')


########
# Main #
########

#
# Prepare file system
#
print "preparing directories..."
if not os.path.isdir("../data/"):
	os.makedirs("../data/")
if not os.path.isdir("../data/alignments"):
	os.makedirs("../data/alignments")
if not os.path.isdir("../data/phylogenies"):
	os.makedirs("../data/phylogenies")

# #
# # Make MSAs w/ ClustalW
# #
# print "Generating multiple sequence alignments of orthogroups..."
# vb = get_obp_ids("../data/vb_agam_obp_geneIDs.txt")
# for gene in vb:
# 	print "\t"+gene
# 	gen_multiple_sequence_alignment(gene)


#
# Convert alignment to FASTA
#
print "Convertig MSAs to fasta..."
for gene in vb:
	print "\t"+gene
	alignment_clustal_to_fasta("../data/alignments/"+gene+".aln") # converted alignments: ../data/alignments/*_fasta.aln

#
# Run FastTree to generate Phylogenies 
#
print "Generating Phylogenies..."
for gene in vb:
	print "\t"+gene
	gen_phylogeny("../data/alignments/"+gene+"_fasta.aln") # phylogenies: ../data/phylogenies/*

