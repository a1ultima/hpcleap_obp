###########
# Imports #
###########

from collections import Counter
from Bio import SeqIO
from splinter import Browser
from splinter import exceptions
import pdb
import os

from time import sleep

#############
# Functions #
#############

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

def get_orthologue_seqs(gene_ids,source="web"):

	""" Get orthologue sequences from VB for each gene id

	Args:
		gene_ids = ['AGAP012658','AGAP012714','AGAP008278','AGAP008279','AGAP028120','AGAP008284','AGAP008282','AGAP008283','AGAP008281','AGAP008280','AGAP013182','AGAP003309','AGAP001189','AGAP001189','AGAP002025','AGAP002188','AGAP002905','AGAP002189','AGAP003307','AGAP003307','AGAP012319','AGAP004433','AGAP003306','AGAP005208','AGAP008398','AGAP010409','AGAP012318','AGAP012320','AGAP012321','AGAP012323','AGAP012323','AGAP012325','AGAP012331','AGAP001409','AGAP001409','AGAP000638','AGAP000638','AGAP000640','AGAP000640','AGAP000641','AGAP000642','AGAP000643','AGAP000644','AGAP000580','AGAP000580','AGAP002190','AGAP002190','AGAP010489','AGAP002191','AGAP002191','AGAP005182','AGAP005182','AGAP009065','AGAP009065','AGAP009402','AGAP009402','AGAP010648','AGAP010648','AGAP010650','AGAP010650','AGAP007289','AGAP009629','AGAP006079','AGAP006079','AGAP006081','AGAP011368','AGAP003530','AGAP001556','AGAP000279','AGAP000278','AGAP012867','AGAP006368','AGAP006278','AGAP002556','AGAP010649','AGAP010649','AGAP011647','AGAP011647','AGAP012322','AGAP012324']

		mode = "web" # or: <path to directory containing fasta data> e.g. "/home/qiime/Desktop/hpcleap_obp/data/agam_orthologues"

	Usage:
		gene_to_orthoFasta = get_orthologue_seqs(['AGAP012658','AGAP012714','AGAP008278','AGAP008279','AGAP028120','AGAP008284','AGAP008282','AGAP008283','AGAP008281','AGAP008280','AGAP013182','AGAP003309','AGAP001189','AGAP001189','AGAP002025','AGAP002188','AGAP002905','AGAP002189','AGAP003307','AGAP003307','AGAP012319','AGAP004433','AGAP003306','AGAP005208','AGAP008398','AGAP010409','AGAP012318','AGAP012320','AGAP012321','AGAP012323','AGAP012323','AGAP012325','AGAP012331','AGAP001409','AGAP001409','AGAP000638','AGAP000638','AGAP000640','AGAP000640','AGAP000641','AGAP000642','AGAP000643','AGAP000644','AGAP000580','AGAP000580','AGAP002190','AGAP002190','AGAP010489','AGAP002191','AGAP002191','AGAP005182','AGAP005182','AGAP009065','AGAP009065','AGAP009402','AGAP009402','AGAP010648','AGAP010648','AGAP010650','AGAP010650','AGAP007289','AGAP009629','AGAP006079','AGAP006079','AGAP006081','AGAP011368','AGAP003530','AGAP001556','AGAP000279','AGAP000278','AGAP012867','AGAP006368','AGAP006278','AGAP002556','AGAP010649','AGAP010649','AGAP011647','AGAP011647','AGAP012322','AGAP012324'],source="web")

	"""

	gene_to_orthoFasta = {} 	# store FASTA

	if source == "web":
		print "\tGathering OBP orthogroups from the: WEB"
		#
		# visit vectorbase, grabbing orthogroup fasta files
		#
		url_prefix = "https://www.vectorbase.org/Anopheles_gambiae/Gene/Compara_Ortholog/PepSequence?_format=Text;db=core;g="
		br = Browser()
		for gene in gene_ids:
			url = url_prefix+gene 
			br.visit(url)
			fasta = br.find_by_css("body > pre:nth-child(1)").first.text
			gene_to_orthoFasta[gene] = fasta 
			sleep(2)
		#
		# Write to file 
		#
		print "\tWriting FASTAs to: ../data/agam_orthologues/*"
		for gene in gene_to_orthoFasta.keys():
			fo = open("../data/agam_orthologues/"+gene+".txt","w")
			fo.write(gene_to_orthoFasta[gene])
			fo.close()
	else: 
		#
		# Grab fasta files from gene ids if available
		#
		# 	file names are given by: source + geneid + ".txt", read iteratively
		print "\tGathering OBP orthogroups from: "+source
		for gene in gene_ids:
			path = source + "/" + gene + ".fasta"
			try:
				fi = open(path,"r")
			except IOError:
				raise IOError('The file: '+path+' does not exist, try: source="web" ...')
			lines = fi.readlines()
			fasta = "".join(lines)
			gene_to_orthoFasta[gene] = fasta
			fi.close()

	return gene_to_orthoFasta

def filter_out_non_anopheline(geneid):

	""" Given a FASTA file, filters out all orthologs that are non-anopheline (see "species_keep")

	Args: 
		geneid = "AGAP002191"  
	Retuns:
		None 
	Usage:
		filter_out_non_anopheline("AGAP002191")

	"""

	species_keep = ['Anopheles_maculatus', 'Anopheles_merus', 'Anopheles_coluzzii', 'Anopheles_darlingi', 'Anopheles_funestus', 'Anopheles_dirus', 'Anopheles_stephensi', 'Anopheles_atroparvus', 'Anopheles_melas', 'Anopheles_gambiae', 'Culex_quinquefasciatus', 'Anopheles_arabiensis', 'Anopheles_minimus', 'Aedes_aegypti', 'Anopheles_albimanus', 'Anopheles_epiroticus', 'Anopheles_sinensis', 'Anopheles_christyi', 'Anopheles_farauti', 'Anopheles_quadriannulatus', 'Anopheles_culicifacies']

	#
	# Read fasta sequences
	#
	print "Reading orthologous FASTA sequences...."
	print "\tFind non-filetered FASTA sequences in: ../data/agam_orthologues/*"

	fi = open("../data/agam_orthologues/"+geneid+".fasta")

	fasta_sequences = SeqIO.parse(fi,'fasta')
	#fasta_sequences = SeqIO.to_dict(SeqIO.parse(fi,'fasta'))

	#
	# Remove "Non-Anopheline" orthologues (i.e. not in the 20 genomes paper)
	#

	sequences = list(fasta_sequences)

	sequences_anopheline = []
	species_full = []
	species_filtered = []

	for seq in sequences:
		# get the species of this orthologue
		gene_species = seq.name.split("-PA_")
		if len(gene_species)<2:
			gene_species = seq.name.split("-PB_")	
		species = gene_species[1]
		species_full.append(species)
		# collect if species is anopheline
		if species in species_keep:
			species_filtered.append(species)
			sequences_anopheline.append(seq)

	#
	# write fasta seq
	#
	fo = open("../data/agam_orthologues_filtered/"+geneid+".fasta", "w")

	SeqIO.write(sequences_anopheline, fo, "fasta")

	fi.close()
	fo.close()

def test_gene_overlap(gene_to_orthoFasta):

	""" Take the output of get_orthologue_seqs() to test which of the orthologous groups of OBPs share gene ids 

	Args:
		gene_to_orthoFasta = get_orthologue_seqs(vb,source="../data/agam_orthologues_filtered")
	Returns:
		n_overlapping_orthogroups 	# a count of overlapping groups of orthologous obps, by their gene ids
	Usage:
		n_overlapping_orthogroups = test_gene_overlap(gene_to_orthoFasta)


	"""

	gene_to_orthoGenes = {}

	for geneId in gene_to_orthoFasta.keys():

		genes = []
		fi = open("../data/agam_orthologues/"+geneId+".fasta","r")
		while True: 
			line = fi.readline()
			if line=="":
				break
			if ">" in line:
				gene_species = line.split("-PA_")
				if len(gene_species)<2:
					gene_species = line.split("-PB_")
				gene = gene_species[0][1:]
				genes.append(gene)
		fi.close()
		gene_to_orthoGenes[geneId] = genes

	count = 0

	for a,i in enumerate(gene_to_orthoGenes.keys()):
		for b,j in enumerate(gene_to_orthoGenes.keys()):
			intersection = set(gene_to_orthoGenes[i]) & set(gene_to_orthoGenes[j])
			if a != b:
				if len(intersection)>0:
					count +=1
					print i + ":::" + j + "::"+str(count)+":"+str(intersection)


	print "There are "+str(count/2)+" orthologous groups of OBPs that share the same geneIds..."

	return count/2

########
# Main #
########

#
# Prepare directories
#
print "preparing directories..."
if not os.path.isdir("../data/"):
	os.makedirs("../data/")
if not os.path.isdir("../data/agam_orthologues"):
	os.makedirs("../data/agam_orthologues")
if not os.path.isdir("../data/agam_orthologues_filtered"):
	os.makedirs("../data/agam_orthologues_filtered")

#
# Read VB Agam OBP gene Ids
#
#mobp=['AGAP000278','AGAP000279','AGAP001556','AGAP012714','AGAP012867','AGAP006368','AGAP003530','AGAP013182','AGAP012658','AGAP012659','AGAP007283','AGAP006759','AGAP012324','AGAP012322','AGAP002556','AGAP007282','AGAP007281','AGAP009629','AGAP006760','AGAP006074','AGAP011368','AGAP011367','AGAP006081','AGAP006080','AGAP006079','AGAP006078','AGAP006077','AGAP006076','AGAP010489','AGAP006075','AGAP007286','AGAP007287','AGAP007289','AGAP010650','AGAP010648','AGAP009402','AGAP009065','AGAP005182','AGAP002191','AGAP001409','AGAP002190','AGAP000580','AGAP000641','AGAP000643','AGAP000642','AGAP000644','AGAP000640','AGAP000638','AGAP010649','AGAP011647','AGAP003306','AGAP012331','AGAP012325','AGAP012323','AGAP012321','AGAP012320','AGAP012319','AGAP012318','AGAP010409','AGAP008398','AGAP005208','AGAP003309','AGAP004433','AGAP003307','AGAP002189','AGAP002905','AGAP002188','AGAP002025','AGAP001189']
vb = get_obp_ids("../data/vb_agam_obp_geneIDs.txt")

#
# Get FASTA sequences of each gene's orthologues
#
print "gathering VB Agam OBP orthologous FASTA sequences..."
#gene_to_orthoFasta = get_orthologue_seqs(vb,source="web") 	# via VB web bot
gene_to_orthoFasta = get_orthologue_seqs(vb,source="../data/agam_orthologues") # via file

#
# Filter out non-anopheline species
#
print "Removing non-anopheline orthologues..."
print "\tUnfiltered FASTAs: ../data/agam_orthologues/*"
print "\tFiltered FASTAs: ../data/agam_orthologues_filtered/*"

for gene in vb:
	print "\t"+gene 
	filter_out_non_anopheline(gene)


#
# Check if Gene IDs appear >1 file
#
print "checking for overlapping genes between orthogroups..."

gene_to_orthoFasta = get_orthologue_seqs(vb,source="../data/agam_orthologues_filtered")

n_overlapping_orthogroups = test_gene_overlap(gene_to_orthoFasta)

#
# Tally species of each orthogroup
#

gene_to_orthoSpecies = {}

for gene in gene_to_orthoFasta.keys():

	print gene

	gene_to_orthoSpecies[gene] = []

	fasta_lines = gene_to_orthoFasta[gene].split("\n")

	for line in fasta_lines:

		if ">" in line:

			print "\t"+line

			gene_to_orthoSpecies[gene].append(line[15:])

	gene_to_orthoSpecies[gene] = Counter(gene_to_orthoSpecies[gene])

#
# Print non-1:1:1-... orthologues
#

for gene in gene_to_orthoSpecies.keys():

	if not sum(gene_to_orthoSpecies[gene].values())/len(gene_to_orthoSpecies[gene].values()) == 1:

		print gene
