###########
# Imports #
###########

from splinter import Browser
from splinter import exceptions
import pdb
import os

from time import sleep

#############
# Functions #
#############

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
		for gene in gene_to_orthoFasta.keys():
			fo = open("../data/agam_orthologues/"+gene+".txt","w")
			fo.write(gene_to_orthoFasta[gene])
			fo.close()
	else: 
		#
		# Grab fasta files from gene ids if available
		#
		# file names are given by: source + geneid + ".txt", read iteratively
		for gene in gene_ids:
			path = source + "/" + gene + ".txt"
			try:
				fi = open(path,"r")
			except IOError:
				raise IOError('The file: '+path+' does not exist, try: source="web" ...')
			lines = fi.readlines()
			fasta = "".join(lines)
			gene_to_orthoFasta[gene] = fasta
			fi.close()

	return gene_to_orthoFasta

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

#
# Get FASTA sequences of each gene's orthologues
#
print "gathering orthogroup fasta sequences..."
vb=['AGAP012658','AGAP012714','AGAP008278','AGAP008279','AGAP028120','AGAP008284','AGAP008282','AGAP008283','AGAP008281','AGAP008280','AGAP013182','AGAP003309','AGAP001189','AGAP001189','AGAP002025','AGAP002188','AGAP002905','AGAP002189','AGAP003307','AGAP003307','AGAP012319','AGAP004433','AGAP003306','AGAP005208','AGAP008398','AGAP010409','AGAP012318','AGAP012320','AGAP012321','AGAP012323','AGAP012323','AGAP012325','AGAP012331','AGAP001409','AGAP001409','AGAP000638','AGAP000638','AGAP000640','AGAP000640','AGAP000641','AGAP000642','AGAP000643','AGAP000644','AGAP000580','AGAP000580','AGAP002190','AGAP002190','AGAP010489','AGAP002191','AGAP002191','AGAP005182','AGAP005182','AGAP009065','AGAP009065','AGAP009402','AGAP009402','AGAP010648','AGAP010648','AGAP010650','AGAP010650','AGAP007289','AGAP009629','AGAP006079','AGAP006079','AGAP006081','AGAP011368','AGAP003530','AGAP001556','AGAP000279','AGAP000278','AGAP012867','AGAP006368','AGAP006278','AGAP002556','AGAP010649','AGAP010649','AGAP011647','AGAP011647','AGAP012322','AGAP012324']
#mo=['AGAP000278','AGAP000279','AGAP001556','AGAP012714','AGAP012867','AGAP006368','AGAP003530','AGAP013182','AGAP012658','AGAP012659','AGAP007283','AGAP006759','AGAP012324','AGAP012322','AGAP002556','AGAP007282','AGAP007281','AGAP009629','AGAP006760','AGAP006074','AGAP011368','AGAP011367','AGAP006081','AGAP006080','AGAP006079','AGAP006078','AGAP006077','AGAP006076','AGAP010489','AGAP006075','AGAP007286','AGAP007287','AGAP007289','AGAP010650','AGAP010648','AGAP009402','AGAP009065','AGAP005182','AGAP002191','AGAP001409','AGAP002190','AGAP000580','AGAP000641','AGAP000643','AGAP000642','AGAP000644','AGAP000640','AGAP000638','AGAP010649','AGAP011647','AGAP003306','AGAP012331','AGAP012325','AGAP012323','AGAP012321','AGAP012320','AGAP012319','AGAP012318','AGAP010409','AGAP008398','AGAP005208','AGAP003309','AGAP004433','AGAP003307','AGAP002189','AGAP002905','AGAP002188','AGAP002025','AGAP001189']

gene_to_orthoFasta = get_orthologue_seqs(vb,source="../data/agam_orthologues")

#
# Check if Gene IDs appear >1 file
#
print "checking for overlapping genes between orthogroups..."
gene_to_orthoGenes = {}

for geneId in gene_to_orthoFasta.keys():

	genes = []
	fi = open("../data/agam_orthologues/"+geneId+".txt","r")
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


for i in gene_to_orthoFasta.keys():
	if len(gene_to_orthoFasta[i]) < 1:
		print i
