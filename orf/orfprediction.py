import sys,os 
import logging

def run_orfpredictor(config):
	'''Run all steps of ORF prediction'''
	print "---"
	print "Starting ORF prediction"
	logging.INFO("Starting ORF prediction")
	create_directories = create_orf_directory(config)
	if create_directories == True:
		logging.INFO("Created ORF directory")
	else:
		logging.ERROR("Couldn't create ORF Directory")
		print "Couldn't create ORF Directory. Please check you have write-permissions"
		sys.exit(1)
	for organism in config["INPUT"]: # iterate over all organisms & run prediction for each
		run_blast(config,organism)
	
def create_orf_directory(config):
	'''Create ORF-Directory'''
	if os.path.isdir(config["OUTPUT"]["folder"]+"orfs/"):
		return True
	else:
		try:
			os.makedirs(config["OUTPUT"]["folder"]+"orfs/")
			return True
		except:
			return False

def run_blast(config,organism):
	'''Blast input-file against database'''
	request = "blastall -p blastx -m 8 -e 1e-5 -d " + 
			config["ORFPREDICTION"]["database"] +
			" -i " + organism["assembly_fasta"] + 
			" -o " + config["OUTPUT"]["folder"]+"orfs/"+organism["prefix"]+"-blast.csv"
	print request 
