import sys,os 
import logging
import subprocess

def run_orfpredictor(config):
	'''Run all steps of ORF prediction'''
	print "---"
	print "Starting ORF prediction"
	logging.info("Starting ORF prediction")
	create_directories = create_orf_directory(config)
	if create_directories == True:
		logging.info("Created ORF directory")
	else:
		logging.error("Couldn't create ORF Directory")
		print "Couldn't create ORF Directory. Please check you have write-permissions"
		sys.exit(1)
	for organism in config["INPUT"]: # iterate over all organisms & run prediction for each
		print "----"
		run_blast(config,config["INPUT"][organism])
		remove_duplicates(config,config["INPUT"][organism])
		predict_orfs(config,config["INPUT"][organism])
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
	print "Blast for "+organism["prefix"]+" is running"
	request = "blastall -p blastx -m 8 -e 1e-5 -d " 
	request = request + config["ORFPREDICTION"]["database"]
	request = request + " -i " + organism["assembly_fasta"]
	request = request + " -o " + config["OUTPUT"]["folder"]+"orfs/"+organism["prefix"]+"-blast.csv"
	try:
		return_value = subprocess.call(request, shell=True)
		if return_value < 0:
			logging.error("Blast for "+organism["prefix"]+
				" was terminated with signal "+str(return_value))
			print "Blast for "+organism["prefix"]+" was terminated: ",str(return_value)
			print >>sys.stderr
			sys.exit(1)

		elif return_value > 0:
			logging.error("Blast for "+organism["prefix"]+
				" failed with return value "+str(return_value))
			print "Blast for "+organism["prefix"]+ " returned: ", str(return_value)
			print >>sys.stderr
			sys.exit(1)

		elif return_value == 0:
			logging.info("Blast for "+organism["prefix"]+" finished successfully")
			print "Blast for "+organism["prefix"]+ " finished successfully"
	except OSError, e:
		logging.error("Blast for "+organism["prefix"]+" failed: "+e)
		print  >>sys.stderr, "Execution failed: ", e
		sys.exit(1)

def remove_duplicates(config,organism):
	'''Remove duplicate hits out of blast-results'''
	blast_handle = open(config["OUTPUT"]["folder"]+"orfs/"+organism["prefix"]+"-blast.csv","r")
	out_handle = open(config["OUTPUT"]["folder"]+"orfs/"+organism["prefix"]+"-filteredblast.csv","w")
	queries = []
	for line in blast_handle:
		line_array = line.strip().split("\t")
		if line_array[0] in queries:
			next
		else:
			out_handle.write(line)
			queries.append(line_array[0])
	blast_handle.close()
	out_handle.close()
	logging.info("Made blast-results for "+organism["prefix"]+ " unique.")

def predict_orfs(config,organism):
	'''Predict ORF-sequences through ORFPredictor'''
	print "Gathering ORFs for "+organism["prefix"]
	logging.info("Starting to gather ORFs for "+organism["prefix"])
	request = "python " + os.getcwd()
	request = request +"/orf/orfpredictor/ORFPREDICTORRR.py -i "
	request = request + config["OUTPUT"]["folder"] + "orfs/"
	request = request + organism["prefix"]+ "-filteredblast.csv -j " 
	request = request + organism["assembly_fasta"]
	request = request + " -t " + config["ORFPREDICTION"]["minimum_length"]
	request = request + " > " + config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"]
	request = request + "-orfs.tsv" 
	return_value = subprocess.call(request, shell=True)
	if return_value == 0:
		logging.info("Gathered ORFs for "+organism["prefix"])
		print "Got ORFs for "+organism["prefix"]
	else:
		logging.error("Gathering ORFs for "+organism["prefix"]+" ended with error")
		print "Encountered an error while getting ORFs for "+organism["prefix"]
		sys.exit(1)
