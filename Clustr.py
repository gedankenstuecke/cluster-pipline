#!/usr/bin/python
from configobj import ConfigObj
import os,sys,getopt
import logging
from orf import orfpredictor
from cluster import orthomcl

def get_parameters(arguments):
	'''Check whether all parameters are present'''

	out_args = {}
	arguments.pop(0)
	optlist, arguments = getopt.getopt(arguments, 'hc:')

	for arg,opt in optlist:
		if arg == "-c":
			out_args["config_file"] = opt
		elif arg == "-h":
			print "Welcome to the Cluster-Generator"
			print "---"
			print "This pipeline does the following stuff:"
			print "1. Finds the ORFs in your contigs of multiple species"
			print "2. Clusters homolog contigs of the different species together"
			print "3. Reduces clusters to only keep ortholog sequences"
			print "4. Aligns AA those clusters which include 1 sequence of each species"
			print "5. Aligns the nucleotides according to the protein alignments while keeping the alignments codon-sensitive"
			print "---"
			print "Please specify the config file using './Clustr.py -c config.cfg'."
			print "If you don't provide a config-file the standard 'config.cfg' will be used"
			sys.exit(1)
		else:
			print "Please use '-h' to get information on how to use the pipeline"
			sys.exit(1)
	if out_args.has_key("config_file"):
		pass
	else:
		out_args["config_file"] = "config.cfg"
	
	return out_args

def read_config(arguments):
	'''Read Config'''
	config = ConfigObj(arguments["config_file"])
	if config["OUTPUT"]["folder"][-1] != "/":
		config["OUTPUT"]["folder"] = config["OUTPUT"]["folder"] + "/"
	if config["ORTHOMCL"]["location"][-1] != "/":
		config["ORTHOMCL"]["location"] = config["ORTHOMCL"]["location"] + "/"
	if config["CLUSTER"]["folder"][-1] != "/":
		config["CLUSTER"]["folder"] = config["CLUSTER"]["folder"] + "/"
	return config

def check_path(config):
	'''Check whether output-folder exists'''
	if os.path.isdir(config["OUTPUT"]["folder"]):
		return True
	else:
		try:
			os.makedirs(config["OUTPUT"]["folder"])
			return True
		except:
			return False

def set_logging(config):
	'''Set up logging-system for the pipeline'''
	if config["OUTPUT"]["logging"] == "True":
		logging.basicConfig(filename=config["OUTPUT"]["folder"]+"output.log",
				    filemode='w', 
				    format='%(levelname)s-%(asctime)s: %(message)s',
				    level=logging.DEBUG)
	else:
		logging.basicConfig(format='%(levelname)s-%(asctime)s: %(message)s',
				    level=logging.DEBUG)
	logging.info("Set up logging")

def main():
	'''Run ALL the stuff!'''
	print ""
	print "--------------------"
	print "Welcome to ClusterPy"
	print "--------------------"
	print ""
	parameters = get_parameters(sys.argv)
	config = read_config(parameters)	
	if check_path(config) != True:
		print "Could not create/open folder for output-data"
		print "Please make sure you have write-permissions for it"
		sys.exit(1)
	set_logging(config)	
	orfpredictor.run_orfpredictor(config)	
	orthomcl.runOrthoMCL(config)

main()
