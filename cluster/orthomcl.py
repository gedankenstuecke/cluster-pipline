import sys,os
import subprocess
import logging

def clusterOrfs(config):
	print "foo"

def createOrthoDir(config):
	'''Create Directory for Cluster-Results'''
	if os.path.isdir(config["OUTPUT"]["folder"]+"cluster/"):
		return True
	else: 
		try:
			os.makedirs(config["OUTPUT"]["folder"]+"cluster/")
			return True
		except:
			return False


def createOrthoConfig(config):
	'''create sql-config for orthoMcl'''
	config_handle = open(config["OUTPUT"]["folder"]+"cluster/orthomcl.cfg","w")
	for entry,value in config["ORTHOMCL"].items():
		config_handle.write(entry+"="+value+"\n")
	config_handle.close()
	logging.info("Created OrthoMCL-config-file")

def setupDatabase(config):
	database_name = config["ORTHOMCL"]["dbConnectString"].split(":")[-1]
	user = config["ORTHOMCL"]["dbLogin"]
	password = config["ORTHOMCL"]["dbPassword"]
	request = "echo 'DROP DATABASE "+database_name+";' |mysql -u "+user+" --password="+password
	try: 
		return_value = subprocess.call(request, shell=True)
		if return_value == 1:
			logging.info("No orthomcl-database to drop")
			print "Don't worry, this happens if no previous database was available"
		elif return_value == 0:
			logging.info("Dropped orthomcl-database")
	except OSError, e:
		logging.error("Dropping database failed: "+e)
		print >>sys.stderr, "Execution failed: ",e
		sys.exit(1)

	request = "echo 'CREATE DATABASE "+database_name+";' |mysql -u "+user+" --password="+password

	try:
		return_value = subprocess.call(request, shell=True)
		if return_value == 1:
			logging.error("Error while creating database")
			sys.exit(1)
		elif return_value == 0:
			logging.info("Created database for orthomcl")
	except OSError, e:
		logging.error("Creating database failed: "+e)
		print >>sys.stderr, "Execution failed: ",e
		sys.exit(1) 

def runOrthoMCL(config):
	'''Run all steps of clustering'''
	print "---"
	print "Start clustering of ORFs"
	logging.info("Start clustering of ORFs")
	create_dir = createOrthoDir(config)
	if create_dir == True:
		logging.info("Created Cluster-Directory")
	else:
		logging.error("Couldn't create Directory")
		print "Couldn't create directory. Please check you have write-permissions"
		sys.exit(1)
	createOrthoConfig(config)
	setupDatabase(config)	
