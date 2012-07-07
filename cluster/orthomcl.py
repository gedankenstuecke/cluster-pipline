import sys,os
import subprocess
import logging

def createOrthoDir(config):
    '''Create Directory for Cluster-Results'''
    if os.path.isdir(config["OUTPUT"]["folder"]+"cluster/"):
        return True
    else:
        try:
            os.makedirs(config["OUTPUT"]["folder"]+"cluster/")
            os.makedirs(config["OUTPUT"]["folder"]+"cluster/compliantFasta/")
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
    '''create database & install db-schema for orthomcl'''
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

    request = config["ORTHOMCL"]["location"]+"orthomclInstallSchema "+config["OUTPUT"]["folder"]+"cluster/orthomcl.cfg "
    request = request + config["OUTPUT"]["folder"]+"cluster/schema_install.log"
    try:
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't install db-schema")
            print "Couldn't install db-schema"
            sys.exit(1)
        else:
            logging.info("Installed DB-schema")
    except OSError, e:
        logging.error("Couldn't install db-schema")
        print "Couldn't install db-schema"
        sys.exit(1)

def adjustFasta(config):
    '''create orthomcl-compliant fasta-files'''
    print "Adjusting FASTAs for OrthoMCL"
    for organism_name in config["INPUT"]:
        organism = config["INPUT"][organism_name]
        logging.info("Adjusting fasta for "+organism["prefix"])
        request = config["ORTHOMCL"]["location"]+"orthomclAdjustFasta "
        request = request + organism["prefix"] + " "
        request = request + config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"]+"-aa-orfs.fasta "
        request = request + organism["unique_field"]
        try:
            return_value = subprocess.call(request, shell=True)
            if return_value != 0:
                logging.error("Some error occured while adjusting fasta for "+organism["prefix"])
                print "Some error occured during adjusting fasta-files for OrthoMCL"
                sys.exit(1)
            else:
                logging.info("Adjusted fasta for "+organism["prefix"])
                request = "mv "+organism["prefix"]+".fasta "
                request = request + config["OUTPUT"]["folder"] + "cluster/compliantFasta/"
                subprocess.call(request, shell=True)
        except OSError, e:
            logging.error("Error while adjusting FAST for "+organism["prefix"])
            print >>sys.stderr, "Execution failed: ",e
            sys.exit(1)
    print "Adjusted FASTAs"

def filterFasta(config):
    '''Filter Files for further use in OrthoMCL'''
    logging.info("Start filtering compliant fastas for OrthoMCL")
    print "Start filtering FASTAs for OrthoMCL"
    request = config["ORTHOMCL"]["location"]+"orthomclFilterFasta "
    request = request + config["OUTPUT"]["folder"] + "cluster/compliantFasta/ "
    request = request + config["ORTHOMCL"]["minimum_length"] + " "
    request = request + config["ORTHOMCL"]["max_percent_stop"]
    try:
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Some error occured during OrthoMCL filtering")
            print "Some error occured during OrthoMCL filtering"
            sys.exit(1)
        else:
            logging.info("Filtered FASTAs for OrthoMCL")
            request = "mv goodProteins.fasta "
            request = request + config["OUTPUT"]["folder"] + "cluster/"
            subprocess.call(request, shell=True)
            request = "mv poorProteins.fasta "
            request = request + config["OUTPUT"]["folder"] + "cluster/"
            subprocess.call(request, shell=True)
    except OSError, e:
        logging.error("Some error occured during OrthoMCL filtering")
        print "Some error occured during OrthoMCL filtering"
        sys.exit(1)

def createBlastDatabase(config):
    '''Create Blast-DB out of goodProteins.fasta'''
    logging.info("Create blastdb out of goodProteins.fasta")
    request = config["ORTHOMCL"]["makeblastdb"] + " -dbtype prot -in "
    request = request + config["OUTPUT"]["folder"] + "cluster/goodProteins.fasta"
    try:
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't create blastdb out of goodProteins.fasta")
            print "Error while creating blast-database for OrthoMCL"
            sys.exit(1)
        else:
            logging.info("Created blastdb out of goodProteins.fasta")
    except OSError, e:
        logging.error("Couldn't create blastdb out of goodProteins.fasta")
        print "Error while creating blast-database for OrthoMCL"
        sys.exit(1)

def runOrthoBlast(config):
    '''Run All-Versus-All Blast for OrthoMCL'''
    request = "blastall -i "
    request = request + config["OUTPUT"]["folder"] + "cluster/goodProteins.fasta -m 8 -d "
    request = request + config["OUTPUT"]["folder"] + "cluster/goodProteins.fasta -F 'm S' "
    request = request + "-v 100000 -b 100000 -e 1e-5 -p blastp -o "
    request = request + config["OUTPUT"]["folder"] + "cluster/all_vs_all_blast.csv"
    try:
        logging.info("Starting all-vs-all-blast")
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't perform all-vs-all blast")
            print "Couldn't perform all-vs-all blast"
            sys.exit(1)
        else:
            logging.info("Performed all-vs-all blast")
    except OSError, e:
        logging.error("Couldn't perform all-vs-all blast")
        print "Couldn't perform all-vs-all blast"
        sys.exit(1)

def BlastParser(config):
    '''Re-format blast-output for OrthoMCL'''
    request = config["ORTHOMCL"]["location"]+"orthomclBlastParser "
    request = request + config["OUTPUT"]["folder"] + "cluster/all_vs_all_blast.csv "
    request = request + config["OUTPUT"]["folder"]+"cluster/compliantFasta"
    request = request + "> "+config["OUTPUT"]["folder"]+"cluster/similarSequences.txt"
    try:
        logging.info("Starting reformatting of blast-results")
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't perform reformatting")
            print "Couldn't perform reformatting"
            sys.exit(1)
        else:
            logging.info("Performed reformatting")
    except OSError, e:
        logging.error("Couldn't perform reformatting")
        print "Couldn't perform reformatting"
        sys.exit(1)

def LoadBlast(config):
    '''Load data into the orthomcl-database'''
    request = config["ORTHOMCL"]["location"] + "orthomclLoadBlast "
    request = request + config["OUTPUT"]["folder"]+"cluster/orthomcl.cfg "
    request = request + config["OUTPUT"]["folder"]+"cluster/similarSequences.txt"
    try:
        logging.info("Loading reformatted blast-results into database")
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't load data into database")
            print "Couldn't load data into database"
            sys.exit(1)
        else:
            logging.info("Loaded data into database")
    except OSError, e:
        logging.error("Couldn't load data into database")
        print "Couldn't load data into database"
        sys.exit(1)

def createPairs(config):
    '''Create pairings in database'''
    request = config["ORTHOMCL"]["location"] + "orthomclPairs "
    request = request + config["OUTPUT"]["folder"]+"cluster/orthomcl.cfg "
    request = request +  config["OUTPUT"]["folder"]+"cluster/ortho_pairs.log cleanup=no"
    try:
        logging.info("Pairing results")
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't pair data")
            print "Couldn't pair data"
            sys.exit(1)
        else:
            logging.info("Paired data")
    except OSError, e:
        logging.error("Couldn't pair data")
        print "Couldn't pair data"
        sys.exit(1)

def dumpPairs(config):
    '''Dump pairs out of database'''
    request = config["ORTHOMCL"]["location"] + "orthomclDumpPairsFiles "
    request = request + config["OUTPUT"]["folder"]+"cluster/orthomcl.cfg"
    try:
        logging.info("Pairing results")
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't dump data")
            print "Couldn't dump data"
            sys.exit(1)
        else:
            logging.info("Dumped data")
            request = "mv mclInput " + config["OUTPUT"]["folder"]+"cluster/"
            subprocess.call(request, shell=True)
            if os.path.isdir(config["OUTPUT"]["folder"]+"cluster/pairs"):
                request = "rm -r "+config["OUTPUT"]["folder"]+"cluster/pairs"
                subprocess.call(request, shell=True)
            request = "mv pairs/ "+config["OUTPUT"]["folder"]+"cluster/pairs"
            subprocess.call(request, shell=True)
    except OSError, e:
        logging.error("Couldn't dump data")
        print "Couldn't dump data"
        sys.exit(1)

def runMcl(config):
    '''run mcl'''
    request = "mcl "+config["OUTPUT"]["folder"]+"cluster/mclInput --abc -I 1.5 -o " + config["OUTPUT"]["folder"]+"cluster/mclOutput"
    try:
        logging.info("Running mcl")
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Encountered error during mcl-run")
            print "Encountered error during mcl-run"
            sys.exit(1)
        else:
            logging.info("Ran mcl")
    except OSError, e:
        logging.error("Encountered error during mcl-run")
        print "Encountered error during mcl-run"
        sys.exit(1)

def mclToGroups(config):
    '''turn mcl-results into groups'''
    request = config["ORTHOMCL"]["location"]+"orthomclMclToGroups "+ config["ORTHOMCL"]["result_prefix"] + " "
    request = request + config["ORTHOMCL"]["start_number"] + " < " + config["OUTPUT"]["folder"]+"cluster/mclOutput > "
    request = request + config["OUTPUT"]["folder"]+"cluster/groups.txt"
    try:
        logging.info("Running mcl to Groups")
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Failure during mcl-to-groups")
            print "An error occured during creation of the groups"
            sys.exit(1)
        else:
            logging.info("Exported final orthomcl-clusters")

    except OSError, e:
        logging.error("Failure during mcl-to-groups")
        print "An error occured during creation of the groups"
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
    adjustFasta(config)
    filterFasta(config)
    createBlastDatabase(config)
    runOrthoBlast(config)
    BlastParser(config)
    LoadBlast(config)
    createPairs(config)
    dumpPairs(config)
    runMcl(config)
    mclToGroups(config)
    print "---"
    print "Finished all clustering steps"
    print "---"
    logging.info("Finished all clustering steps")
