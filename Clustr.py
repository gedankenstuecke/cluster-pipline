#!/usr/bin/python
from configobj import ConfigObj
import os,sys,getopt
import logging
import datetime
from orf import orfpredictor
from cluster import orthomcl
from filtering import filter
from trimming import alignmenttrim

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

def logStatistics(config):
    '''Print statistics to log-file'''
    logging.info("##############")
    logging.info("# STATISTICS #")
    logging.info("##############")
    logging.info("")
    logging.info("## ORFs ##")
    output_base = config["OUTPUT"]["folder"]
    for organism in config["INPUT"]:
        organism_hash = config["INPUT"][organism]
        if organism_hash["skip_orfs"] != "True":
            orf_file = output_base + "orfs/" + organism_hash["prefix"] + "-orfs.tsv"
        else:
            orf_file = output_base + "orfs/" + organism_hash["prefix"] + "-nucleotide-orfs.fasta"
        if organism_hash["skip_orfs"] == "True":
            orf_counter = returnNumberOfLines(orf_file,">")
            logging.info("Skipped prediction. ORFs for "+organism_hash["prefix"]+": "+str(orf_counter))
        else:
            logging.info("Did prediction: ORFs for "+organism_hash["prefix"]+": "+str(orf_counter))
    logging.info("")
    logging.info("## Clustering ##")
    logging.info("Number of total cluster: " + str(returnNumberOfLines(output_base + "cluster/groups.txt")))
    logging.info("Number of cluster with each species present: "+ str(returnNumberOfLines(output_base + "cluster/groups_filtered.txt")))
    logging.info("Number of cluster without paralogs: " + str(returnNumberOfLines(output_base + "cluster/paralog-free-clusters.csv")))
    logging.info("")
    logging.info("## Alignments ##")
    logging.info("Number of protein alignments: "+str(returnNumberOfFiles(output_base + "cluster/protein_alignments/")))
    logging.info("Number of nucleotide alignments: "+str(returnNumberOfFiles(output_base + "cluster/nucleotide_alignments/")))
    if config["TRIMMING"]["trim"] == "True":
        logging.info("")
        logging.info("## Trimming ##")
        logging.info("Number of trimmed alignments: " + str(returnNumberOfFiles(output_base + "trimming/trimmed_alignments/")))
    logging.info("##############")
    logging.info("# END OF LOG #")
    logging.info("##############")

def returnNumberOfLines(filename,filter_char=None):
    '''Count number of lines in file'''
    handle = open(filename,"r")
    counter = 0
    for line in handle:
        if filter_char is None:
            counter += 1
        elif line[0] == filter_char:
            counter += 1 
    return counter 

def returnNumberOfFiles(folder):
    print folder
    number_of_files = len([name for name in os.listdir(folder) if os.path.isfile(folder+name)]) 
    print number_of_files
    return number_of_files
    
def main():
    '''Run ALL the stuff!'''
    print ""
    print "--------------------"
    print "Welcome to ClusterPy"
    print "--------------------"
    print ""
    start = str(datetime.datetime.now())
    parameters = get_parameters(sys.argv)
    config = read_config(parameters)
    if check_path(config) != True:
        print "Could not create/open folder for output-data"
        print "Please make sure you have write-permissions for it"
        sys.exit(1)
    set_logging(config)
    orfpredictor.run_orfpredictor(config)
    orthomcl.runOrthoMCL(config)
    filter.runFiltering(config)
    if config["TRIMMING"]["trim"] == "True":
        alignmenttrim.trimNucleotideAlignments(config)
    logStatistics(config)
    stop = str(datetime.datetime.now())
    logging.info("Finished all steps")
    print "--"
    print "Job started at: "+start
    print "Job finished at: "+stop
    print "--"
    print ""
    print "----------------------"
    print "- Finished all steps -"
    print "----------------------"
    print ""
main()

