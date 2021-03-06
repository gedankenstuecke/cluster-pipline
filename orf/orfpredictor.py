import sys,os
import logging
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
        if config["INPUT"][organism]["skip_orfs"] == "False":
            print "----"
            run_blast(config,config["INPUT"][organism])
            remove_duplicates(config,config["INPUT"][organism])
            predict_orfs(config,config["INPUT"][organism])
            extract_orfs(config,config["INPUT"][organism])
            if config["SNP"]["call_snps"] == "True":
                extractSnpOrfs(config,config["INPUT"][organism])
        else:
            print "----"
            movePeptides(config,config["INPUT"][organism])
            if config["SNP"]["call_snps"] == "True":
                moveSnps(config,config["INPUT"][organism])

def extractSnpOrfs(config,organism):
    '''Extract ORFs out of SNP-fastas'''
    orf_handle = open(config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"]+"-orfs.tsv")
    snp_handle = open(config["OUTPUT"]["folder"] + "snps/" + organism["prefix"]+"_snps.fasta")
    snp_contigs = SeqIO.to_dict(SeqIO.parse(snp_handle,"fasta"))
    snp_handle.close()
    orfs = orf_handle.readlines()
    orf_handle.close()
    sequences = {}
    for line in orfs:
        if line[0] != "#":
            line_array = line.split("\t")
            start = int(line_array[2])
            stop = int(line_array[3])
            print line_array[0]
            print "start: "+str(start)
            print "stop: "+str(stop)
            if start > stop:
                sequences[line_array[0]] = snp_contigs[line_array[0]].reverse_complement()[stop:start+1]
            else:
                sequences[line_array[0]] = snp_contigs[line_array[0]][start:stop+1]
    snp_orf_handle = open(config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"]+"-nt-snps.fasta","w")
    out_sequences = []
    for key,value in sequences.items():
        out_sequences.append(SeqRecord(Seq(str(value.seq)),id=key,description=""))
    SeqIO.write(out_sequences,snp_orf_handle,"fasta")
    snp_orf_handle.close()

def moveSnps(config,organism):
    '''Move Snps if already in ORF-format'''
    print "Moving SNP-Fasta for "+ organism["prefix"] +" to ORF-Location"
    request = "cp " + config["OUTPUT"]["folder"] + "snps/" + organism["prefix"] + "_snps.fasta "
    request = request + config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"] + "-nt-snps.fasta"
    try:
        return_value = subprocess.call(request,shell=True)
        if return_value != 0:
            print "Couldn't copy nucleotide-SNP-file for "+organism["prefix"]
            logging.error("Couldn't copy nucleotide-SNP-file for "+organism["prefix"])
            sys.exit(1)
        else:
            print "Copied SNPs-nucleotide file for "+organism["prefix"]
            logging.info("Copied SNPs-nucleotide file for "+organism["prefix"])
    except OSError, e:
        print "Couldn't copy nucleotide-SNP-file for "+organism["prefix"]
        logging.error("Couldn't copy nucleotide-SNP-file for "+organism["prefix"])
        sys.exit(1)

def movePeptides(config,organism):
    '''Move peptide to fit with structure of pipeline'''
    print "Copy files for "+ organism["prefix"]
    request = "cp " + organism["assembly_fasta"] + " " + config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"] + "-nucleotide-orfs.fasta"
    try:
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't copy nucleotide-file for "+organism["prefix"])
            print "Couldn't copy nucleotide-file for "+organism["prefix"]
            sys.exit(1)
        else:
            logging.info("Copied nucleotide-file for "+organism["prefix"])
    except OSError, e:
        logging.error("Couldn't copy nucleotide-file for "+organism["prefix"])
        print "Couldn't copy nucleotide-file for "+organism["prefix"]
        sys.exit(1)

    request = "cp " + organism["peptides_fasta"] + " " + config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"] + "-aa-orfs.fasta"
    try:
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't copy aa-file for "+organism["prefix"])
            print "Couldn't copy aa-file for "+organism["prefix"]
            sys.exit(1)
        else:
            logging.info("Copied aa-file for "+organism["prefix"])
    except OSError, e:
        logging.error("Couldn't copy aa-file for "+organism["prefix"])
        print "Couldn't copy aa-file for "+organism["prefix"]
        sys.exit(1)


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

def extract_orfs(config,organism):
    '''Extract nucleotide & AA-sequences out of ORF-results'''
    print "Extracting sequences"
    logging.info("Extracting nt & aa sequences out of ORF-results")
    orf_handle = open(config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"]+"-orfs.tsv")
    nucleotide_handle = open(config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"]+"-nucleotide-orfs.fasta","w")
    aa_handle = open(config["OUTPUT"]["folder"] + "orfs/" + organism["prefix"]+"-aa-orfs.fasta","w")

    for line in orf_handle:
        if line[0] != "#":
            line_array = line.strip().split("\t")
            identifier = line_array[0]
            nt_sequence = line_array[4]
            aa_sequence = line_array[5]
            nucleotide_handle.write(">"+identifier+"\n")
            nucleotide_handle.write(nt_sequence+"\n")
            aa_handle.write(">"+identifier+"\n")
            aa_handle.write(aa_sequence+"\n")
    orf_handle.close()
    aa_handle.close()
    nucleotide_handle.close()
    print "Extracted sequences"
    logging.info("Extracted sequences for "+organism["prefix"])
