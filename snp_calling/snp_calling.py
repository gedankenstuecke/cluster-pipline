import subprocess
import pysam
import os,sys
import logging
from Bio import SeqIO
from string import upper
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

codes = {
        "AC":"M",
        "AG":"R",
        "AT":"W",
        "CG":"S",
        "CT":"Y",
        "GT":"K",
        "ACG":"V",
        "ACT":"H",
        "AGT":"D",
        "CGT":"B",
        "ACGT":"N",
        "N":"N",
}

def readSam(filename):
    sam = pysam.Samfile(filename,"rb")
    return sam

def readContigs(filename):
    handle = open(filename,"r")
    contigs = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return contigs

def iterateContigs(contigs,samfile,minimum_coverage=10,minimum_frequency=0.2):
    out_contigs = {}
    contigs_with_snps = 0
    number_of_snps = 0
    for contig_name,contig in contigs.items():
        found_snp = False
        sequence_list = list(str(contig.seq))
        iterator = samfile.pileup(contig_name)
        for pileupcolumn in iterator:
            position = pileupcolumn.pos
            number_of_reads = pileupcolumn.n
            if number_of_reads > minimum_coverage:
                genotypes = callSnp(pileupcolumn,minimum_frequency,number_of_reads)
                if len(genotypes) > 1:
                    genotype_list = list(genotypes)
                    genotype_list.sort()
                    replacement_base = codes["".join(genotype_list)]
                    sequence_list[position] = replacement_base
                    found_snp = True
                    number_of_snps += 1
        out_contigs[contig_name] = "".join(sequence_list)
        if found_snp == True:
            contigs_with_snps += 1
    print "Number of SNPs: " + str(number_of_snps)
    print "Number of Contigs with SNPs: "+ str(contigs_with_snps)
    return out_contigs

def callSnp(pileupcolumn,minimum_frequency,number_of_reads):
    bases = {"A":0,"C":0,"G":0,"T":0}
    out_bases = ""
    for read in pileupcolumn.pileups:
        if "ACGT".find(read.alignment.seq[read.qpos]) != -1:
            bases[upper(read.alignment.seq[read.qpos])] += 1

    for base in bases:
        if float(bases[base])/number_of_reads > minimum_frequency:
            out_bases = out_bases + base
    if out_bases == "":
        out_bases = "N"
    return out_bases

def dumpSnpSequences(out_contigs,out_file):
    handle = open(out_file,"w")
    sequences = []
    for key,value in out_contigs.items():
        sequences.append(SeqRecord(Seq(value),id=key,description=""))
    SeqIO.write(sequences,handle,"fasta")
    handle.close()

def createSnpDirectory(config):
    '''Create SNP-Directory'''
    if os.path.isdir(config["OUTPUT"]["folder"]+"snps/"):
        return True
    else:
        try:
            os.makedirs(config["OUTPUT"]["folder"]+"snps/")
            return True
        except:
            return False

def moveSnpFile(config,organism):
    '''Move already done SNP file to fit structure of pipeline'''
    print "Copy SNP file for "+organism["prefix"]
    request = "cp " + organism["snp_file"] + " " + config["OUTPUT"]["folder"] + "snps/" + organism["prefix"]+"_snps.fasta"
    try:
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't copy nucleotide-file for "+organism["prefix"])
            print "Couldn't copy nucleotide-file for "+organism["prefix"]
            sys.exit(1)
        else:
            logging.info("Copied SNP-fasta for "+organism["prefix"])
    except OSError, e:
        logging.error("Couldn't copy nucleotide-file for "+organism["prefix"])
        print "Couldn't copy nucleotide-file for "+organism["prefix"]
        sys.exit(1)

def runSnpCalling(config):
    '''
    Predict SNPs for all specified pairs of contigs/bam-files
    Requires contigs used as reference sequence and an indexed
    and sorted BAM file. To create one use samtools
    '''
    logging.info("Start calling SNPs")
    create_directory = createSnpDirectory(config)
    if create_directory == True:
        logging.info("Created SNP directory")
    else:
        logging.error("Couldn't create SNP directory")
        print "Couldn't create SNP directory. Please check that you've write-permissions"
        sys.exit(1)
    for organism in config["INPUT"]:
        if config["INPUT"][organism]["skip_snps"] != "True":
            try:
                out_file = config["OUTPUT"]["folder"]+"snps/"
                out_file = out_file + config["INPUT"][organism]["prefix"]+"_snps.fasta"
                print "----"
                print "call SNPs for "+config["INPUT"][organism]["prefix"]
                contigs = readContigs(config["INPUT"][organism]["assembly_fasta"])
                sam = readSam(config["INPUT"][organism]["bam_file"])
                out_contigs = iterateContigs(contigs=contigs,
                                            samfile=sam,
                                            minimum_coverage=int(config["SNP"]["coverage"]),
                                            minimum_frequency=float(config["SNP"]["frequency"]))
                dumpSnpSequences(out_contigs=out_contigs,
                                            out_file=out_file)
            except:
                logging.error("SNP Calling: An error occured. Have a look at the output")
                print "Some error occured during SNP-parsing."
                sys.exc_info()[0] 
                sys.exit(1)
        else:
            moveSnpFile(config,config["INPUT"][organism])            
