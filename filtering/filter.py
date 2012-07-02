import sys,os
import logging
import subprocess

def checkNumber(config):
    '''Check that all species are represented in cluster'''
    cluster_handle = open(config["OUTPUT"]["folder"]+"cluster/groups.txt","r")
    cluster_hash = {}
    counter = 0
    for line in cluster_handle:
        temp_hash = {}
        line_array = line.strip().split(" ")
        for organism in config["INPUT"]:
            prefix = config["INPUT"][organism]["prefix"]
            if line.find(prefix) != -1:
                temp_hash[prefix] = "there"
        if len(temp_hash) == len(config["INPUT"]):
            cluster_hash[line_array[0][:-1]] = line_array[1:]
        else:
            counter += 1
    cluster_handle.close()
    cluster_output = open(config["OUTPUT"]["folder"]+"cluster/groups_filtered.txt","w")
    for key,values in cluster_hash.items():
        output = key + "("+str(len(values))+" genes,"+str(len(config["INPUT"]))+" taxa):\t"
        for value in values:
            output = output + value+"("+value.split("|")[0]+") "
        output.strip()
        output = output + "\n"
        cluster_output.write(output)
    cluster_handle.close()

def createOrthoGG(config):
    '''Create Ortho-GG-File to make data compatible with scripts for old orthomcl-version'''
    ortho_file = open(config["OUTPUT"]["folder"]+"cluster/ortho.gg","w")
    for organism in config["INPUT"]:
        organism_file = open(config["OUTPUT"]["folder"]+"cluster/compliantFasta/"+config["INPUT"][organism]["prefix"]+".fasta")
        out_line = config["INPUT"][organism]["prefix"]+": "
        for line in organism_file:
            if line[0] == ">":
                out_line = out_line + line[1:-1] + " "
        out_line = out_line[:-1]+ "\n"
        ortho_file.write(out_line)
        organism_file.close()
    ortho_file.close()

def createClusters(config):
    '''Create cluster-fasta-files'''
    request = "python " + config["CLUSTER"]["folder"] + "06-build-clusters.py "
    request = request + config["OUTPUT"]["folder"] + "cluster/groups_filtered.txt "
    request = request + config["OUTPUT"]["folder"] + "cluster/goodProteins.fasta"
    try:
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            logging.error("Couldn't create FASTAs for clusters")
            print "Couldn't create FASTAs for clusters"
            sys.exit(1)
        else:
            logging.info("Created FASTAs out of clusters")
            print "Created FASTAs out of clusters"
    except OSError, e:
        logging.error("Couldn't create FASTAs for clusters")
        print "Couldn't create FASTAs for clusters"
        sys.exit(1)

def moveParalogCluster(config):
    '''Move paralog-results to result-folder'''
    if os.path.isdir(config["OUTPUT"]["folder"]+"cluster/paralog_clusters/"):
        try:
            request = "rm "+config["OUTPUT"]["folder"]+"cluster/paralog_clusters/*"
            subprocess.call(request,shell=True)
        except:
            logging.error("Couldn't move clusters")
            print "Couldn't move clusters"
            sys.exit(1)
    else:
        try:
            os.makedirs(config["OUTPUT"]["folder"]+"cluster/paralog_clusters/")
        except:
            print "Couldn't move clusters"
            logging.error("Couldn't move clusters")
            sys.exit(1)
    try:
        request = "mv "+config["ORTHOMCL"]["result_prefix"]+"*.*fasta "
        request = request + config["OUTPUT"]["folder"]+"cluster/paralog_clusters"
        subprocess.call(request,shell=True)        
        request = "mv "+config["ORTHOMCL"]["result_prefix"]+"*.ids "
        request = request + config["OUTPUT"]["folder"]+"cluster/paralog_clusters"
        subprocess.call(request,shell=True)
    except:    
            print "Couldn't move clusters"
            logging.error("Couldn't move clusters")
            sys.exit(1)

def filterClusters(config):
    '''Filter paralog sequences out of clusters'''
    request = "python "+ config["CLUSTER"]["folder"] + "07-remove-paralogs.py "
    request = request + config["OUTPUT"]["folder"] + "cluster/paralog_clusters/ "
    request = request + config["OUTPUT"]["folder"] + "cluster/goodProteins.fasta "
    request = request + config["OUTPUT"]["folder"] + "cluster/ortho.gg"
    try:
        return_value = subprocess.call(request, shell=True)
        if return_value != 0:
            print "Couldn't remove paralogs"
            logging.error("Couldn't remove paralogs")
            sys.exit(1)
        else:
            request = "mv paralog-free-clusters/ "+ config["OUTPUT"]["folder"]+"cluster/"
            subprocess.call(request, shell=True)
            request = "mv noparalogs.orthomcl.out " + config["OUTPUT"]["folder"] +"cluster/paralog-free-clusters.csv"
            subprocess.call(request, shell=True)
            logging.info("Removed paralogs")
            print "Removed paralogs"
    except OSError, e:
        print "Couldn't remove paralogs"
        logging.error("Couldn't remove paralogs")
        sys.exit(1)

def alignProteins(config):
    if os.path.isdir(config["OUTPUT"]["folder"] + "cluster/protein_alignments"):
        try:
            request = "rm " + config["OUTPUT"]["folder"]+"cluster/protein_alignments/*"
            subprocess.call(request,shell=True)
        except:
            logging.error("Couldn't align proteins")
            print "Couldn't align proteins"
            sys.exit(1)
    else:
        try:
            os.makedirs(config["OUTPUT"]["folder"]+"cluster/protein_alignments/")
        except:
            logging.error("Couldn't align proteins")
            print "Couldn't align proteins"
            sys.exit(1)

    files = os.listdir(config["OUTPUT"]["folder"]+"cluster/paralog-free-clusters/")
    fastas = []
    for f in files:
        if f.find(".ufasta") != -1:
            fastas.append(f)

    for fasta in fastas:
        request = os.getcwd() + "/filtering/muscle3.8.31 -in "
        request = request + config["OUTPUT"]["folder"]+"cluster/paralog-free-clusters/" 
        request = request + fasta + " -out "
        request = request + config["OUTPUT"]["folder"]+"cluster/protein_alignments/"
        request = request + fasta.replace(".ufasta","protein_alignment.fasta")
        subprocess.call(request, shell=True)
        

def runFiltering(config):
    '''Run all steps of filtering & creating final clusters'''
    checkNumber(config)
    createOrthoGG(config)
    createClusters(config)
    moveParalogCluster(config)
    filterClusters(config)
    alignProteins(config)
    
def createFakeConfig():
    config = {}
    config["OUTPUT"] = {}
    config["OUTPUT"]["folder"] = "/bastian/hanno/orthomcl_rip_pig/"
    config["INPUT"] = {}

    config["INPUT"]["ORGANISM1"] = {}
    config["INPUT"]["ORGANISM2"] = {}
    config["INPUT"]["ORGANISM1"]["prefix"] = "pig"
    config["INPUT"]["ORGANISM2"]["prefix"] = "rip"
    config["CLUSTER"]["folder"] = "/bastian/cluster_pipeline/filtering/"
    return config

# config = createFakeConfig()
# checkNumber(config)
# createOrthoGG(config)
