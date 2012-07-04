# Clustering pipeline
A small collection of scripts to perform clustering of transcriptome-sequences.

## Goals
* Predict open reading frames (ORFs) from transcriptome contigs using blast-hits as guidance
* Cluster orthologous ORFs of different species together
* Create peptide & nucleotide-alignments from those clusters

## Dependencies
* (Bio)python for the pipeline itself
* OrthoMCL v2.0, which needs 
** MySQL
** blastall-package
* Perl to run the protein2nucleotide-alignment-script

## Usage
In short:
* Run Clustr.py and wait

In detail:
* Clustr.py will try to find a file named "config.cfg" to get all the settings, an example is _config.cfg.example_ (surprise!). You can also specify the config-file using the _"-c"_ parameter while starting: `Clustr.py -c /some/config/somewhere.cfg`
* In the config-file you will have to provide (looking at the example should give you an idea): 
** The organisms/assemblies you want to cluster
** You can choose to skip the ORF prediction if you already have protein- & nucleotide-sequences for the ORFs
*** In this case: Specify the nucleotide-ORF-fasta as *assembly_fasta* and the protein-ORF-fasta as *peptides_fasta*
** If you want to perform the ORF prediction: Where the blastx-supported protein-blast-database is located you want to use for the ORF-prediction
** Where your makeblastdb-binary and the OrthoMCL-binaries are located
** The credentials for your MySQL-database (write & delete-access is required)
** Some output-details for OrthoMCL (name of the resulting clusters etc.)
** After setting up your config: Run Clustr.py

## To-Do
* Implement further trimming of resulting alignments
* Implement analysis of resulting alignments

## Known bugs
* Sometimes the nucleotide-sequences, given by the ORF-prediction, are longer than the corresponding protein-sequences. Due to this the protein2nucleotide-alignment won't work. Probably a bug in the orf-prediction-scripts. 
