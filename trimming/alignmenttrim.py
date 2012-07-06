from Bio import AlignIO
import sys
import os
import getopt
import subprocess
import logging

def read_alignment(infile):
    '''Read Alignment-Fasta'''
    alignment = AlignIO.read(open(infile),"fasta")
    return alignment


def find_longest_start_gap(sequence):
    '''Find longest Gap-sequence at sequence-start'''
    gap_counter = 1
    while gap_counter < len(sequence):
        if sequence.find("---"*gap_counter) == 0:
            gap_counter += 1
        else:
            break
    gap_counter = gap_counter-1
    return gap_counter*3


def find_longest_end_gap(sequence):
    '''Find longest Gap-sequence at sequence-end'''
    gap_counter = 1
    while gap_counter < len(sequence):
        if sequence[-3*gap_counter:] == "---"*gap_counter:
            gap_counter += 1
        else:
            break
    gap_counter = gap_counter-1
    return gap_counter*3


def row_iterator(alignment):
    '''Iterate over all rows'''
    start = 0
    end = 0
    for row in alignment:
        temp_start = find_longest_start_gap(str(row.seq))
        temp_end = find_longest_end_gap(str(row.seq))
        if temp_start > start:
            start = temp_start
        if temp_end > end:
            end = temp_end
    return [start,end]

def trim(config,infile):
    '''Trim sequences accordingly'''
    alignment = read_alignment(config["OUTPUT"]["folder"]+"cluster/nucleotide_alignments/"+infile)
    delimitors = row_iterator(alignment)
    if delimitors[1] != 0:
        trimmed_alignment = alignment[:,delimitors[0]:-delimitors[1]]
    else:
        trimmed_alignment = alignment[:,delimitors[0]:]
    removed_bases(alignment,delimitors)
    return trimmed_alignment

def removed_bases(alignment,delimitors):
    for row in alignment:
        print "Removed bases in front of "+row.id+": "+str(len(row.seq[:delimitors[0]])-row.seq[:delimitors[0]].count("-"))
        print "Removed bases in end of "+row.id+": "+str(len(row.seq[-delimitors[1]:])-row.seq[-delimitors[1]:].count("-"))

def remove_interim_gaps(alignment,cutoff):
    '''Remove gaps which are located inside of alignment'''
    counter = 0
    hit_counter = 0
    new_alignment = []
    gap_counter = 0
    while counter < len(alignment[0]):
        base_pairing = alignment[:,counter]
        if base_pairing.find("-") == -1:
            hit_counter += 1
            gap_counter = 0
        elif gap_counter < 3:
            gap_counter += 1
            hit_counter += 1
        else:
            #CHECK FOR %3 in LINE 89 AND 97
            if hit_counter >= cutoff:
                if new_alignment == []:
                    print "-- n=1 --"
                    print "HIT:"
                    print "Counter: "+str(counter)
                    print "Hit Counter. "+str(hit_counter)
                    print "Gap Counter: "+str(gap_counter)

                    overhead = alignment[:,counter-hit_counter:counter-gap_counter].get_alignment_length()%3
                    print alignment[:,counter-hit_counter+overhead:counter-gap_counter]
                    print "Overhead: "+ str(overhead)
                    new_alignment = alignment[:,counter-hit_counter+overhead:counter-gap_counter]
                else:
                    print "-- n>1 --"
                    print "HIT:"
                    print "Counter: "+str(counter)
                    print "Hit Counter. "+str(hit_counter)
                    print "Gap Counter: "+str(gap_counter)
                    overhead = alignment[:,counter-hit_counter:counter-gap_counter].get_alignment_length()%3
                    print "Overhead: "+ str(overhead)
                    new_alignment = new_alignment + alignment[:,counter-hit_counter+overhead:counter-gap_counter]
            hit_counter = 0
            gap_counter = 0
        counter += 1
    if hit_counter >= cutoff:
        if new_alignment == []:
            print "-- START & END --"
            print "HIT:"
            print "Counter: "+str(counter)
            print "Hit Counter. "+str(hit_counter)
            print "Gap Counter: "+str(gap_counter)
            overhead = alignment[:,counter-hit_counter:counter-gap_counter].get_alignment_length()%3
            print alignment[:,counter-hit_counter+overhead:counter-gap_counter]
            print "Overhead: "+ str(overhead)
            new_alignment = alignment[:,counter-hit_counter+overhead:counter-gap_counter]
        else:
            print "-- END --"
            print "HIT:"
            print "Counter: "+str(counter)
            print "Hit Counter. "+str(hit_counter)
            print "Gap Counter: "+str(gap_counter)
            overhead = alignment[:,counter-hit_counter:counter-gap_counter].get_alignment_length()%3
            print "Overhead: "+ str(overhead)
            new_alignment = new_alignment + alignment[:,counter-hit_counter+overhead:counter-gap_counter]
    hit_counter = 0
    
    if new_alignment == []:
        new_alignment = alignment
    
    return new_alignment

def trim_directory(config):
    '''create directory for trimmed output'''
    output_folder = config["OUTPUT"]["folder"]+"trimming/trimmed_alignments/"
    if os.path.isdir(output_folder):
        request = "rm -rf "+output_folder.replace("/trimmed_alignments/","")
        print request
        subprocess.call(request,shell=True)
        os.makedirs(output_folder)
    else:
        os.makedirs(output_folder)

def alignment_writer(alignment_to_write,config,alignment_name):
    '''write trimmed alignment'''
    output_folder = config["OUTPUT"]["folder"]+"trimming/trimmed_alignments/"
    if alignment_to_write != []:
        trim_length = config["TRIMMING"]["trim_length"]
        alignment_to_write = alignment_to_write[:,int(trim_length):-int(trim_length)]

        if alignment_to_write.get_alignment_length()%3 == 0:
            print "Length Alignment % 3: "+str(alignment_to_write.get_alignment_length()%3)
            AlignIO.write(alignment_to_write, output_folder+alignment_name, "fasta")
        else:
            print "ERROR: "+name

def get_file_list(config):
    files = os.listdir(config["OUTPUT"]["folder"]+"cluster/nucleotide_alignments/")
    return files

def trimNucleotideAlignments(config):
    '''Trim nucleotide-alignments according to config.cfg'''
    logging.info("Start trimming of nucleotide sequences")
    print "---"
    print "Start trimming of nucleotide sequences"
    trim_directory(config)
    alignments = get_file_list(config)
    for alignment in alignments:
        if os.path.getsize(config["OUTPUT"]["folder"] + "cluster/nucleotide_alignments/" + alignment) > 0:
            trimmed_alignment = trim(config,alignment)
            new_alignment = remove_interim_gaps(trimmed_alignment,config["TRIMMING"]["min_length_cutoff"])
            alignment_writer(new_alignment,config,alignment)
    print "---"
    print "Finished trimming"
    print "---"
    logging.info("Finished trimming of sequences")
