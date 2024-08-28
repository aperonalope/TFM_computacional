"""
Created on April 2024
Identification of circular RNAs from the reads of a SAM file and an annotation file in GTF format
@author: Alvaro Perona (aperonalope@alumni.unav.es)
"""

import argparse
import os
from utils import _step_
import pandas as pd
import numpy as np
import re

def main(input_csv: str, input_gtf, suffix: str):
    # Obtain full directory and file name without extensions (.fastq.gz, .fastq, .fq.gz, .fq....) 
    save_dir, file_root = _step_(input_csv, "")

    # Open the input csv file, the input annotation file, and the output csv file
    inputfile = open(input_csv, "r")
    input2 = pd.read_csv(input_gtf, header=None, sep=" ")
    outputfile = open(os.path.join(save_dir, file_root + "_resumen.csv"), "wb")

    # Add the header to the output file
    strings = "read_length counts read_length chr_position chr gene CIGARS\n"
    outputfile.write(strings.encode())

    # Separate the columns of the input annotation file into numpy vectors
    chr_gtf = input2.iloc[:, 0].to_numpy()
    ini_gtf = input2.iloc[:, 1].to_numpy()
    fin_gtf = input2.iloc[:, 2].to_numpy()
    gen_gtf = input2.iloc[:, 4].to_numpy()

    # We read the first line of the file
    l1 = inputfile.readline().rstrip("\n")
    
    #We initiate the value of position
    position = -1000
    
    # We start an infinite loop
    while True:
        # We initialize the cigar list and a counter that will turn one at the last perfect match of a locus.
        soft_pe = 0
        cigar = []
        
        # This condition is to check if the file is empty, and if so, we break the loop. This ensures that the infitite while truee loop breaks
        if l1 == "":
            break
        

        # We start a second loop
        while True:
            
            
            
            # We save the l1 line in l2 to be able to compare two consecutive lines
            l2 = l1

            #We extract information of the read into variables.
            # The read name contains the rotation, the read name, and the count number, separated by "_"
            read_name = l2.split(" ")[0] + "_" + l2.split(" ")[1] + "_" + l2.split(" ")[2]

            # The position is the position of the read in the chromosome
            position = l2.split(" ")[3]

            # We store the cigar string of the current line in the cigar list
            cigar.append(l2.split(" ")[4])

            # We read the next line
            l1 = inputfile.readline().rstrip("\n")

            # If the line does not have the expected number of fields, break the loop
            if len(l1.split(" ")) != 7:
                break
            
            
            # This condition is to break the appending of cigars, as we consider that if the conditions are met, the new line does not correspond to the same circular RNA locus
            #The conditions are: change of more than 10 nucleotides between two consecutive lines or that the read name(without rotations) changes.
            if (abs(int(l1.split(" ")[3]) - int(position))) >= 10 or l1.split(" ")[1] != read_name.split("_")[1]:
                
                soft_pe += 1
                #We store the information of this read in new variables (probably not needed, but too afraid to change it for the time being)
                rotation_perfect = read_name.split("_")[0]
                name_perfect = read_name.split("_")[1]
                count = read_name.split("_")[2]
                cigar_bueno = l2.split(" ")[4]
                sequence_perfe = l2.split(" ")[6]
                chromosome = l2.split(" ")[5]
                
                
                
                #We initiallize the list of anotated genes as NA so that if none is found, NA is written in the output csv
                
                anotated_genes = "NA"
                
                
                # If the variable is 1, we store the data of the variables in a record in the output csv and we break the innermost loop
                if soft_pe >= 1:
                    genes = []
                    #We parse the anotation file vectors until it finds an index where the postion of our reed is bigger than the initial position of a gene and smaller than the final position.

                    for k in np.arange(0, len(input2)):
                        if chromosome == chr_gtf[k]:
                            if int(position) >= ini_gtf[k] and int(position) <= fin_gtf[k]:
                                genes.append(gen_gtf[k])

                    if len(genes) == 0:
                        anotated_genes = "NA"
                    elif len(genes) == 1:
                        anotated_genes = genes[0]
                    else:
                    #If more than one gene falls in those coordinates we store all of them separated by a "-"
                        anotated_genes = "-".join(genes)
                    strings = f"{rotation_perfect}{name_perfect} {count} {len(sequence_perfe)} {position} {chromosome} {anotated_genes} {sequence_perfe} {cigar_bueno}\n"
                    outputfile.write(strings.encode())
                break
    
    # We close the input and output files
    inputfile.close()
    outputfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Identifying circulars')
    parser.add_argument('-i', '--READS', type=str, required=True, help='Select the path where the input modified SAM file is located')
    parser.add_argument('-a', '--ANNOTATION', type=str, required=True, help='Select the path where the modified anotation file is located')
    parser.add_argument('-s', '--SUFFIX', type=str, required=False, default="", help='Select the suffix for saving output files with')
    args = parser.parse_args()

    main(args.READS, args.ANNOTATION, args.SUFFIX)

