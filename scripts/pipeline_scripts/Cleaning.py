
"""
Created on April 2024
Grouping of the reads with the same name
@author: Alvaro Perona (aperonalope@alumni.unav.es)
"""


import argparse
import os
import numpy as np

from utils import rotate, _step_


def main(file:str, suffix:str=""):
    # We obtain the full directory and file name without extensions of the input file
    save_dir, file_root = _step_(file, suffix)
    
    # We open the input and output csv
    input_csv = open(file, "r")
    output_csv = open(os.path.join(save_dir, file_root + "_final"+".csv"),"w")

    #We write the header into the output file
    header="Sequence"+" "+"Counts"+" "+"Length"+" "+"Sequence"+" "+"Matches_ID"+" "+"nºMatches"+"\n"
    output_csv.write(header)
    
    #We create an empty diccionary
    diccionary={}
    diccionary_sequencias={}
    #We skip the header of the input file and read a line
    fila=input_csv.readline().rstrip("\n").split(" ")
    fila=input_csv.readline().rstrip("\n").split(" ")
    
    #We create an infinite loop that eands when there are no more lines
    while True:
        if fila[0]=="":
            break
        #The keys of the dictionary are the concatenation of the rotation, read name, count and sequence of a line, separated by a "_"
        strings=fila[0].split("@")[1]+"_"+fila[1]+"_"+fila[2]+"_"+fila[6]
        gene_id=fila[5]
        
        #If the key is already in the dictionary, we append the gene id to that keys values
        
        if strings in diccionary:
            diccionary[strings].append(gene_id)
        #If the read is not a key in the dictionary, we create a new key with the gene id as value
        else:
            diccionary[strings]=[gene_id]
        #We read a new line
        fila=input_csv.readline().rstrip("\n").split(" ")
        

    #We create a list with the keys of the dictionary
    read_names=list(diccionary.keys())
    #We iterate over those keys
    for k in np.arange(0,len(read_names)):
        
        #we separate the fields of the key into strings separated by a space
        #!We can substitute this by simply replacing the "_" by a space
        linea= read_names[k].split("_")[0]+" "+read_names[k].split("_")[1]+" "+read_names[k].split("_")[2]+" "+read_names[k].split("_")[3]+" "
        
        #We extract the values of the key
        
        genes=diccionary[read_names[k]]
        #We initialize the count of genes to 0
        count_genes=0
        linea=linea
        #We iterate over the genes of that key, and we add them to the string separated by a ":" and add 1 to the count of genes
        for gene_id in genes:
            count_genes=count_genes+1
            linea=linea+":"+gene_id
        
        #We also append a space and the count of genes
        linea=linea+" "+str(count_genes)+"\n"
        print(linea)
        #We write the line into the output file
        output_csv.write(linea)
    #We close the input and output files
    print()
    input_csv.close()
    output_csv.close()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Perform shifting')
    parser.add_argument('-i', '--INPUT_FASTQ', type=str, required=True, help='Select the path where the input SAM file is located')
    parser.add_argument('-s', '--SUFFIX', type=str, required=False, default="", help='Select the suffix for saving output files with')
    args = parser.parse_args()

    main(args.INPUT_FASTQ, args.SUFFIX)
    
