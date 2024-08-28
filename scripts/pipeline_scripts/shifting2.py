
"""
Created on Mon Sep 12 06:37:07 2022
Adaptación de CirSeq a reads de 150 nt
@author: Rocío Marugán (rmaruganpin@alumni.unav.es)
"""

import gzip
import argparse
import os

from utils import rotate, _step_


def main(fastq_file:str, suffix:str=""):

    save_dir, file_root = _step_(fastq_file, suffix)
    
    # Obtención de todas las permutaciones de la unidad repetitiva
    infile = gzip.open(fastq_file, "r")
    outfile = gzip.open(os.path.join(save_dir, file_root + ".fastq.gz"),"wb")

    while True:
        SequenceID = infile.readline().decode().split("_")[1:]
        SequenceID = "_".join(SequenceID)
        if SequenceID == "":
            break
    
        Sequence = infile.readline().decode()
        Sequence = Sequence.rstrip("\n") 
        EmptyLine = infile.readline()
        QualityScores = infile.readline().decode()
        QualityScores = QualityScores.rstrip("\n")
        
        Sequence = rotate(Sequence, SequenceID, EmptyLine, QualityScores, outfile)
        
    infile.close()
    outfile.close()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Perform shifting')
    parser.add_argument('-i', '--INPUT_FASTQ', type=str, required=True, help='Select the path where the input SAM file is located')
    parser.add_argument('-s', '--SUFFIX', type=str, required=False, default="", help='Select the suffix for saving output files with')
    args = parser.parse_args()

    main(args.INPUT_FASTQ, args.SUFFIX)
    
