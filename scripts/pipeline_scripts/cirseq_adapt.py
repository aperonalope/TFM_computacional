"""
Created on July of 2024
@author: Alvaro Perona (aperonalope@hpclogin.uanv.es)
"""
import argparse
import numpy as np
from scipy.stats import mode
import gzip
import os
from pathlib import Path

from utils import rotate, _step_


# Obtención de todas las permutaciones de la unidad repetitiva

#We specify the read length
SEQUENCE_LENGTH = 150
#We specify the number of nucleotides of the splitter by which we will split the read to idenity repetitive units
splitter=15

def main(fastq_file:str, suffix:str):
    
    #We open the input and output files
    save_dir, file_root = _step_(fastq_file, suffix)
    infile = gzip.open(fastq_file, "rb")

    outfile1 = gzip.open(os.path.join(save_dir, file_root + ".fastq.gz"),"wb")
    outfile3 = gzip.open(os.path.join(save_dir, file_root + "_rotated.fastq.gz"),"wb")
    while True:
        
        
        #We break the infinite loop when the line is empty
        if SequenceID == "":
            break
        #We store each line of a FASTQ read into a separate variable
        SequenceID = infile.readline().decode().rstrip("\n")
        Sequence = infile.readline().decode().rstrip("\n") 
        EmptyLine = infile.readline()
        QualityScores = infile.readline().decode().rstrip("\n")
        
        
        
        counter=0
        pos=0
        colgao="nada"
        #Identificación de todos los posibles sets de repeticiones y guardado de su longitud
        for i in np.arange(0, splitter):
            
            Substrings = Sequence.split(Sequence[i:i+splitter])
            
            #If the splitting found the splitter at two points  it will have at least 3 elements in the list of substrings
            #? Maybe len substrings mayor que uno will suffice
            if counter==0 and (len(Substrings)>2):
                
                #We change the counter so that once a valid splitting is found, otherones are not tried by the for loop
                counter=counter+1
                #This variable is the string of nucleotides of the splitter
                sp=Sequence[i:i+splitter]
                
                #We find the consensus sequence as the sum of the splitter plus the remaining nucleotides in the second element of the list
                ConsensusSequence=Sequence[i:splitter+len(Substrings[1])]
                #The initial position where the splitter was found (for clarity)
                pos=i
                #The nucleotides at the last element of the substring list
                colgao=Substrings[-1]
                
                #Number of substrings
                numberiter=len(Substrings)
                #We initialize the hypothethical read
                suppopsedly=""
                #We initialize a counter and a for loop that iterates over the substrings
                conta=1
                for string in Substrings:
                    #This are the posible actions when initializing the hypothetical read
                    if conta==1:
                        #If the string is empy it means that the splitter was found at the initial nucleotide, so we dont do nothing
                        if string=="":
                            suppopsedly=suppopsedly
                        #If it is not empty we sum take as many nucleoitdes as letters in this first string from the last nucleotides of the consensus sequence
                        else:
                            
                            suppopsedly=suppopsedly+ConsensusSequence[-len(string):]
                    #This are the posible actions when handling the last part of the hypothethical read. There are a number of different scenarios and we have conditionals to handle all of them
                    elif conta==numberiter:
                        if string=="":
                            
                            suppopsedly=suppopsedly+sp
                        elif len(string)==splitter:
                            if len(ConsensusSequence)<splitter*2:
                                
                                suppopsedly=suppopsedly+ConsensusSequence+ConsensusSequence[:splitter*2-len(ConsensusSequence)]
                            else:
                                suppopsedly=suppopsedly+ConsensusSequence[:splitter*2]
                        elif (len(string)+splitter)>len(ConsensusSequence):
                            suppopsedly=suppopsedly+ConsensusSequence+ConsensusSequence[:(len(string)+splitter)-len(ConsensusSequence)]
                        else:
                            
                            suppopsedly=suppopsedly+ConsensusSequence[:len(string)+splitter]
                    #When constucting the middle section we just append the consensus sequence
                    else:
                        
                        suppopsedly=suppopsedly+ConsensusSequence
                    #We increase the counter so that we know at which point of the hypothetical read construction we are in
                    conta=conta+1
                    
                   
        if colgao==Sequence[pos+splitter:pos+splitter+len(colgao)] and ConsensusSequence.count(ConsensusSequence[0])!=len(ConsensusSequence):
            
            
            
            
            #We compare the hypothetical read and the actual read
            if suppopsedly==Sequence:
                
                #Puede haber unidades repetitivas menores de 15 nucleotidos. No interesan de momento. Comporbamos si es e caso y si lo es las guardamos en un output aparte para ser estudiadas mas adleante
                comprobacion=ConsensusSequence[0:round(len(ConsensusSequence)/3)-1]
                if ConsensusSequence.split(comprobacion)!=["",ConsensusSequence[round(len(ConsensusSequence)/3)-1:]]:
                    g=0
                    
                    if all(element=="" for element in ConsensusSequence.split(comprobacion+ConsensusSequence.split(comprobacion)[1]))==True:
                        g=1
                        RepeatLength=len(ConsensusSequence)
                        SequenceID = SequenceID + ":" + str(RepeatLength) + "\n"
                        ConsensusSequence=ConsensusSequence + "\n"
                        QualityScores = QualityScores[:RepeatLength] + "\n"
                        outfile3.write(SequenceID.encode()) 
                        outfile3.write(ConsensusSequence.encode()) 
                        outfile3.write(EmptyLine) 
                        outfile3.write(QualityScores.encode())
                    else:
                        
                        RepeatLength=len(ConsensusSequence)
                        SequenceID=SequenceID+":"+str(RepeatLength)+"\n"
                        ConsensusSequence=ConsensusSequence + "\n"
                        QualityScores = QualityScores[:RepeatLength] + "\n"
                        
                        outfile1.write(SequenceID.encode())
                        outfile1.write(ConsensusSequence.encode())
                        outfile1.write(EmptyLine)
                        outfile1.write(QualityScores.encode())
                #Si la unidad repetitiva, not tiene  subunidades repetitivas dentro de ella, las guardamos.
                else:
                    
                    RepeatLength=len(ConsensusSequence)
                    SequenceID = SequenceID + ":" + str(RepeatLength) + "\n"
                    ConsensusSequence=ConsensusSequence + "\n"
                    QualityScores = QualityScores[:RepeatLength] + "\n"
                    
                    outfile1.write(SequenceID.encode())
                    outfile1.write(ConsensusSequence.encode())
                    outfile1.write(EmptyLine)
                    outfile1.write(QualityScores.encode())






    infile.close()
    outfile1.close()
    
    
    
    
    outfile3.close()






if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Algorithm to detect repetitive units in reads')
    parser.add_argument('-i', '--INPUT_FASTQ', type=str, required=True, help='Select the path where the input FASTQ file is located')
    parser.add_argument('-s', '--SUFFIX', type=str, required=True, help='Select the suffix for saving output files with')
    args = parser.parse_args()

    main(args.INPUT_FASTQ, args.SUFFIX)
