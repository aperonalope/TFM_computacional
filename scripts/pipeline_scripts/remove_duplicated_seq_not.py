import gzip
import argparse
import os

from utils import _step_

def main(sam_file:str, fastq_file:str, suffix:str=""):
    save_dir, file_root = _step_(sam_file, suffix)

    # Abre el archivo SAM y almacena los nombres de secuencias que deseas eliminar
    listaID=[]
    
    with open(sam_file, "r") as sam:

        for line in sam:
            
            
            if line.startswith("rotation"):
                f=line.split("_")[1]
                listaID.append(f)
    listaID=list(set(listaID))
    # Abre el archivo FASTQ comprimido en gzip y un nuevo archivo FASTQ comprimido en gzip para escribir las secuencias filtradas
    anterior=["vivae_lvino"]
    with gzip.open(fastq_file, "rt") as fastq, gzip.open(os.path.join(save_dir, file_root + "_f.fastq.gz"), "wt") as filtered_fastq:
        buffer = []  # Almacena temporalmente las 4 líneas de cada secuencia
        for line in fastq:
            buffer.append(line)
            if len(buffer)==4:
                header=buffer[0].split("_")[1]
                if header!= anterior[0].split("_")[1]:
                    anterior=buffer
                    if header in listaID:
                        filtered_fastq.writelines(buffer)
                buffer=[]
                



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Remove duplicates')
    parser.add_argument('-i', '--INPUT_SAM', type=str, required=True, help='Select the path where the input SAM file is located')
    parser.add_argument('-f', '--INPUT_FASTQ', type=str, required=True, help='Select the path where the input FASTQ file is located')
    parser.add_argument('-s', '--SUFFIX', type=str, required=False, default="", help='Select the suffix for saving output files with')
    args = parser.parse_args()

    main(args.INPUT_SAM, args.INPUT_FASTQ, args.SUFFIX)
