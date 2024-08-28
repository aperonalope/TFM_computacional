import gzip
import argparse
import os

from utils import _step_

def main(sam_file:str, fastq_file:str, suffix:str=""):

    save_dir, file_root = _step_(sam_file, suffix)

    # Abre el archivo SAM y almacena los nombres de secuencias que deseas eliminar
    sequences_to_remove = set()
    with open(sam_file, "r") as sam:
        for line in sam:
            if not line.startswith("@"):
                fields = line.split("\t")
                qname = fields[0]  # Obtenemos el campo QNAME
                sequences_to_remove.add(qname)

    # Abre el archivo FASTQ comprimido en gzip y un nuevo archivo FASTQ comprimido en gzip para escribir las secuencias filtradas
    with gzip.open(fastq_file, "rt") as fastq, gzip.open(os.path.join(save_dir, file_root + "_f.fastq.gz"), "wt") as filtered_fastq:
        buffer = []  # Almacena temporalmente las 4 líneas de cada secuencia
        num = 0

        for line in fastq:
            buffer.append(line)
            
            # Cuando hemos almacenado las 4 líneas de una secuencia
            if len(buffer) == 4:
                # Modificamos el encabezado para eliminar la parte final
                fastq_header = buffer[0].strip().rsplit(' ', 1)[0][1:]
                
                # Comprobamos si el encabezado coincide con los nombres de secuencias a eliminar
                if fastq_header in sequences_to_remove:
                    num += 1
                # Si no coincide, escribimos las líneas en el archivo filtrado
                else:
                    filtered_fastq.writelines(buffer)
                
                buffer = []  # Reiniciamos el buffer para la siguiente secuencia
                
    print(f"{num} sequences removed")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Remove duplicates')
    parser.add_argument('-i', '--INPUT_SAM', type=str, required=True, help='Select the path where the input SAM file is located')
    parser.add_argument('-f', '--INPUT_FASTQ', type=str, required=True, help='Select the path where the input FASTQ file is located')
    parser.add_argument('-s', '--SUFFIX', type=str, required=False, default="", help='Select the suffix for saving output files with')
    args = parser.parse_args()

    main(args.INPUT_SAM, args.INPUT_FASTQ, args.SUFFIX)
