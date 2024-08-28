import re
import argparse
import os

from utils import _step_

def main(sam_file:str, suffix:str=""):

    save_dir, file_root = _step_(sam_file, suffix)

    # Abre el archivo SAM y almacena las secuencias en las tres categorías
    f = open(sam_file, "r")

    outfile_no_softclip = open(os.path.join(save_dir, file_root + "_no_softclip.sam"), "w")
    outfile_one_softclip = open(os.path.join(save_dir, file_root + "_one_softclip.sam"), "w")
    outfile_both_softclip = open(os.path.join(save_dir, file_root + "_both_softclip.sam"), "w")

    totalcount = 0
    no_softclip_count = 0
    one_softclip_count = 0
    both_softclip_count = 0
    softclipmin = 2  # Mínima longitud de softclip requerida

    while True:
        line = f.readline()
        record = line.strip('\n').split()
        if not record:
            break
        if record[0][0] != "@":
            totalcount += 1
            cigar = record[5]
            beginning_eval = "S" in cigar[1:4]
            # Verificamos si "S" está presente tanto al principio como al final de CIGAR
            if beginning_eval and cigar[-1] == "S":
                both_softclip_count += 1
                outfile_both_softclip.write(line)
            elif beginning_eval or cigar[-1] == "S":
                one_softclip_count += 1
                outfile_one_softclip.write(line)
            else:
                no_softclip_count += 1
                outfile_no_softclip.write(line)

    f.close()
    outfile_no_softclip.close()
    outfile_one_softclip.close()
    outfile_both_softclip.close()

    print("Total de secuencias procesadas:", totalcount)
    print("Secuencias sin softclip:", no_softclip_count)
    print("Secuencias con softclip en un lado:", one_softclip_count)
    print("Secuencias con softclip en ambos lados:", both_softclip_count)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='CirSeq adaptation for 150-nt-long reads')
    parser.add_argument('-i', '--INPUT_SAM', type=str, required=True, help='Select the path where the input SAM file is located')
    parser.add_argument('-s', '--SUFFIX', type=str, required=False, default="", help='Select the suffix for saving output files with')
    args = parser.parse_args()

    main(args.INPUT_SAM, args.SUFFIX)

