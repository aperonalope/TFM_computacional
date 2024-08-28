#SCRIPT PARA ELIMINAR CONTAR Y ELIMINAR SECUENCIAS IGUALES

from collections import Counter
from scipy.stats import mode
import gzip
import argparse
import os

import pandas as pd

from Bio.Seq import Seq

from utils import _step_


#Función para realizar la Barrow Wheeler Transform, que genera un tag que nos permitirá buscar secuencias repetidas más fácilmente
def bwt(a):
    words = list(a)
    bwt = []

    for i in range(len(words)):
        word = a[-1] + a[:-1]
        new = ''.join(word)
        a = new
        bwt.append(new)

    sort = sorted(bwt)

    output = []
    for el in sort:
        last = el[-1]
        output.append(last)

    return "".join(output)


def main(fastq_file:str, suffix:str=""):

    save_dir, file_root = _step_(fastq_file, suffix)

    infile = gzip.open(fastq_file, "r")
 
    sequence_list = []
    sequence_dict = {}
    id_list = [] 
    qscores_list = []

    #Vamos leyendo el fichero de las secuencias detectadas como repetitivas y creamos una lista con los tags generados por la BWT
    #También generamos un diccionario con cada secuencia y su correspondiente tag de BWT para luego poder realizar el paso inverso

    while True:
        SequenceID = infile.readline().decode().rstrip("\n").split(" ")[0]

        if SequenceID == "":
            break
            
        #leemos cada linea y guardamos las variables
        Sequence = infile.readline().decode().rstrip("\n") 
        EmptyLine = infile.readline()
        QualityScores = infile.readline().decode().rstrip("\n")
        
        #Guardamos la BWT de la secuencia y añadimos un registro más al diccionario con secuencia y tag
        result = bwt(Sequence)
        sequence_list.append(result)
        sequence_dict[result] = Sequence
        
        #Guardamos los IDs y los QualityScores de cada secuencia de manera ordenada
        id_list.append(SequenceID)
        qscores_list.append(QualityScores)

    infile.close()


    #Dado que las secuencias solo tienen el mismo tag si son iguales o son iguales pero con un origen diferente, para eliminar  
    #duplicados y contarlos bastará con usar la función Counter sobre la lista de tags
    sequence_counts = Counter(sequence_list)

    #Generar fastq

    result_list = []
    output_file = gzip.open(os.path.join(save_dir, file_root + "_filtered.fastq.gz"), "wt")

    for sequence in sequence_counts:
        #Seleccionamos el SequenceID correspondiente a la primera posición que se halla al buscar cada secuencia en la lista del 
        #contador
        pos = sequence_list.index(sequence)
        id_seq = id_list[pos]
        qscore = qscores_list[pos]
        
        #Obtenemos la cuenta asociada al tag
        count = sequence_counts[sequence]
        
        #Obtenemos la secuencia repetitiva a la que se asocia el tag
        seq = sequence_dict[sequence]
        
        record = "{}_{}\n{}\n+\n{}\n".format(id_seq,count,seq, qscore)    
        output_file.write(record)
        result_list.append((count, seq, id_seq, qscore ))
    output_file.close()


    # Crear un DataFrame a partir de result_list
    df = pd.DataFrame(result_list, columns=['Count', 'Sequence', 'SequenceID', 'QScore'])

    # Exportar el DataFrame a un archivo Excel
    df.to_excel(os.path.join(save_dir, file_root + '_filtered.xlsx'), index=False)

    n_most_common = 10
    #print(f"Most {n_most_common} common sequences: {[s[0] for s in sequence_counts.most_common(n_most_common)]}")
    #print(f"Number of sequences: {len(sequence_counts)}")
    #print(f"Number of reads: {sum(sequence_counts.values())}")
    #print(f"Average reads per seqeunce: {sum(sequence_counts.values())/len(sequence_counts)}")


    #Buscar coincidencias con el complemento invertido y eliminarlas
    #Hacemos listas auxiliares para luego poder eliminar los duplicados y mantener el orden

    secuencias = [elemento[1] for elemento in result_list]
    contadores = [elemento[0] for elemento in result_list]
    ids = [elemento[2] for elemento in result_list]
    qscores = [elemento[3] for elemento in result_list]
    secuencias_bwt = list(sequence_counts.keys())

    for i, secuencia in enumerate(secuencias):

        secuencia_obj = Seq(secuencia)
        #Hallamos el complemento invertido
        complemento_invertido = secuencia_obj.reverse_complement()
        
        #Realizamos la BWT
        complemento_invertido = bwt(complemento_invertido)
        
        if complemento_invertido in secuencias_bwt:
            coincidencias = [idx for idx, x in enumerate(secuencias_bwt) if x == complemento_invertido] 
            
            # Eliminar las coincidencias de la lista en orden inverso para no alterar las posiciones de la lista
            for index in reversed(coincidencias):
                # Sumamos los duplicados que vamos a borrar al contador de la secuencia en cuestión
                contadores[i] = contadores[i] + contadores[index]
                
                # Eliminamos los duplicados de las listas
                secuencias_bwt.pop(index)            
                secuencias.pop(index)
                contadores.pop(index)
                qscores.pop(index)
                ids.pop(index)
                    

    data = {
        'secuencias': secuencias,
        'contadores': contadores,
        'ids': ids,
        'qscores': qscores
    }

    # Crear el DataFrame a partir del diccionario
    df2 = pd.DataFrame(data)

    # Guardar el DataFrame en un archivo FASTQ (Con etiqueta "final")
    with gzip.open(os.path.join(save_dir, file_root + '_final.fastq.gz'), 'wt') as f:
        for _, row in df2.iterrows():
            f.write(f'{row["ids"]}_{row["contadores"]}\n')
            f.write(f'{row["secuencias"]}\n')
            f.write('+\n')
            f.write(f'{row["qscores"]}\n')

    # Exportar el DataFrame a un archivo Excel
    df2.to_excel(os.path.join(save_dir, file_root + '_final.xlsx'), index=False)



if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Filter duplicated reads')
    parser.add_argument('-i', '--INPUT_FASTQ', type=str, required=True, help='Select the path where the input SAM file is located')
    parser.add_argument('-s', '--SUFFIX', type=str, required=False, default="", help='Select the suffix for saving output files with')
    args = parser.parse_args()

    main(args.INPUT_FASTQ, args.SUFFIX)

