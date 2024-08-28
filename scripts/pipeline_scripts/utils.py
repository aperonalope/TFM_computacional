
import os


#Useful functions

def rotate(Sequence, ID, EmptyLine, QualityScores, outfile):
    for i in range(len(Sequence)):
        SequenceID_2 = "@" + "rotation" + str(i)+"_"+ ID
        outfile.write(SequenceID_2.encode())

        Sequence_2 = Sequence[-i:] + Sequence[:-i] + "\n"
        outfile.write(Sequence_2.encode())

        outfile.write(EmptyLine)
        
        QualityScores_2 = QualityScores + "\n"
        outfile.write(QualityScores_2.encode())


# Call functions at every step of the pipeline to manage file names and saving location

def get_dir(file_path):
    return os.path.dirname(os.path.abspath(file_path))

def get_file_root(file_path, suffix, sep="_"):
    if suffix != "":
        return "".join([os.path.basename(file_path).split(".")[0], sep, suffix])
    else:
        return os.path.basename(file_path).split(".")[0]

def _step_(file_path, suffix):
    save_dir = get_dir(file_path)
    file_root = get_file_root(file_path, suffix)

    return save_dir, file_root

    
