import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib_venn import venn3

sequences=list()
csvs=[f for f in os.listdir('./') if f.endswith('final.csv')]
for file in csvs:
    df=pd.read_csv(file,delimiter=' ')
    df=df.drop_duplicates(subset=df.columns[0],keep='first')
    df=df[df.columns[3]].tolist()
    sequences.append(df)
set1=set(sequences[0])
set2=set(sequences[1])
set3=set(sequences[2])
in_three=set1 & set2 & set3

only_in_set1 = set1 - set2 - set3

# Elements only in set2
only_in_set2 = set2 - set1 - set3

# Elements only in set3
only_in_set3 = set3 - set1 - set2

# Elements in set1 and set2 but not in set3
in_set1_and_set2_not_set3 = (set1 & set2) - set3

# Elements in set1 and set3 but not in set2
in_set1_and_set3_not_set2 = (set1 & set3) - set2

# Elements in set2 and set3 but not in set1
in_set2_and_set3_not_set1 = (set2 & set3) - set1
print(len(in_three))
print(len(only_in_set1))
print(len(only_in_set2))
print(len(only_in_set3))
print(len(in_set1_and_set2_not_set3))
print(len(in_set1_and_set3_not_set2))
print(len(in_set2_and_set3_not_set1))
print(csvs)
venn3(subsets=(len(only_in_set1),len(only_in_set2),len(in_set1_and_set2_not_set3),len(only_in_set3),len(in_set1_and_set3_not_set2),len(in_set2_and_set3_not_set1),len(in_three)),set_labels=('KO 24h', 'WT 0h', 'WT 24h'))

plt.savefig('./venn_diagram_total.png',dpi=1000)
