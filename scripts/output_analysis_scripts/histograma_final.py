import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

csvsss = [f for f in os.listdir('./') if f.endswith('.csv')]
colors = ['b', 'g', 'r']
counter = 0

for file in csvsss:
    csv = pd.read_csv(file, sep=' ')
    print(file)
    csv_unique = csv.drop_duplicates(subset=csv.columns[0], keep='first')
    
    read_name = csv_unique.iloc[:,0].to_numpy()
    counts = csv_unique.iloc[:,1].to_numpy()
    lengths = csv_unique.iloc[:,2].to_numpy()
    UNIQ_LENGTHS = np.unique(lengths)

    summed_coun = np.zeros(len(UNIQ_LENGTHS))

    for i, cou in enumerate(UNIQ_LENGTHS):
        summed_coun[i] = np.sum(counts[lengths == cou])

    print(summed_coun)
    plt.figure(figsize=(12,8))
    plt.bar(UNIQ_LENGTHS, summed_coun, alpha=0.5, color=colors[counter], label=file.split("_")[1] + " " + file.split("_")[2])
    
    # Increased font sizes by 100%
    plt.xlabel('Lengths', fontsize=24)
    plt.ylabel('Summed Counts', fontsize=24)
    plt.ylim(0, 400)
    plt.title('Summed Counts by Length', fontsize=28)
    plt.legend(loc='upper left', fontsize=18)
    
    # Increase tick size
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    

    plt.savefig('barplot_countss_by_length' + file.split("_")[1] + file.split("_")[2] + '.png',dpi=600)
    plt.clf()
    counter = counter + 1



