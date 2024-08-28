import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

csvs= [f for f in os.listdir('./') if f.endswith('final.csv')]
#print
for files in csvs:
    df=pd.read_csv(files,delimiter=" ")
    df=df.drop_duplicates(subset=df.columns[0],keep='first')
    n_matches=df['Counts'].value_counts().sort_index().to_numpy()
    print(n_matches)
    etiquetas=['1','2','3','4','5','6','>6']
    print(n_matches.tolist())
    print(np.sum(n_matches[6:].tolist()))
    results=n_matches[0:6].tolist()+[np.sum(n_matches[6:]).tolist()]
    print(results)
    total=np.sum(results)
    percentages=[f'{(c/total)*100 :0.1f}%' for c in results]
    colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0','#ffb3e6','#c4e17f']
    labels=[f'{num}({per})' for num,per in zip(etiquetas,percentages)]
    print(percentages)
    plt.figure(figsize=(12,8))
    plt.pie(results, colors=colors,startangle=140)
    plt.axis('equal')
    plt.axis('equal')
    plt.legend(labels, bbox_to_anchor=(0.84, 0.7),title="nº counts", fontsize=18,title_fontsize=24)
    plt.title(files.split('_')[1]+' '+files.split('_')[2],fontsize=28)

    plt.savefig('pie_chart_counts'+'_'+files.split('_')[1]+'_'+files.split('_')[2]+'.png',dpi=1000)
    plt.clf()
