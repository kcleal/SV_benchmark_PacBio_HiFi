import pandas as pd
import json
import matplotlib.pyplot as plt
import seaborn as sns

data = []
for caller in ('sniffles', 'cuteSV', 'delly', 'dysgu', 'NGSEP.vcf_SVsLongReads', 'svim'):
    for n in range(44, 50): 
        with open(f'truvari_{n}_{caller}/summary.json', 'r') as f:
            dta = json.load(f)
        d = {'caller': caller.replace('.vcf_SVsLongReads',''), 'sample': n, 'TP': None}
        d.update({k: v for k, v in dta.items() if k in ('precision', 'recall', 'f1', 'gt_concordance', 'TP-base', 'FP', 'FN')})
        d['TP'] = d['TP-base']
        del d['TP-base']
        d['depth'] = pd.read_csv(f'SRR103822{n}.cov.tsv', sep='\t')['meandepth'].iloc[0]
        data.append(d)
    with open(f'truvari_all_{caller}/summary.json', 'r') as f:
         dta = json.load(f)
    d = {'caller': caller.replace('.vcf_SVsLongReads',''), 'sample': 'merged', 'TP': None}
    d.update({k: v for k, v in dta.items() if k in ('precision', 'recall', 'f1', 'gt_concordance', 'TP-base', 'FP', 'FN')})
    d['TP'] = d['TP-base']
    del d['TP-base']
    d['depth'] = pd.read_csv(f'all.cov.tsv', sep='\t')['meandepth'].iloc[0]
    data.append(d)

df = pd.DataFrame().from_records(data)
cols = ['caller', 'TP', 'FP', 'FN', 'precision', 'recall', 'f1', 'gt_concordance']
print('Merged runs:')
print(df[df['sample'] == 'merged'][cols].round(4).to_markdown(index=False))
print('Single run average:')
print(df[df['sample'] != 'merged'][cols].groupby('caller', sort=False).mean().round(3).to_markdown())
print('Single runs:')
cols.append('depth')
print(df[df['sample'] != 'merged'][cols].round(4).to_markdown(index=False))

sns.set_style("whitegrid")
plt.figure(figsize=(8, 4))
plt.subplots_adjust(right=0.75)
ax = sns.scatterplot(df, x="recall", y="precision", hue="caller", style="sample", s=100)
ax.legend(bbox_to_anchor=(1.03, 0.5), loc='center left')
plt.savefig("benchmark_result.png")


plt.figure(figsize=(8, 4))
plt.subplots_adjust(right=0.75)
ax = sns.scatterplot(df[df['sample'] != 'merged'], x="depth", y="f1", hue="caller", style="sample", s=100)
ax.legend(bbox_to_anchor=(1.03, 0.5), loc='center left')
plt.ylim(0.86, 0.96)
plt.savefig("benchmark_result_f1.png")
plt.show()

