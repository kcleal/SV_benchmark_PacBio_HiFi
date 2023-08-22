import pandas as pd
import json
import matplotlib.pyplot as plt
import seaborn as sns

data = []
for caller in ('sniffles', 'cuteSV', 'delly', 'dysgu', 'NGSEP.vcf_SVsLongReads'):
    for n in range(44, 50):
        with open(f'truvari_{n}_{caller}/summary.json', 'r') as f:
            dta = json.load(f)
        d = {'caller': caller.replace('.vcf_SVsLongReads',''), 'sample': n, 'TP': None}
        d.update({k: v for k, v in dta.items() if k in ('precision', 'recall', 'f1', 'gt_concordance', 'TP-base', 'FP', 'FN')})
        d['TP'] = d['TP-base']
        del d['TP-base']
        data.append(d)
    if not os.path.exists(f'truvari_all_{caller}/summary.json'):
        continue
    with open(f'truvari_all_{caller}/summary.json', 'r') as f:
         dta = json.load(f)
    d = {'caller': caller.replace('.vcf_SVsLongReads',''), 'sample': 'all', 'TP': None}
    d.update({k: v for k, v in dta.items() if k in ('precision', 'recall', 'f1', 'gt_concordance', 'TP-base', 'FP', 'FN')})
    d['TP'] = d['TP-base']
    del d['TP-base']
    data.append(d)

df = pd.DataFrame().from_records(data)
cols = ['caller', 'TP', 'FP', 'FN', 'precision', 'recall', 'f1', 'gt_concordance']
print('Merged runs:')
print(df[df['sample'] == 'all'][cols].round(4).to_markdown(index=False))
print('Single run average:')
print(df[df['sample'] != 'all'][cols].groupby('caller', sort=False).mean().round(3).to_markdown())
print('Single runs:')
print(df[df['sample'] != 'all'][cols].round(4).to_markdown(index=False))

sns.set_style("whitegrid")
plt.figure(figsize=(7, 4))
plt.subplots_adjust(right=0.8)
ax = sns.scatterplot(df, x="recall", y="precision", hue="caller", style="sample", s=100)
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.xlim(0.75, 1); plt.ylim(0.92, 0.98)
plt.savefig("benchmark_result.png")
plt.show()