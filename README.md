:bar_chart: SV Benchmark using PacBio HiFi data
===============================================

This is a reproducible benchmark of structural variant (SV) callers using PacBio HiFi reads. The truth set is described in detail here:

[A robust benchmark for detection of germline large deletions and insertions. Zook et al., 2020. Nature Biotechnology](https://www.nature.com/articles/s41587-020-0538-8)


![plot](./benchmark_result.png)
![plot](./benchmark_result_f1.png)

Single run average
------------------
  
| caller   |      TP |      FP |       FN |   precision |   recall |    f1 |   gt_concordance |
|:---------|--------:|--------:|---------:|------------:|---------:|------:|-----------------:|
| sniffles | 8659.83 | 564     |  981.167 |       0.939 |    0.898 | 0.918 |            0.936 |
| cuteSV   | 8490.17 | 571.833 | 1150.83  |       0.937 |    0.881 | 0.908 |            0.92  |
| delly    | 7864.5  | 297     | 1776.5   |       0.964 |    0.816 | 0.883 |            0.911 |
| dysgu    | 8955.83 | 481.5   |  685.167 |       0.949 |    0.929 | 0.939 |            0.929 |
| NGSEP    | 8742.5  | 575.833 |  898.5   |       0.938 |    0.907 | 0.922 |            0.926 |
| svim     | 8838.33 | 726.167 |  802.667 |       0.924 |    0.917 | 0.92  |            0.841 |

Merged runs
-----------
| caller   |   TP |   FP |   FN |   precision |   recall |     f1 |   gt_concordance |
|:---------|-----:|-----:|-----:|------------:|---------:|-------:|-----------------:|
| sniffles | 9393 |  589 |  248 |      0.941  |   0.9743 | 0.9573 |           0.9772 |
| cuteSV   | 9367 |  548 |  274 |      0.9447 |   0.9716 | 0.958  |           0.9863 |
| delly    | 9320 |  625 |  321 |      0.9372 |   0.9667 | 0.9517 |           0.9827 |
| dysgu    | 9421 |  538 |  220 |      0.946  |   0.9772 | 0.9613 |           0.9728 |
| NGSEP    | 9300 |  624 |  341 |      0.9371 |   0.9646 | 0.9507 |           0.9682 |
| svim     | 9346 |  699 |  295 |      0.9304 |   0.9694 | 0.9495 |           0.9766 |


Reads were from PacBio Sequel II HiFi. Six runs were tested individually (~8X coverage), or after merging together (~51X coverage). SV callers tested were as follows:

- [sniffles v2.2.0](https://github.com/fritzsedlazeck/Sniffles)

- [cuteSV v2.0.3](https://github.com/tjiangHIT/cuteSV)

- [delly v1.1.6](https://github.com/dellytools/delly)

- [dysgu v1.6.0](https://github.com/kcleal/dysgu)

- [NGSEP v4.3.2](https://github.com/NGSEP/NGSEPcore)

- [svim v2.0.0](https://github.com/eldariont/svim)


For benchmarking [truvari v4.0.0](https://github.com/ACEnglish/truvari) was used with parameters `-r 1000 --passonly -p 0 --dup-to-ins`

Run the ``run_single_sample.sh`` script to map and call SVs for each sequencing run. Wait for that to complete then run the ``run_merged.sh`` script. 

### Requirements:

- ~ 350 Gb space, 64 gb Ram, 8 cores
- Docker / Singularity (Required for Mac, optional for Linux)
- Job scheduler of some kind (slurm used here)


## Setup environment

```
mkdir benchmark && cd benchmark
docker run -it --memory="64g" --mount src="${PWD}",target=/results,type=bind condaforge/mambaforge
mamba update conda -y && cd results
```
Note, you may need to set the memory and swap space manually using Docker Desktop on Mac.

Install tools:

```
mamba create -c bioconda -c conda-forge -n bench python=3.9 awscli samtools=1.17 bcftools minimap2=2.26 sniffles=2.2.0 cuteSV=2.0.3 truvari=4.0.0 delly=1.1.6 -y
conda activate bench
pip install dysgu==1.6.0 svim==2.0.0
wget https://github.com/NGSEP/NGSEPcore/releases/download/v4.3.2/NGSEPcore_4.3.2.jar
```

## Grab datasets

Reference genome:
```
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh37/hs37d5.fa.gz
gunzip hs37d5.fa.gz && samtools faidx hs37d5.fa
```

PacBio reads found online at SRA or EBI under accession [PRJNA586863](https://www.ebi.ac.uk/ena/browser/view/PRJNA586863):
```
for i in {44..49}; do wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/0${i}/SRR103822${i}/SRR103822${i}_subreads.fastq.gz; done
```

SV truth set:
```
ftp=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6
wget ${ftp}/HG002_SVs_Tier1_v0.6.bed
wget ${ftp}/HG002_SVs_Tier1_v0.6.vcf.gz
wget ${ftp}/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
```

## Run SV callers

Map data and call SVs for each run:

```
for i in {44..49}; do \

$(cat > job${i}.sh << EOF
#!/bin/bash
#SBATCH --job-name=joball
#SBATCH --partition=compute
#SBATCH --output=out.${i}.%J
#SBATCH --error=err.${i}.%J
#SBATCH --time=0-24:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

minimap2 -t8 -a -x map-hifi hs37d5.fa SRR103822${i}_subreads.fastq.gz | \
  samtools view -bh - | \
  samtools sort -o SRR103822${i}.mm2.bam -
samtools index -@8 SRR103822${i}.mm2.bam

samtools coverage -r 1 SRR103822${i}.mm2.bam > SRR103822${i}.cov.tsv

sniffles --threads 8 --input SRR103822${i}.mm2.bam --vcf HG002_${i}.pacbio.sniffles.vcf

dysgu call --mode pacbio --procs 8 -x --clean hs37d5.fa wd_${i} SRR103822${i}.mm2.bam -o HG002_${i}.pacbio.dysgu.vcf

mkdir wd_cuteSV_${i}
cuteSV -t 8 -s 2 --genotype SRR103822${i}.mm2.bam hs37d5.fa HG002_${i}.pacbio.cuteSV.vcf wd_cuteSV_${i}

delly lr -g hs37d5.fa SRR103822${i}.mm2.bam > HG002_${i}.pacbio.delly.vcf

java -jar NGSEPcore_4.3.2.jar SingleSampleVariantsDetector -runOnlySVs -runLongReadSVs -i SRR103822${i}.mm2.bam -r hs37d5.fa -o HG002_${i}.pacbio.NGSEP.vcf

svim alignment --sample SRR103822${i} --tandem_duplications_as_insertions --interspersed_duplications_as_insertions wd_svim_SRR103822${i} SRR103822${i}.mm2.bam hs37d5.fa
bcftools view -i 'QUAL >= 2' wd_svim_SRR103822${i}/variants.vcf |  bcftools sort -Ov -o HG002_${i}.pacbio.svim.vcf -

EOF
)
  # change this to 'bash job${i}' if running without slurm
  sbatch job${i}.sh
done
```

Merge runs and call SVs:

```
$(cat > joball.sh << EOF
#!/bin/bash
#SBATCH --job-name=joball
#SBATCH --partition=compute
#SBATCH --output=out.all.%J
#SBATCH --error=err.all.%J
#SBATCH --time=0-24:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

samtools merge -@8 -o all.bam *mm2.bam; samtools index -@8 all.bam

samtools coverage -r 1 all.bam > all.cov.tsv

dysgu call --mode pacbio --procs 8 -x hs37d5.fa wd_all all.bam -o HG002_all.pacbio.dysgu2.vcf

sniffles --threads 8 --input all.bam --vcf HG002_all.pacbio.sniffles.vcf

mkdir wd_cuteSV_all
cuteSV -t 8 --genotype all.bam hs37d5.fa HG002_all.pacbio.cuteSV.vcf wd_cuteSV_all

delly lr -g hs37d5.fa all.bam > HG002_all.pacbio.delly.vcf

java -jar -Xmx64G NGSEPcore_4.3.2.jar SingleSampleVariantsDetector -runOnlySVs -runLongReadSVs -i all.bam -r hs37d5.fa -o HG002_all.pacbio.NGSEP.vcf

svim alignment --sample all --tandem_duplications_as_insertions --interspersed_duplications_as_insertions wd_svim_all all.bam hs37d5.fa
bcftools view -i 'QUAL >= 10' wd_svim_all/variants.vcf |  bcftools sort -Ov -o HG002_all.pacbio.svim.vcf -

EOF
)
# run 'bash joball.sh' if running without slurm
sbatch joball.sh
```


## Benchmark

Run truvari:
```
callers=( "sniffles" "cuteSV" "dysgu" "delly" "NGSEP.vcf_SVsLongReads" "svim" )

for i in {44..49}
do
  for name in "${callers[@]}"
  do
    bgzip -f HG002_${i}.pacbio.${name}.vcf
    tabix -f HG002_${i}.pacbio.${name}.vcf.gz
    truvari bench -f hs37d5.fa \
                  -b HG002_SVs_Tier1_v0.6.vcf.gz \
                  --includebed HG002_SVs_Tier1_v0.6.bed \
                  -c HG002_${i}.pacbio.${name}.vcf.gz \
                  --passonly -r 1000 -p 0 --dup-to-ins \
                  -o truvari_${i}_${name}
  done
done


for name in "${callers[@]}"
do
  bgzip -f HG002_all.pacbio.${name}.vcf
  tabix -f HG002_all.pacbio.${name}.vcf.gz
  truvari bench -f hs37d5.fa \
                -b HG002_SVs_Tier1_v0.6.vcf.gz \
                --includebed HG002_SVs_Tier1_v0.6.bed \
                -c HG002_all.pacbio.${name}.vcf.gz \
                --passonly -r 1000 -p 0 --dup-to-ins \
                -o truvari_all_${name}
done
```

Run the plotting script:

```python3 plot_benchmark.py```

Individual samples
------------------
| caller   |   TP |   FP |   FN |   precision |   recall |     f1 |   gt_concordance |   depth |
|:---------|-----:|-----:|-----:|------------:|---------:|-------:|-----------------:|--------:|
| sniffles | 8569 |  556 | 1072 |      0.9391 |   0.8888 | 0.9132 |           0.9309 |  7.6135 |
| sniffles | 8608 |  550 | 1033 |      0.9399 |   0.8929 | 0.9158 |           0.931  |  8.0543 |
| sniffles | 8607 |  588 | 1034 |      0.9361 |   0.8927 | 0.9139 |           0.9367 |  8.4205 |
| sniffles | 8723 |  571 |  918 |      0.9386 |   0.9048 | 0.9214 |           0.9413 |  8.8387 |
| sniffles | 8747 |  559 |  894 |      0.9399 |   0.9073 | 0.9233 |           0.9375 |  9.2231 |
| sniffles | 8705 |  560 |  936 |      0.9396 |   0.9029 | 0.9209 |           0.9393 |  9.0572 |
| cuteSV   | 8396 |  539 | 1245 |      0.9397 |   0.8709 | 0.904  |           0.9136 |  7.6135 |
| cuteSV   | 8424 |  550 | 1217 |      0.9387 |   0.8738 | 0.9051 |           0.9099 |  8.0543 |
| cuteSV   | 8435 |  587 | 1206 |      0.9349 |   0.8749 | 0.9039 |           0.9195 |  8.4205 |
| cuteSV   | 8575 |  588 | 1066 |      0.9358 |   0.8894 | 0.912  |           0.9244 |  8.8387 |
| cuteSV   | 8571 |  583 | 1070 |      0.9363 |   0.889  | 0.9121 |           0.9261 |  9.2231 |
| cuteSV   | 8540 |  584 | 1101 |      0.936  |   0.8858 | 0.9102 |           0.9234 |  9.0572 |
| delly    | 7600 |  289 | 2041 |      0.9634 |   0.7883 | 0.8671 |           0.9314 |  7.6135 |
| delly    | 7693 |  272 | 1948 |      0.9659 |   0.7979 | 0.8739 |           0.9271 |  8.0543 |
| delly    | 7781 |  287 | 1860 |      0.9644 |   0.8071 | 0.8788 |           0.8905 |  8.4205 |
| delly    | 8002 |  311 | 1639 |      0.9626 |   0.83   | 0.8914 |           0.8869 |  8.8387 |
| delly    | 8098 |  310 | 1543 |      0.9631 |   0.84   | 0.8973 |           0.9143 |  9.2231 |
| delly    | 8013 |  313 | 1628 |      0.9624 |   0.8311 | 0.892  |           0.9141 |  9.0572 |
| dysgu    | 8848 |  433 |  793 |      0.9533 |   0.9177 | 0.9352 |           0.925  |  7.6135 |
| dysgu    | 8861 |  462 |  780 |      0.9504 |   0.9191 | 0.9345 |           0.9216 |  8.0543 |
| dysgu    | 8939 |  488 |  702 |      0.9482 |   0.9272 | 0.9376 |           0.9264 |  8.4205 |
| dysgu    | 9015 |  493 |  626 |      0.9481 |   0.9351 | 0.9416 |           0.9322 |  8.8387 |
| dysgu    | 9062 |  505 |  579 |      0.9472 |   0.9399 | 0.9436 |           0.933  |  9.2231 |
| dysgu    | 9010 |  508 |  631 |      0.9466 |   0.9346 | 0.9406 |           0.9339 |  9.0572 |
| NGSEP    | 8632 |  564 | 1009 |      0.9387 |   0.8953 | 0.9165 |           0.9149 |  7.6135 |
| NGSEP    | 8644 |  560 |  997 |      0.9392 |   0.8966 | 0.9174 |           0.9211 |  8.0543 |
| NGSEP    | 8706 |  588 |  935 |      0.9367 |   0.903  | 0.9196 |           0.925  |  8.4205 |
| NGSEP    | 8827 |  579 |  814 |      0.9384 |   0.9156 | 0.9269 |           0.9291 |  8.8387 |
| NGSEP    | 8837 |  582 |  804 |      0.9382 |   0.9166 | 0.9273 |           0.9338 |  9.2231 |
| NGSEP    | 8809 |  582 |  832 |      0.938  |   0.9137 | 0.9257 |           0.9306 |  9.0572 |
| svim     | 8754 |  678 |  887 |      0.9281 |   0.908  | 0.9179 |           0.8133 |  7.6135 |
| svim     | 8754 |  720 |  887 |      0.924  |   0.908  | 0.9159 |           0.8284 |  8.0543 |
| svim     | 8784 |  720 |  857 |      0.9242 |   0.9111 | 0.9176 |           0.834  |  8.4205 |
| svim     | 8925 |  721 |  716 |      0.9253 |   0.9257 | 0.9255 |           0.8491 |  8.8387 |
| svim     | 8930 |  773 |  711 |      0.9203 |   0.9263 | 0.9233 |           0.8615 |  9.2231 |
| svim     | 8883 |  745 |  758 |      0.9226 |   0.9214 | 0.922  |           0.8584 |  9.0572 |
