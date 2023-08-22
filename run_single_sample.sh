#!/bin/bash

# Map and call SVs, configure scheduler header as needed:
for i in {44..49}; do \

$(cat > job${i}.sh << EOF
#!/bin/bash
#SBATCH --job-name=job${i}
#SBATCH --partition=compute
#SBATCH --output=out${i}.%J
#SBATCH --error=err.${i}.%J
#SBATCH --time=0-6:00
#SBATCH --mem-per-cpu=6G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

minimap2 -t8 -a -x map-hifi hs37d5.fa SRR103822${i}_subreads.fastq.gz | \
  samtools view -bh - | \
  samtools sort -o SRR103822${i}.mm2.bam -
samtools index -@8 SRR103822${i}.mm2.bam

samtools coverage -r 1 SRR103822${i}.mm2.bam > SRR103822${i}.cov.tsv

sniffles --threads 8 --input SRR103822${i}.mm2.bam --vcf HG002_${i}.pacbio.sniffles.vcf

dysgu call --mode pacbio --procs 8 --clean hs37d5.fa wd_${i} SRR103822${i}.mm2.bam -o HG002_${i}.pacbio.dysgu.vcf

mkdir wd_cuteSV_${i}
cuteSV -t 8 -s 2 --genotype SRR103822${i}.mm2.bam hs37d5.fa HG002_${i}.pacbio.cuteSV.vcf wd_cuteSV_${i}

delly lr -g hs37d5.fa SRR103822${i}.mm2.bam > HG002_${i}.pacbio.delly.vcf

java -jar NGSEPcore_4.3.2.jar SingleSampleVariantsDetector -runOnlySVs -runLongReadSVs -i SRR103822${i}.mm2.bam -r hs37d5.fa -o HG002_${i}.pacbio.NGSEP.vcf

EOF
)
  # change this to 'bash job${i}' if running without slurm
  sbatch job${i}.sh
done
