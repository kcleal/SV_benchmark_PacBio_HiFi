#!/bin/bash

# Configure job scheduler as needed

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

# run 'bash joball.sh' if not using slurm
sbatch joball.sh
