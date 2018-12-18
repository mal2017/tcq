cp /Volumes/ottlab/resources/index/hisat2/GRCm38/mm10* ./

mkdir idx

mv mm10* idx/


## R1 SENSE

hisat2 -x idx/mm10 -1 SRR6262104_1.fastq.gz -2 SRR6262104_2.fastq.gz --add-chrname --rna-strandness FR --mp 4,2 -p 8 | \
	samtools view -b -h > timelapse.fr.bam



samtools sort -@ 4 -o timelapse.fr.srt.bam timelapse.fr.bam


samtools index timelapse.fr.srt.bam

## R1 ANTISENSE (simulated)

hisat2 -x idx/mm10 -1 SRR6262104_2.fastq.gz -2 SRR6262104_1.fastq.gz --add-chrname --rna-strandness RF --mp 4,2 -p 12 | \
	samtools view -b -h > timelapse.rf.bam

samtools sort -@ 4 -o timelapse.rf.srt.bam timelapse.rf.bam

samtools index timelapse.rf.srt.bam