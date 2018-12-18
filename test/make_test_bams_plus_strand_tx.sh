################################################################################
## R1 sense library
## Make tests from some forward reads around the myc locus in mm10
## reads from TimelapseSeq paper, clontech pico v2
################################################################################

samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.531777' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/insertion_forward.r1s.myc.bam


samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.15211333' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/deletion_forward.r1s.myc.bam

samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.12992419' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/intron_forward.r1s.myc.bam

samtools index test/insertion_forward.r1s.myc.bam
samtools index test/deletion_forward.r1s.myc.bam
samtools index test/intron_forward.r1s.myc.bam


################################################################################
## R1 sense library
## Make tests from some reverse reads around the max locus in mm10
## reads from TimelapseSeq paper, clontech pico v2
################################################################################

samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.21038055' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/insertion_rev.r1s.myc.bam


samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.17760503' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/deletion_rev.r1s.myc.bam

samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.11460943' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/intron_rev.r1s.myc.bam

samtools index test/insertion_rev.r1s.myc.bam
samtools index test/deletion_rev.r1s.myc.bam
samtools index test/intron_rev.r1s.myc.bam



################################################################################
## R1 antisense library
## Make tests from some forward reads around the myc locus in mm10
## reads from TimelapseSeq paper, clontech pico v2
################################################################################

samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.531777' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/insertion_forward.r1as.myc.bam


samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.15211333' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/deletion_forward.r1as.myc.bam

samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.12992419' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/intron_forward.r1as.myc.bam

samtools index test/insertion_forward.r1as.myc.bam
samtools index test/deletion_forward.r1as.myc.bam
samtools index test/intron_forward.r1as.myc.bam


################################################################################
## R1 antisense library
## Make tests from some reverse reads around the myc locus in mm10
## reads from TimelapseSeq paper, clontech pico v2
################################################################################

samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.21038055' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/insertion_rev.r1as.myc.bam


samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.17760503' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/deletion_rev.r1as.myc.bam

samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr15:61979978-61999209 | \
  grep -e 'SRR6262104\.11460943' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/intron_rev.r1as.myc.bam

samtools index test/insertion_rev.r1as.myc.bam
samtools index test/deletion_rev.r1as.myc.bam
samtools index test/intron_rev.r1as.myc.bam
