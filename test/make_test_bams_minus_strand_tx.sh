################################################################################
## R1 sense library
## Make tests from some forward reads around the myc locus in mm10
## reads from TimelapseSeq paper, clontech pico v2
################################################################################

samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.8034746' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/insertion_forward.r1s.brd2.bam


samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.8839025' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/deletion_forward.r1s.brd2.bam

samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.17081434' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/intron_forward.r1s.brd2.bam

samtools index test/insertion_forward.r1s.brd2.bam
samtools index test/deletion_forward.r1s.brd2.bam
samtools index test/intron_forward.r1s.brd2.bam


################################################################################
## R1 sense library
## Make tests from some reverse reads around the brd2 locus in mm10
## reads from TimelapseSeq paper, clontech pico v2
################################################################################

samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.11222972' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/insertion_rev.r1s.brd2.bam


samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.5118896' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/deletion_rev.r1s.brd2.bam

samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.fr.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.13082535' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/intron_rev.r1s.brd2.bam

samtools index test/insertion_rev.r1s.brd2.bam
samtools index test/deletion_rev.r1s.brd2.bam
samtools index test/intron_rev.r1s.brd2.bam



################################################################################
## R1 antisense library
## Make tests from some forward reads around the myc locus in mm10
## reads from TimelapseSeq paper, clontech pico v2
################################################################################

samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.8034746' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/insertion_forward.r1as.brd2.bam


samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.8839025' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/deletion_forward.r1as.brd2.bam

samtools view -h -F 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.17081434' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/intron_forward.r1as.brd2.bam

samtools index test/insertion_forward.r1as.brd2.bam
samtools index test/deletion_forward.r1as.brd2.bam
samtools index test/intron_forward.r1as.brd2.bam


################################################################################
## R1 antisense library
## Make tests from some reverse reads around the brd2 locus in mm10
## reads from TimelapseSeq paper, clontech pico v2
################################################################################

samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.11222972' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/insertion_rev.r1as.brd2.bam


samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.5118896' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/deletion_rev.r1as.brd2.bam

samtools view -h -f 0x10 /Users/christopherott/Desktop/rerun_sharon_fastq/timelapse.rf.srt.bam chr17:34111690-34123092 | \
  grep -e 'SRR6262104\.13082535' -e '\@SQ' -e 'VN' | \
  samtools sort | \
  samtools view -b > ~/Documents/projects/tcq/test/intron_rev.r1as.brd2.bam

samtools index test/insertion_rev.r1as.brd2.bam
samtools index test/deletion_rev.r1as.brd2.bam
samtools index test/intron_rev.r1as.brd2.bam
