hisat2 -p 4 -x ./hisat2_grch38/grch38/genome \
    -1 subamostras_5M/SRR37459106_1.fastq.gz \
    -2 subamostras_5M/SRR37459106_2.fastq.gz \
    | samtools sort -@ 4 -o aligned/SRR37459106.sorted.bam

hisat2 -p 4 -x ./hisat2_grch38/grch38/genome \
    -1 subamostras_5M/SRR37459107_1.fastq.gz \
    -2 subamostras_5M/SRR37459107_2.fastq.gz \
    | samtools sort -@ 4 -o aligned/SRR37459107.sorted.bam

hisat2 -p 4 -x ./hisat2_grch38/grch38/genome \
    -1 subamostras_5M/SRR37459108_1.fastq.gz \
    -2 subamostras_5M/SRR37459108_2.fastq.gz \
    | samtools sort -@ 4 -o aligned/SRR37459108.sorted.bam

hisat2 -p 4 -x ./hisat2_grch38/grch38/genome \
    -1 subamostras_5M/SRR37459109_1.fastq.gz \
    -2 subamostras_5M/SRR37459109_2.fastq.gz \
    | samtools sort -@ 4 -o aligned/SRR37459109.sorted.bam

hisat2 -p 4 -x ./hisat2_grch38/grch38/genome \
    -1 subamostras_5M/SRR37459110_1.fastq.gz \
    -2 subamostras_5M/SRR37459110_2.fastq.gz \
    | samtools sort -@ 4 -o aligned/SRR37459110.sorted.bam

hisat2 -p 4 -x ./hisat2_grch38/grch38/genome \
    -1 subamostras_5M/SRR37459111_1.fastq.gz \
    -2 subamostras_5M/SRR37459111_2.fastq.gz \
    | samtools sort -@ 4 -o aligned/SRR37459111.sorted.bam

samtools index aligned/SRR37459106.sorted.bam
samtools index aligned/SRR37459107.sorted.bam
samtools index aligned/SRR37459108.sorted.bam
samtools index aligned/SRR37459109.sorted.bam
samtools index aligned/SRR37459110.sorted.bam
samtools index aligned/SRR37459111.sorted.bam

samtools flagstat aligned/SRR37459106.sorted.bam > aligned/SRR37459106.flagstat.txt
samtools flagstat aligned/SRR37459107.sorted.bam > aligned/SRR37459107.flagstat.txt
samtools flagstat aligned/SRR37459108.sorted.bam > aligned/SRR37459108.flagstat.txt
samtools flagstat aligned/SRR37459109.sorted.bam > aligned/SRR37459109.flagstat.txt
samtools flagstat aligned/SRR37459110.sorted.bam > aligned/SRR37459110.flagstat.txt
samtools flagstat aligned/SRR37459111.sorted.bam > aligned/SRR37459111.flagstat.txt
