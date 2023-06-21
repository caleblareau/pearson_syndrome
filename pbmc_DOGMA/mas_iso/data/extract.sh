samtools view PT1_rep1.scisoseq.mapped.bam | grep AACGTTATGTAGCAGG > deletion_PT1_rep1.txt
samtools view PT1_rep1.scisoseq.mapped.bam | grep AACGTTATCGTCACAG > intact1_PT1_rep1.txt
samtools view PT1_rep1.scisoseq.mapped.bam | grep CTATAGTTGTAGCAGG > intact2_PT1_rep1.txt

samtools view PT1_rep2.scisoseq.mapped.bam | grep AACGTTATGTAGCAGG > deletion_PT1_rep2.txt
samtools view PT1_rep2.scisoseq.mapped.bam | grep AACGTTATCGTCACAG > intact1_PT1_rep2.txt
samtools view PT1_rep2.scisoseq.mapped.bam | grep CTATAGTTGTAGCAGG > intact2_PT1_rep2.txt
