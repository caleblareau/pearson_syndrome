
bsub -q big python3 ../pythondev/scSplit/scSplit count -v Pearson-Mix-5.HQ.vcf -i pearson/cell-line-mix/Pearson-Mix-5_v12-mtMask/outs/possorted_bam.bam -b pearson/cell-line-mix/Pearson-Mix-5_v12-mtMask/outs/filtered_peak_bc_matrix/barcodes.tsv -r ref_mix5.csv -a alt_mix5.csv -c hg19_common.txt -o pearson/demuxlet
bsub -q big python3 ../pythondev/scSplit/scSplit count -v Pearson-Mix-3.HQ.vcf -i pearson/cell-line-mix/Pearson-Mix-3_v12-mtMask/outs/possorted_bam.bam -b pearson/cell-line-mix/Pearson-Mix-3_v12-mtMask/outs/filtered_peak_bc_matrix/barcodes.tsv -r ref_mix3.csv -a alt_mix3.csv -c hg19_common.txt -o pearson/demuxlet




