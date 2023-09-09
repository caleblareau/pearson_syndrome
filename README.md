# Pearson Syndrome single-cell genomics analysis repository
_Caleb Lareau_ 

This repository contains code needed to reproduce analyses and figures associated with this work:

Lareau, _et al._ *Single-cell multi-omics of mitochondrial DNA disorders reveals dynamics of purifying selection across human immune cells*. 2023 **Nature Genetics** https://www.nature.com/articles/s41588-023-01433-8.

This table coordinates the different analysis folders (in addition to extended data figures)

```
Figure 1 | cell_line_scatac,mtDNA_deletion_sims,
Figure 2 | pbmc_scatac,pbmc_scrna,pbmc_DOGMA
Figure 3 | tcell_culture
Figure 4 | melas_kss_cpeo
Figure 5 | bmmnc_scatac,cd34_scatac,pt3_chr7_del_scatac
Figure 6 | bmmnc_asap
Figure 7 | erythroid_culture_scatac,erythroid_culture_scrna
```

## Setup
To best use this resource, we recommend pairing with large data files (that are not compatible with github as they exceed 100Mb). These files are available from the [Open Science Framework](https://osf.io/bnxjh/).

Once one downloads the zip archieve from OSF (~60 Gb), place the extracted folder in the same directory as this repository named `pearson_large_data_files` (as shown below). This will enable running custom code to reproduce items in the output folders. 

```
.
├── pearson_large_data_files
│   ├── input
│   └── output
└── pearson_syndrome (this repository)
    ├── README.md
```

## Software
The single-cell large mitochondrial DNA genotyping workflow is available in the [mgatk-del module, documented here](https://github.com/caleblareau/mgatk/wiki/Large-deletion-calling-and-heteroplasmy-estimation). 


Questions? Raise an issue on this repository

<br><br>
