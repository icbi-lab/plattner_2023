# plattner_2022
Source code to reproduce figures from Plattner et al. (Manuscript submitted)


## Obtain the data
All preprocessed data required to run these analyses are available from [Zenodo](https://doi.org/10.5281/zenodo.7015015).

```bash
cd ../../../data
wget https://zenodo.org/record/7015015/files/clinical_data_organoids.tsv?download=1
wget -O exomeseq.zip https://zenodo.org/record/7015015/files/exomeseq.zip?download=1
wget -O imaging.zip https://zenodo.org/record/7015015/files/imaging.zip?download=1
wget -O phosphoproteomics.zip  https://zenodo.org/record/7015015/files/phosphoproteomics.zip?download=1
wget -O proteomics.zip https://zenodo.org/record/7015015/files/proteomics.zip?download=1
wget -O rnaseq.zip https://zenodo.org/record/7015015/files/rnaseq.zip?download=1
wget -O scrnaseq.zip https://zenodo.org/record/1015015/files/scrnaseq.zip?download=1
unzip *.zip
```

## Contact
If you have questions regarding the analyses, please use the [issue tracker](https://github.com/icbi-lab/plattner_2022/issues).
