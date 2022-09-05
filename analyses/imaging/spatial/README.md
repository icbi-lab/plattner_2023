# Spatial analyses of IMC

All software dependencies to run the IMC analysis python script are declared in `IMCanalysis.yml`.
You can create a conda environment from that file using

```bash
conda env create -f IMCanalysis.yml
```

## Obtain the data
All preprocessed data required to run these analyses are available from [Zenodo](https://doi.org/10.5281/zenodo.7015015).

```bash
cd ../../../data
wget -O imaging.zip https://zenodo.org/record/7015015/files/imaging.zip?download=1
unzip imaging.zip
cd analyses/imaging/spatial
```

## Run the analyses
To run the IMC voronoi analysis from the `analyses` run:

```bash
conda activate IMCanalysis

mkdir -p ../../../results/IMCspatial

./spatialAnalysis_voronoi.py \
    --result_dir ../../../results/IMCspatial \
    --data_dir ../../../data/imaging/Single_cell_data \
    --num_cores 4 \
    --legend True \
    --skip_phenotype Unknown CD57+_cells CD38+_cells


./plotVoronoi.py \
    --result_dir ../../../results/IMCspatial \
    --data_dir ../../../data/imaging/Single_cell_data

./cohort_zscore_heatmap.py  \
    --result_dir ../../../results/IMCspatial \
    --data_dir ../../../results/IMCspatial \
    --s_count \
    --direction attract \
    --discrete_cm \
    --skip_phenotype Unknown CD57+_cells CD38+_cells

./microaggregates_PD1_PDL1.py \
    --result_dir ../../../results/IMCspatial \
    --data_dir ../../../data/imaging/PD1_PDL1 \
    --num_cores 4

```


## Contact
If you have questions regarding the analyses, please use the [issue tracker](https://github.com/icbi-lab/plattner_2022/issues).

