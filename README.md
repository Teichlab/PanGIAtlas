# megagut
Integration of single cell gut datasets. 

Core datasets (remapped using STAR v2.7.9a, human reference Cell Ranger 2020-A using [STARsolo pipeline v1.0](https://github.com/cellgeni/STARsolo)): 
- Williams et al., Cell, 2021, [DOI](https://doi.org/10.1016/j.cell.2021.05.013)
- Caetano et al., eLife, 2021, [DOI](https://doi.org/10.7554/eLife.62810)
- Huang et al., Cell, 2019, [DOI](https://doi.org/10.1016/j.cell.2019.10.027)
- Jaeger et al., Nat. Comm., 2021, [DOI](https://doi.org/10.1038/s41467-021-22164-6)
- Elmentaite et al., Dev. Cell, 2020, [DOI](https://doi.org/10.1016/j.devcel.2020.11.010)
- Wang et al., J Exp Med, 2020; [DOI](https://doi.org/10.1084/jem.20191130)
- Kinchen et al., Cell, 2018, [DOI](https://doi.org/10.1016/j.cell.2018.08.067)
- Parikh et al., Nature, 2019, [DOI](https://doi.org/10.1038/s41586-019-0992-y)
- James et al., Nat. Immunol., 2020, [DOI](https://doi.org/10.1038/s41590-020-0602-z)
- Lee et al., Nature Genetics, 2020, [DOI](https://doi.org/10.1038/s41588-020-0636-z)
- Uzzan et al., Nat Med, 2022, [DOI](https://doi.org/10.1038/s41591-022-01680-y)
- Che et al., Cell Discov., 2021, [DOI](https://doi.org/10.1038/s41421-021-00312-y)
- Domínguez Conde et al., Science, 2022, [DOI](https://doi.org/10.1126/science.abl5197)
- Elmentaite et al., Nature,  2021, [DOI](https://doi.org/10.1038/s41586-021-03852-1)
- Holloway et al., Cell Stem Cell, 2020, [DOI](https://doi.org/10.1016/j.stem.2020.11.008)
- Yu et al., Cell, 2021, [DOI](https://doi.org/10.1016/j.cell.2021.04.028)
- Li et al., Nat Immunol., 2019, [DOI](https://doi.org/10.1038/s41590-018-0294-9!)
- Martin et al., Cell, 2019 [DOI](https://doi.org/10.1016/j.cell.2019.08.008)
- He et al., Genome Biology, 2020 [DOI](https://doi.org/10.1186/s13059-020-02210-0)
- Jeong et al., 2021 Clin. Cancer Res. 2021 [DOI](https://doi.org/10.1158/1078-0432.CCR-21-0792)

Extended atlas (original count matrices used): 
- Nowicki-Osuch et al., Science, 2021, [DOI](https://doi.org/10.1126/science.abd1449)
- Kumar et al., Cancer Discovery, 2022, [DOI](https://doi.org/10.1158/2159-8290.CD-21-0683)
- Corridoni et al., Nature Medicine, 2020, [DOI](https://doi.org/10.1038/s41591-020-1003-4)
- Smillie et al., Cell, 2019, [DOI](https://doi.org/10.1016/j.cell.2019.06.029)
- Pelka et al., Cell, 2021, [DOI](https://doi.org/10.1016/j.cell.2021.08.003)
- Qian et al., Cell Research, 2020, [DOI](https://doi.org/10.1038/s41422-020-0355-0)
- Uhlitz et al., EMBO Molecular Med, 2021, [DOI](https://doi.org/10.15252/emmm.202114123)
- Becker et al., Nat Genet, 2022, [DOI](https://doi.org/10.1038/s41588-022-01088-x)
- Luoma et al, Cell 2021 [DOI](https://doi.org/10.1016/j.cell.2020.06.001)
- Drokhlyansky et al., Cell, 2020, [DOI](https://doi.org/10.1016/j.cell.2020.08.003)
- Cao et al., Science, 2020, [DOI](https://doi.org/10.1126/science.aba7721)
- Gao et al., Nat. Cell. Biol., 2018, [DOI](https://doi.org/10.1038/s41556-018-0105-4)
- Kong et al., Immunity, 2023, [DOI](https://doi.org/10.1016/j.immuni.2023.01.002)
- Devlin et al., Gastroenterology, 2021 [DOI](https://doi.org/10.1053/j.gastro.2020.12.030)
- Boland et al., Science Immunology, 2020 [DOI](DOI: 10.1126/sciimmunol.abb4432)
- The Tabula Sapiens Consortium, Stephen R Quake July 2021 [DOI](https://www.science.org/doi/10.1126/science.abl4896)


## Contents

- `metadata`: this folder contains metadata used for analysis (small data objects)
- `src`: this folder contains notebooks and scripts
- `thrash`: old code snippets no longer in use, kept for OCD

## Set-up conda environment

1. Create a new conda environment

```
conda create --name megagut python=3.9
conda activate megagut
```

2. Install R and R dependencies in the environment

```
conda install conda-forge::r-base==4.0.5 bioconda::bioconductor-edger==3.32.1 rpy2==3.4.2 conda-forge::r-statmod==1.4.37
```


