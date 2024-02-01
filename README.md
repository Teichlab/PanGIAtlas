# Pan-GI integration
Integration of single cell datasets across the gastrointestinal tract. 

<p align="center">

<img src="https://user-images.githubusercontent.com/77395759/226679884-3d5616db-13bc-49a2-9e73-ccd25d0c91db.png" width=50% height=50%>

</p>

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
- He et al., Genome Biol., 2020 [DOI](https://doi.org/10.1186/s13059-020-02210-0)
- Jeong et al., Clin. Cancer Res., 2021 [DOI](https://doi.org/10.1158/1078-0432.CCR-21-0792)
- Kim et al., NPJ Precis. Oncol., 2022 [DOI](https://doi.org/10.1038/s41698-022-00251-1)
- Madisoon et al., Genome Biol., 2019 [DOI](https://doi.org/10.1186/s13059-019-1906-x)
- Pagella et al., iScience, 2021 [DOI](https://doi.org/10.1016/j.isci.2021.102405)
- Chen et al., J. Dent. Res., 2022 [DOI](https://doi.org/10.1177/002203452210760)
- Costa-da-Silva et al., iScience, 2022 [DOI](https://doi.org/10.1016/j.isci.2021.103592)

Extended atlas (original count matrices used): 
- Smillie et al., Cell, 2019, [DOI](https://doi.org/10.1016/j.cell.2019.06.029)
- Kong et al., Immunity, 2023, [DOI](https://doi.org/10.1016/j.immuni.2023.01.002)
- Garrido-Trigo et al., Nat Comms, 2023 [DOI](https://doi.org/10.1038/s41467-023-40156-6)



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


