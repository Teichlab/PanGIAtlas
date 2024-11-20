# Pan-GI Cell Atlas
Integration of human single cell datasets across the gastrointestinal tract. 

Welcome to the Human Pan-Gastrointestinal Cell Atlas GitHub repository! Here you will find notebooks and scripts related to the [publication](https://doi.org/10.1038/s41586-024-07571-1) along with a [tutorial](https://github.com/Teichlab/PanGIAtlas/blob/main/notebooks/PanGI_atlas_example_notebook.ipynb) of how to map new data to the atlas. You can explore and download the atlas at [gutcellatlas.org](https://www.gutcellatlas.org/pangi.html). 

<p align="center">

![WebPortal_image-01](https://github.com/Teichlab/PanGIAtlas/assets/77395759/7f2d683f-572e-475a-82e6-b21f465a0d9d)

</p>

## Atlas background
For a comprehensive explanation of the atlas, please check out the [paper](https://doi.org/10.1038/s41586-024-07571-1). In brief, the Pan-GI atlas brings together 25 single-cell RNA sequencing datasets from across the gastrointestinal (GI) tract and the human lifespan. The resource encompasses a healthy reference atlas of ~1.1 million cells from GI samples of 189 healthy controls and an extended atlas with an additional 500k cells from 5 GI diseases. The data was processed uniformly using a newly developed automated quality control approach ([scAutoQC](https://github.com/Teichlab/sctk)), described in the paper. The atlas was utilised to investigate cellular changes in inflammatory intestinal diseases (coeliac disease, ulcerative colitis and Crohn’s disease) and identified metaplastic cells with inflammatory signalling circuits. Users of the atlas can investigate their own hypotheses within the existing datasets in the atlas and can map newly generated data to help annotate and contextualise with the Pan-GI atlas.

## How to use the atlas

To explore and download the atlas, please go to [gutcellatlas.org](https://www.gutcellatlas.org/pangi.html). 

**Browse:** There is a cellxgene browser on gutcellatlas.org to explore the main atlas of 1.6 million cells. Currently under construction are cellxgene browsers for each cell lineage, in order to browse the data in more detail. Please check back for updates on when they are available, likely early 2025, in the meantime you can use the cellxgene [selecting and subsetting cells function](https://cellxgene.cziscience.com/docs/04__Analyze%20Public%20Data/4_1__Hosted%20Tutorials).

**Download:** There are multiple download links on the website for the core atlas, the extended atlas and the lineage specific objects (.h5ad and .rds). There are also links for the scANVI-based reference models to map and annotate your own data to the atlas and CellTypist models to annotate new data.

**Map and annotate:** In order to map and annotate newly generated data you will need to download objects and/or models, depending on the research question and method. To map and annotate using scANVI-based reference models, please see our [tutorial](https://github.com/Teichlab/PanGIAtlas/blob/main/notebooks/PanGI_atlas_example_notebook.ipynb). To annotate using CellTypist, please see existing tutorials on the [CellTypist website](https://www.celltypist.org/tutorials). Any single-cell data can be mapped and annotated using the atlas, however for the most accurate results it is recommended to use scAutoQC (link to [github page](https://github.com/Teichlab/sctk) and [tutorial](https://teichlab.github.io/sctk/notebooks/automatic_qc.html)).

## Any questions?
Please submit an issue to this GitHub repository.


## Contents

- `notebooks`: this folder contains all notebooks
- `scripts`: this folder contains scripts

## Datasets
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



