## African Demography Review
This repository provides a Snakemake workflow used for the data analysis presented in Pfennig et al. 2023, *Evolutionary genetics and admixture in African populations*, under review at *GBE*.

#### Datasets

We combine and harmonize genomic datasets from"<br> 
- 20 European, Middle Eastern, and African populations from the Simons Genome Diversity Project (SGDP)<sup>1</sup>
- 69 African populations from Crawford et al.(2019)<sup>2</sup> and Scheinfeldt et al. (2019)<sup>3</sup> (These datasets were merged by Matt Hansen (unpublished work) and kindly provided to us)
- 17 Sudanese populations from Hollfelder et al. (2017)<sup>4</sup>
- 8 North African populations for which city information allowed to obtain geographic coordinates. For the populations from Libya and Egypt, the center of the respectibve country was taken (see `data/coordinates_arauna_et_al.txt`) from Arauna et al. (2017)<sup>5</sup>
- 14 Sahelian populations from Fortes-Lima et al. (2022)<sup>6</sup>

Populations with more than 2 samples are downsampled to 2 indivdiuals to avoid effects from variable sample sizes. If a population was included in more than one dataset, only samples from one dataset were included for this population. We merge the dataset on biallelic common SNPs. Conflicting SNPs are strand flipped once and SNPs that remain in conflict are excluded. A minor allele frequency threshold of 0.05 and genotyping rate per SNP threshold of 0.1 are applied (These parameters can be changed in the `config/config.yaml`. Subsequently, we apply LD pruning (max r<sup>2</sup>=0.05) before estimating ancestry proportions using ADMIXTURE<sup>7</sup> and effective migrations rates using FEEMS<sup>8</sup> (Figure 2C). For visualization, ADMIXTURE results are aligned across all Ks using CLUMPAK<sup>9</sup> prior to plotting (Figure S1). The estimated admixture proportions of the best K are shown in Figure 2A and ancestry proportions interpolated using the Kriging method in Figure 2B.

#### Software dependencies
- bcftools v1.13 (path can be updated in `config/config.yaml`)
- All other software is stored in provided conda environment, which are automatically created and used by Snakemake.

#### Data dependencies
- The dataset from Crawford et al.<sup>2</sup> needs to be downloaded manually (including ethnographic and sample information) and put into the `data` directory. The plink prefix of the files can be updated in the `config/config.yaml`.
- The dataset from Fortes-Lima et al.<sup><5</sup> needs to be downloaded separately. Population coordinates are provided in Table S2 in the corresponding publications. The plink prefix can be updated in the `config/config.yaml`.
- Geographic coordinates of the 18 Sudanese populations from Hollfelder et al. (2017)<sup>3</sup> (geographic coordinates are available upon request from the original authors). The path can be updated `config/config.yaml`.
- All other required data is automatically retrieved

#### Comments
We use K=4 as the best K for generating the admixture plot, although it was not strictly the best K. The CV error was plateuing and jumped at K=7. To get the best CV error run `grep CV {results_path}/{prefix}.*.Q`. You can then set the desired K for plotting the `config/config.yaml` file. Note that the Kriging plot of admixture results currently only works for K<=4.

As of now, this Snakemake is workflow **is not written in a general way** that allows to download/merge data from any dataset, but only those mentioned above.

#### References
1. Mallick S, Li H, Lipson M, et al. The Simons Genome Diversity Project: 300 genomes from 142 diverse populations. Nature 538, 201–206 (2016). [https://doi.org/10.1038/nature18964](https://doi.org/10.1038/nature18964)
2. Crawford, N. G., Kelly, D. E., Hansen, M. E. B., et al. (2017) Loci associated with skin pigmentation identified in African populations. Science 358, eaan8433. [10.1126/science.aan8433](10.1126/science.aan8433)
3. Scheinfeldt LB et al. 2019. Genomic evidence for shared common ancestry of East African hunting-gathering populations and insights into local adaptation. Proc. Natl. Acad. Sci. 116:4166–4175. doi: [https://www.pnas.org/doi/full/10.1073/pnas.1817678116](https://www.pnas.org/doi/full/10.1073/pnas.1817678116).
4. Hollfelder N, Schlebusch CM, Günther T, et al. (2017) Northeast African genomic variation shaped by the continuity of indigenous groups and Eurasian migrations. PLoS Genet 13(8): e1006976. [https://doi.org/10.1371/journal.pgen.1006976](https://doi.org/10.1371/journal.pgen.1006976)
5. Arauna LR, Mendoza-Revilla J, Mas-Sandoval A, et al. (2017) Recent Historical Migrations Have Shaped the Gene Pool of Arabs and Berbers in North Africa, Molecular Biology and Evolution, Volume 34, Issue 2, [https://doi.org/10.1093/molbev/msw218](https://doi.org/10.1093/molbev/msw218)
6. Fortes-Lima C, Tříska P, Čížková M, et al. (2022) Demographic and Selection Histories of Populations Across the Sahel/Savannah Belt, Molecular Biology and Evolution, Volume 39, Issue 10, msac209, [https://doi.org/10.1093/molbev/msac209](https://doi.org/10.1093/molbev/msac209)
7. Alexander DH, Novembre J, and Lange K. Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19:1655–1664, 2009.[http://www.genome.org/cgi/doi/10.1101/gr.094052.109](http://www.genome.org/cgi/doi/10.1101/gr.094052.109)
8. Marcus J, Wooseok H, Barber RF, Novembre J (2021) Fast and flexible estimation of effective migration surfaces eLife 10:e61927. [https://doi.org/10.7554/eLife.61927](https://doi.org/10.7554/eLife.61927)
9. Kopelman, N, Mayzel, J, Jakobsson, M, Rosenberg, N, Mayrose, I (2015) CLUMPAK: a program for identifying clustering modes and packaging population structure inferences across K. Molecular Ecology Resources 15(5): 1179-1191, [https://doi.org/10.1111/1755-0998.12387](https://doi.org/10.1111/1755-0998.12387)
