## African Demography Review
This repository provides a Snakemake workflow used for the data analysis presented in Pfennig et al. 2023, *Evolutionary genetics and admixture in African populations*, under review at *GBE*.

#### Datasets

We combine and harmonize genomic datasets from"<br> 
- 20 European, Middle Eastern, and African populations from the Simons Genome Diversity Project (SGDP)<sup>1</sup>
- 69 African populations from Crawford et al. (2017)<sup>2</sup> and Scheinfeldt et al. (2019)<sup>3</sup> (These datasets were merged by Matt Hansen (unpublished work) and kindly provided to us)
- 17 Sudanese populations from Hollfelder et al. (2017)<sup>4</sup>
- 8 North African populations from Arauna et al. (2017)<sup>5</sup>
- 14 Sahelian populations from Fortes-Lima et al. (2022)<sup>6</sup>

Populations with more than 2 samples are downsampled to 2 indivdiuals to avoid effects from variable sample sizes. If a population was included in more than one dataset, only samples from one dataset were included for this population. We merge the dataset on biallelic common SNPs. Conflicting SNPs are strand flipped once and SNPs that remain in conflict are excluded. A minor allele frequency threshold of 0.05 and genotyping rate per SNP threshold of 0.1 and  LD pruning (max r<sup>2</sup>=0.05) were applied (these parameters can be changed in the `config/config.yaml`. Ancestry proportions were estyimated using ADMIXTURE<sup>7</sup> and different values of K. For visualization, ADMIXTURE results are aligned across all Ks using CLUMPAK<sup>9</sup> prior to plotting (Figure 2). Ancestry proportions of the best K were then interpolated in space using the ordinary Kriging method and effective migration rates were estimated using FEEMS<sup>8</sup> (Figure 3).

#### Software dependencies
- bcftools v1.13 (path can be updated in `config/config.yaml`)
- All other software is stored in provided conda environment, which are automatically created and used by Snakemake.

#### Data dependencies
- The dataset from Crawford et al.<sup>2</sup> (dbGaP: [phs001396.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001396.v1.p1)) and Scheinfeldt et al.<sup>3</sup> (dbGaP: [phs001780.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001780.v1.p1)) need to be manually downloaded and merged and put into the `data` directory. The plink prefix of the files can be updated in the `config/config.yaml`. We provide coordinates of populations used in our analysis in `data/coordinates_crawford_et_al_and_scheinfeldt_et_al.txt`.
- The dataset from Fortes-Lima et al.<sup>5</sup> needs to be downloaded separately from [https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12243?query=E-MTAB-12243](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12243?query=E-MTAB-12243) (the dataset was kindly made available by the authors prior to publication...). Population coordinates were extracted from Table S2 in the corresponding publication and are stored in `data/coordinates_fortes_lima_et_al.txt`. The plink prefix can be updated in the `config/config.yaml`.
- Geographic coordinates of the 18 Sudanese populations from Hollfelder et al. (2017)<sup>3</sup> were kindly made available by the original authors. We provide the coordinates in `data/coordinates_hollfelder_et_al.txt`.
- If city information was available for populations in Arauna et al. (2017)<sup>5</sup>, the coordinates of the city were used. For the populations from Libya and Egypt, the center of the respective country was taken (see `data/coordinates_arauna_et_al.txt`). The data file was kindly made available by the origina authors. 
- All other required data is automatically retrieved

#### Comments
The Kriging plot of admixture results currently only works for K<=4, as we observed the lowest CV error at K=4. To get the best CV error across different Ks run `grep CV {results_path}/{prefix}.*.Q`.

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
