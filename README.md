## Admixture with EM algorithm using Squared Extrapolation Methods

### Introduction

This repo implements Admixture using 2 methods

1. E-M algorithm as a naive approach
2. Squared Extrapolation Methods (SQUAREM) to accelerate E-M algorithm

There are 3 tasks to test the implementations

- Step 1: Test on Chr 21 from 1000Genomes dataset on 5 globally distributed populations (ASW, CEU, GWD, PER, PUL)

- Step 2: Test on 1000Genomes dataset, but on less globally diverse, but still distinguishable populations (FIN, IBS, GBR, CEU, TSI)

- Step 3: Test on harmonized 1000Genomes and Human Genome Diversity Project (HGDP) dataset, with Basque, Russian, FIN, IBS, GBR, IBS and TSI populations. We are interested in seeing if admixture display significant different ancestry patterns for Basque people and other Europeans, given their linguistic uniqueness.

### Code structure

- `run_step1.ipynb`, `run_step1.ipynb` and `run_step3.ipynb` are the main 3 notebooks, each for one dataset. The notebook include running the reference ADMIXTURE tools as well as EM and SQUAREM. Benchmarking results are printed by cells. The notebook also include visualization of ancestry disctribution after admixture.
- `admixture_em.py` and `admixture_squarem.py` are the source code of E-M algorithm and SQUAREM method. `plot_admixture.py` visualize the ancestry results.
- `pre-processing_stepx.sh` are pre-processing files. They would read from the `.vcf` and `.tsv` files and perform L-D pruning. You do not need to run these files as the pruned results are already stored in this repo. The original `.vcf` files however, are too big to be uploaded to GitHub.
- `reference_stepx.sh` calls the reference ADMIXTURE tool. They will be called by the main `.ipynb` files.


### Steps

#### Clone the repo

```bash
git clone https://github.com/zhenghuama/admixture-basque.git
```

#### Dependencies
- OS: Linux x86_64
- Python Packages: pandas, numpy, matplotlib, ipykernel
- Optional Packages: plink, admixture

Optional packages are required if the user wish to start from .vcf file and use plink tools to do LD-pruning and/or run the admixture tools as a reference.

We have prepared the post LD-pruning data as intermediate dataset and the admixture program output as reference results. 

We also provide the conda enviroment here.

```bash
conda env create -f environment.yml
conda activate bioinfo
```

This would take a while to resolve all dependencies. It's not necessary to use the given conda enviroment if the dependencies are already satisfied.

#### Run the code

Simply open `run_step1.ipynb`, `run_step2.ipynb`, `run_step3.ipynb` and run all the cells to sequentailly call 3 admixture methods and see the benchmarking and visualization results. The EM algorithm will take a while to finish. Make sure to select the correct python kernel for the notebook (should be conda enviroment bioinfo). The notebook has a `which python` cell to check the python version.

#### Data and log structure
Each step has its own folder under `data_and_resutls`. For each step,
- `datasets` include the `.tsv` file
- `intermediate` stores data required by admixture after L-D pruning
- `reference_results` stores admixture results by the reference ADMIXTURE tool, logs file ends with `.3.0`
- `em_results` stores admixture results by the implemented EM algorithm, logs file ends with `.3.0`
- `squarem_results` stores admixture results by the implemented SQUAREM method, logs file ends with `.3.0`
