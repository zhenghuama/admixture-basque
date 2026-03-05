## Admixture with EM algorithm using Squared Extrapolation Methods

### Introduction

This repo implements Admixture using 2 methods

1. E-M algorithm as a naive approach
2. Squared Extrapolation Methods (SQUAREM) to accelerate E-M algorithm

There are 3 tasks to test the implementations

- Step 1: Test on Chr 21 from 1000Genomes dataset on 5 globally distributed populations (ASW, CEU, GWD, PER, PUL)

- Step 2 (TODO): Test on 1000Genomes dataset, but on less globally diverse populations

- Step 3 (TODO): Test on Human Genome Diversity Project (HGDP) dataset, with a focus on the Basque populations v.s. other European populations.

### Code structure

- `admixture_em.py` and `admixture_squarem.py`: source code of E-M algorithm and SQUAREM method
- `run_em_step1.ipynb` and `run_squarem_step1.ipynb`: Notebooks to run the 2 above methods and visualization.
- `data_and_result/step1/intermediate`: `.bed`, `.bim` and `.fam` files post LD-pruning.


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

The program for steps 1, 2 and 3 are identical. We only put the intermediate dataset and reference results of step 1 for now. (TODO: Put step 2 and step 3 datasets).

Simply open `run_em_step1.ipynb` and `run_squarem_step1.ipynb` and run all the cells to see the visualized admixture results and how it compares with the reference used in class.

Logs will be available at `data_and_results/step1/results_em/step1_admixture.3.0` and `data_and_results/step1/results_squarem/step1_admixture.3.0`. The logs contain insights e.g., convergence speed and if SQUAREM is effecive. TODO: Visualize convergence comparison & improvement by adopting SQUAREM.

### TODOs

1.  Test the programs on datasets for steps 2 and 3.
2.  Visualize convergence of E-M methods v.s. SQUAREM.


