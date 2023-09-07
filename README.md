# scDREAMER-experiments

# scDREAMER
## Overview
**scDREAMER** is a single-cell data integration framework that employs a novel adversarial variational autoencoder for learning lower-dimensional cellular embeddings and a batch classifier neural network for the removal of batch effects. See our preprint below for more details. 

<img src='architecture.pdf'>


## Installation

A stable `pip` installation release for scDREAMER package will be made available shortly. For now, we recommend users to directly clone our stable `main` branch and set `scDREAMER` as the working directory. Creating conda environment using ./ENVIRONMENTS/scDREAMER.yml will install all the dependent packages and libraries. scDREAMER can be set up as follows 

```
git clone https://github.com/Zafar-Lab/scDREAMER.git
cd scDREAMER/Environments
conda env create -f scDREAMER.yml
conda activate scdreamer
```
## What Computational tasks can scDREAMER be used for?

`scDREAMER` suite can be used for:
1. scDREAMER for an unsupervised integration of multiple batches
2. scDREAMER-SUP for a supervised integration across multiple batches
3. scDREAMER-SUP can also be when cell type annotations are missing in the datasets i.e., 10%, 20%, 50%
4. Atlas level and cross-species integration
5. Large datasets with ~1 million cells

## Documentation
read the docs: [https://scdreamer.readthedocs.io/en/latest/](https://scdreamer.readthedocs.io/en/latest/)

## Contributing
In case of any bug reports, enhancement requests, general questions, and other contributions, please create an issue. For more substantial contributions, please fork this repo, push your changes to your fork, and submit a pull request with a good commit message.

## Cite this article
Ajita Shree*, Musale Krushna Pavan*, Hamim Zafar. scDREAMER: atlas-level integration of single-cell datasets using deep generative model paired with adversarial classifier. bioRxiv 2022.07.12.499846; doi: https://doi.org/10.1101/2022.07.12.499846  
\* equally contributed
