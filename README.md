---
title: PhenoGenius
emoji: genie
sdk: streamlit
sdk_version: 1.25.0
app_file: phenogenius_app.py
python_version: 3.8
pinned: true
---

# PhenoGenius

Symptom interaction modeling for precision medicine

## Overview

Symptom interaction model provide a method to standardize clinical descriptions and fully exploit phenotypic data in precision medicine.

This repository contains scripts and files to use PhenoGenius, the phenotype matching system for genetic disease based on this model. **Please try PhenoGenius in the cloud at [https://huggingface.co/spaces/kyauy/PhenoGenius](https://huggingface.co/spaces/kyauy/PhenoGenius).**

If you use PhenoGenius, please cite:
> Yauy et al., Learning phenotypic patterns in genetic disease by symptom interaction modeling. medrXiv (2023). [https://doi.org/10.1101/2022.07.29.22278181](https://doi.org/10.1101/2022.07.29.22278181)

## Install

- Requirements

```bash
python == 3.8 #(pyenv install 3.8)
poetry #(https://python-poetry.org/docs/#installation)
```

- Install dependencies

```bash
poetry install
```

## Use streamlit webapp in your desktop

### Run

```bash
poetry shell
streamlit run phenogenius_app.py
```

## Use command-line client

### Run

```bash
poetry shell
python phenogenius_cli.py --hpo_list HP:0000107,HP:0000108,HP:0001407,HP:0005562 --result_file PKD1.tsv
```
```
Usage: phenogenius_cli.py [OPTIONS]

Options:
  --result_file TEXT  Output file name, default = match.tsv
  --hpo_list TEXT     (Mandatory) List of HPO terms to match, separated with
                      commas
  --gene_list TEXT    (Optional) List of genes in NCBI ID format to match,
                      separated with commas
```

## Explore interactive graphs of symptoms interactions

### Human Phenotype Ontology

Click on the image!
[![HPO](data/graph/onto_image.png)](https://ouestware.gitlab.io/retina/beta/#/graph/?url=https%3A%2F%2Fraw.githubusercontent.com%2Fkyauy%2FPhenoGenius%2Fmain%2Fdata%2Fgraph%2Fontology.gexf&r=v&n=n3453&sa=r&ca=f&st[]=n&st[]=f)

### Groups of interacting symptoms

Click on the image!
[![Groups](data/graph/group_image.png)](https://ouestware.gitlab.io/retina/beta/#/graph/?url=https%3A%2F%2Fraw.githubusercontent.com%2Fkyauy%2FPhenoGenius%2Fmain%2Fdata%2Fgraph%2F390groups.gexf&r=v&n=n16738&sa=r&ca[]=f&ca[]=l&st=f&ls=5oGenius%2Fmain%2Fdata%2Fgraph%2Fontology.gexf&r=v&n=n3453&sa=r&ca=f&st[]=n&st[]=f)

Enjoy !

## License

*PhenoGenius* is licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) for the full license text.

## Misc

*PhenoGenius* is a collaboration of :

[![SeqOne](data/img/logo-seqone.png)](https://seqone.com/)

[![Université Grenoble Alpes](data/img/logo-uga.png)](https://iab.univ-grenoble-alpes.fr/)
