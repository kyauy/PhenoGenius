# PhenoGenius

Symptom interaction modeling for precision medicine

<img src="https://raw.githubusercontent.com/kyauy/PhenoGenius/main/data/img/phenogenius.png" width="200" height="200">

## Overview

Symptom interaction model provide a method to standardize clinical descriptions and fully exploit phenotypic data in precision medicine.

This repository contains scripts and files to use PhenoGenius, the phenotype matching system for genetic disease based on these model. **Please try PhenoGenius in the cloud at [https://tinyurl.com/phenogenius-app](https://tinyurl.com/phenogenius-app).**

If you use PhenoGenius, please cite:
> Yauy et al., Learning phenotypic patterns in genetic disease by symptom interaction modeling. medrXiv (2022). [https://doi.org/10.1101/2022.07.29.22278181](https://doi.org/10.1101/2022.07.29.22278181)

## Use streamlit webapp in your desktop

Requirements

```bash
python > 3.8
poetry
```

### Install dependencies

```bash
poetry install
```

### Run

```bash
poetry shell
streamlit run phenogenius_app.py
```

Enjoy !

## License

*PhenoGenius* is licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE) for the full license text.

## Misc

*PhenoGenius* is a collaboration of :

[![SeqOne](data/img/logo-seqone.png)](https://seqone.com/)

[![Université Grenoble Alpes](data/img/logo-uga.png)](https://iab.univ-grenoble-alpes.fr/)

<a href="https://www.chu-grenoble.fr/content/service-de-genetique-genomique-et-procreation"><img src="https://www.chu-grenoble.fr/sites/all/themes/acti_main/tpl/img/logo.png"></a>
