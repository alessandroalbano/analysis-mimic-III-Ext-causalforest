# analysis-mimic-III-Ext-causalforest

This repository contains the R code for the analysis of outcomes using MIMIC-III discharge summaries.

## Contents

- `code-ASBMI-Albano-DiMaria-Sciandra-Plaia.R`: main analysis script containing data processing, cleaning, and outcome modelling.

## Data

The analysis uses data from the [MIMIC-III Clinical Database](https://physionet.org/content/mimiciii/1.4/).  
Due to data usage restrictions, **we cannot share the datasets directly**. To access MIMIC-III data:

1. Complete the required data use agreement and credentialing at [PhysioNet](https://physionet.org).
2. Download the database (provided as a compressed folder).

### Required tables

For running this analysis, you will need at minimum:

- `NOTEEVENTS.csv`  
- `DIAGNOSES_ICD.csv`

These files should be placed in your working directory before running the script.

## Citation

If you use this code or data pipeline, please cite:

- Albano A, Di Maria C, Sciandra M, Plaia A.  
  *Causal forests for discovering diagnostic language in Electronic Health Records*.  
  Applied Stochastic Models in Business and Industry. Accepted for publication, DOI: https://doi.org/10.1002/asmb.70038.

- **MIMIC-III Database**:  
  Johnson AE, Pollard TJ, Shen L, et al. (2016). MIMIC-III, a freely accessible critical care database. *Scientific Data* 3, 160035. [https://doi.org/10.1038/sdata.2016.35](https://doi.org/10.1038/sdata.2016.35)


## License

This repository contains only R code.  
Data remain under the original MIMIC-III license and cannot be redistributed.

