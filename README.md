# analysis-mimicIII-causalforest

This repository contains the R code for the analysis of outcomes using MIMIC-III discharge summaries.

## Contents

- `XX.R`: main analysis script containing data processing, cleaning, and outcome modelling.

## Data

The analysis uses data from the [MIMIC-III Clinical Database](https://physionet.org/content/mimiciii/1.4/).  
Due to data usage restrictions, **we cannot share the datasets directly**.  
To access MIMIC-III data:

1. Complete the required data use agreement and credentialing at [PhysioNet](https://physionet.org).
2. Download the database (provided as a compressed folder).

### Required tables

For running this analysis, you will need at minimum:

- `NOTEEVENTS.csv`  
- `DIAGNOSES_ICD.csv`

These files should be placed in your working directory before running the script.

## Citation

If you use this code or data pipeline, please cite:

- [Johnson et al., MIMIC-III Critical Care Database, 2016](https://doi.org/10.1038/sdata.2016.35)

## License

This repository contains only code and is released under the MIT License.  
Data remain under the original MIMIC-III license and cannot be redistributed.

