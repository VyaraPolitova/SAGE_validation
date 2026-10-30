# SAGE-validation

Analysis code for evaluating the SAGE cognitive screening tool: descriptive statistics, convergent validity, known-groups validity, internal consistency, and regression models.

## Data

The raw data is not included due to privacy policy.

To run the scripts with your own access to the data, place the file in a `data/` folder at the project root:

    data/SAGE_MMSE_MoCA-data.csv

Open `Thesis.Rproj` in RStudio before running any script, so the working directory is set to the project root and the path resolves correctly.

## Contents

The repository contains one R script per analysis. The scripts are independent and can be run in any order.

| File | Analysis |
|------|----------|
| `Descriptives_tables_type1.R` | Descriptive statistics |
| `Convergent_validity.R` | Convergent validity |
| `Known_groups_validity.R` | Known-groups validity |
| `Internal_consistency.R` | Internal consistency |
| `Regressions.R` | Regression models |
| `Thesis.Rproj` | RStudio project file |

## Requirements

R packages used:

`tidyverse`, `janitor`, `gtsummary`, `gt`, `psych`, `moments`, `Hmisc`, `lavaan`, `semTools`, `semPlot`, `patchwork`, `corrplot`, `conflicted`, `cowplot`, `e1071`, `car`, `effectsize`, `pwr`, `mvnormtest`, `lmtest`, `sandwich`, `biotools`, `boot`, `AER`, `MVN`, `MASS`

Install any missing packages with:

```r
install.packages(c("tidyverse", "janitor", "gtsummary", "gt", "psych", "moments",
                   "Hmisc", "lavaan", "semTools", "semPlot", "patchwork", "corrplot",
                   "conflicted", "cowplot", "e1071", "car", "effectsize", "pwr",
                   "mvnormtest", "lmtest", "sandwich", "biotools", "boot", "AER",
                   "MVN", "MASS"))
```

## How to run

Open `Thesis.Rproj` in RStudio, then run any script on its own. Each one loads the packages it needs.

## License

The code is released under the MIT License. See the `LICENSE` file.

## Citation

If you use this code, please cite:

Politova, V. (2026). *Adapting Self-Administered Gerocognitive Examination
(SAGE) for use as a screening tool for cognitive impairment
and dementia within Bulgarian context.* [Bachelor's thesis, New Bulgarian University].
