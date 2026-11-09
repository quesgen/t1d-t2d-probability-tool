# Prototype tool to estimate the probability of type 1 diabetes in adults with new-onset diabetes

R Shiny source code for the prototype web calculator described in:

> *Prototype Tool to Estimate the Probability of Type 1 Diabetes in Adults with New-Onset Diabetes
> Across Four World Regions.* Manuscript under review at *BMJ Open Diabetes Research & Care*
> (bmjdrc-2025-005904). Full citation to be added on publication.

**Live tool:** https://www.diabetestechnology.org/t1dvt2d

> ⚠️ **Research prototype.** The tool has not been externally validated and is not intended for
> clinical decision-making.

## What the tool does

For an adult with new-onset diabetes, the tool:

1. Estimates a **pre-laboratory probability** of type 1 diabetes (T1D) versus type 2 diabetes (T2D)
   from a region-specific logistic model using age at diagnosis, BMI, sex and parental history of
   diabetes. Regions: Western Europe, South Asia, Eastern Europe and Northern Europe.
2. **Updates** that probability with likelihood ratios for islet autoantibodies (GADA, IA-2A,
   ZnT8A, IAA) and fasting C-peptide, using the log-odds form of Bayes's rule. The C-peptide
   likelihood ratio is the ratio of gamma densities fitted for T1D and T2D.
3. Applies two expert-panel **rule-in criteria** (any two positive autoantibodies, or diabetic
   ketoacidosis with at least one positive autoantibody), reporting a probability of T1D >99%
   alongside the calculated value.

The formulas are given in the Supplementary Material of the paper. Combining likelihood ratios
assumes the tests are conditionally independent given diabetes type; see the paper's Limitations.

## Files

| File | Contents |
|---|---|
| `app.R` | The Shiny application. The region-specific intercepts and coefficients for log(age) and log(BMI) are in the `region_coef` list near the top. |
| `Parameters.R` | Autoantibody likelihood ratios, the gamma parameters for fasting C-peptide, and the pooled log-likelihood ratios for male sex and parental history of diabetes. |

This repository contains application code only. It contains **no patient data**: the model
parameters were estimated from datasets that belong to the contributing investigators and are not
publicly available (see the paper's data availability statement).

## Running locally

```r
install.packages(c("shiny", "ggplot2", "plotly", "viridis", "dplyr", "scales"))
shiny::runApp()
```

## Citation

Please cite the paper (citation to follow on publication) and, for the code, this repository (see
`CITATION.cff`).

## License

MIT — see `LICENSE`. Copyright © 2026 Diabetes Technology Society.
