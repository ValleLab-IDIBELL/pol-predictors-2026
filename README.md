# pol-predictors-2026

Analysis code for the benchmarking of in-silico pathogenicity predictors in the
exonuclease domains (ED) of **POLE** and **POLD1**.

> **AlphaMissense shows improved performance for predicting the pathogenicity of POLE and
> POLD1 exonuclease domain variants: A preliminary benchmarking study with implications
> for cancer variant classification.**
>
> Raúl Marín\*, Julen Viana-Errasti\*, Pilar Mur, Hortensia Rivera, Sean V. Tavtigian,
> Laura Valle.
>
> *Genes & Diseases*, 2026 (in revision).
>
> \*These authors contributed equally.

This repository contains the single script that reproduces every figure and table of the
manuscript from the source variant table, together with the environment needed to run it.

---

## Repository structure

```
pol-predictors-2026/
├── data/
│   └── POLE_POLD1_predictors.xlsx      # input variant table
├── results/                            # generated outputs — NOT tracked
├── POLE_POLD1_predictor_evaluation.R   # the analysis
├── environment.yml                     # conda environment
├── CITATION.cff                        # how to cite this repository
├── LICENSE
└── README.md
```

`results/` is created by the script on first run and is git-ignored: every file in it is
reproducible from `data/` and the script.

---

## Input data

`data/POLE_POLD1_predictors.xlsx` holds one row per variant, in three sheets:

| Sheet | Contents |
|---|---|
| `controls` | Pathogenic and Benign ED variants used as the labelled control set |
| `VUS` | Variants of uncertain significance in the same domains |
| `gnomAD` | gnomAD ED missense variants |

Columns read by the script:

- `REVEL`, `AM`, `BD`, `CADD`, `PAI`, `PAI3D` — predictor scores on their native scales.
  Missing values (`NA`, empty, `n.a.`) are handled per predictor, so a variant lacking one
  score still contributes to the others.
- `Categorization` — `Benign` / `Pathogenic` / `VUS`; the evidence-independent label.
- `Classification` — the specific term behind each `Categorization` (`B`, `LB`, `LP`, `P`,
  `hotVUS`, `LP - somatic`, `P - somatic`). Used only to colour Supplementary Figure 1.
- `GENE` — `POLE` or `POLD1`.

The `gnomAD` sheet needs only the predictor columns and `GENE`.

---

## Requirements

R, with:

| Package | Used for |
|---|---|
| `tidyverse` | data handling and all `ggplot2` figures |
| `readxl` | reading the input workbook |
| `pROC` | ROC curves, AUC, paired DeLong test |
| `boot` | stratified BCa bootstrap CIs for the AUC |
| `Hmisc` | `binconf()` — exact Clopper–Pearson binomial CIs |
| `openxlsx` | writing the result workbooks |

Exact versions are pinned in `environment.yml`. The script prints `sessionInfo()` when it
finishes, so the run's own versions are always recorded alongside the results.

---

## Reproducing the analysis

**1. Clone the repository**

```bash
git clone https://github.com/ValleLab-IDIBELL/pol-predictors-2026.git
cd pol-predictors-2026
```

**2. Recreate the environment**

```bash
conda env create --name pol_predictors --file environment.yml
conda activate pol_predictors
```

**3. Run the script from the repository root**

```bash
Rscript POLE_POLD1_predictor_evaluation.R
```

Paths inside the script are relative to the repository root, so run it from there. All
figures and tables are written to `results/`, and a full report is printed to the console.

---

## Outputs

All written to `results/`.

### Figures (PDF)

| File | Manuscript |
|---|---|
| `Figure1A_predictor_scores_B-P.pdf` | Figure 1A — predictor scores in the P/B controls |
| `Figure1B_ROC_curves.pdf` | Figure 1B — ROC curve of every predictor |
| `Figure1E_VUS_AM_thresholds.pdf` | Figure 1E — VUS AM scores vs the AM ladder |
| `Figure1F_VUS_REVEL_thresholds.pdf` | Figure 1F — VUS REVEL scores vs the REVEL ladder |
| `Figure1G_gnomAD_AM_thresholds.pdf` | Figure 1G — gnomAD AM scores vs the AM ladder |
| `Figure1H_gnomAD_REVEL_thresholds.pdf` | Figure 1H — gnomAD REVEL scores vs the REVEL ladder |
| `SupplementaryFigure1_predictor_scores_B-P_classification.pdf` | Supplementary Figure 1 — controls, coloured by classification subgroup |
| `SupplementaryFigure2_predictor_scores_B-P_gene.pdf` | Supplementary Figure 2 — controls, coloured by gene |
| `SupplementaryFigure3_predictor_scores_VUS.pdf` | Supplementary Figure 3 — predictor scores in the VUS |
| `SupplementaryFigure4_predictor_scores_gnomAD.pdf` | Supplementary Figure 4 — predictor scores in gnomAD |

Figures 1C and 1D are drawn from the result tables below.

### Tables (XLSX)

| File | Sheet | Contents |
|---|---|---|
| `POLE_POLD1_AUC_results.xlsx` | `AUC_results` | AUC, 95% BCa CI and permutation *p* per predictor (Figure 1C) |
| | `pairwiseDeLong_p` | Pairwise paired DeLong *p*-value matrix (Figure 1C) |
| `POLE_POLD1_evidence_rules.xlsx` | `controls` | Sensitivity, specificity, PPV, NPV, yield and conflicts of the AM / REVEL / Hybrid rules on the controls (Figure 1D) |
| `POLE_POLD1_tier_distribution.xlsx` | `SupplTable4_AM` | Supplementary Table 4 — controls, VUS and gnomAD counts per AM evidence interval |
| | `SupplTable5_REVEL` | Supplementary Table 5 — the same for REVEL |

---

## Citation

If you use this code, please cite the paper:

> Marín R\*, Viana-Errasti J\*, Mur P, Rivera H, Tavtigian SV, Valle L. AlphaMissense shows
> improved performance for predicting the pathogenicity of POLE and POLD1 exonuclease
> domain variants: A preliminary benchmarking study with implications for cancer variant
> classification. *Genes & Diseases*. 2026. (in revision)

The DOI will be added here once the paper is published. `CITATION.cff` in this repository
carries the same information in machine-readable form, and GitHub's *Cite this repository*
button reads it directly.

---

## License

Released under the MIT License — see [`LICENSE`](LICENSE). The licence covers the analysis code;
the variant data in `data/` derives from public sources and from the cohorts described in the manuscript, and remains subject to the terms of those sources.

---

## Contact

- **Laura Valle** — [lvalle@idibell.cat](mailto:lvalle@idibell.cat) — Hereditary Cancer
  Program, Catalan Institute of Oncology (ICO), Hereditary Cancer Group, Bellvitge
  Biomedical Research Institute (IDIBELL), L'Hospitalet de Llobregat, Barcelona, Spain.
- Questions about the code or problems reproducing the analysis: please open an
  [issue](https://github.com/ValleLab-IDIBELL/pol-predictors-2026/issues).
