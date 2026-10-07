# Urban-coral-climate-resilience

Scripts, data, outputs, and figures/tables associated with a laboratory experiment conducted at UM-CIMAS/NOAA-AOML (https://www.aoml.noaa.gov/experimental-reef-lab/) to determine whether 'urban corals' found on artificial substrates in the Port of Miami are more resilient to thermal and OA stress compared to their natural reef counterparts.

### Michael Studivan -- <studivanms@gmail.com>
### version: October 8, 2026

---

## Experiment

Fragments of *Orbicella faveolata* and *Siderastrea siderea* from two urban sites (Star Island, MacArthur North) and two reef sites (Emerald Reef, Rainbow Reef) were exposed to four treatments in a factorial design: contemporary or acidified pH, at ambient temperature or under a bleaching (heat) ramp. Several fragments were cut from each parent colony, so parent colony is the unit of replication for site and habitat comparisons throughout.

Treatment codes used in the scripts and outputs:

| Code | Treatment |
|---|---|
| CC | Contemporary pH + ambient temperature (control) |
| LC | Acidified + ambient temperature (OA) |
| CH | Contemporary pH + bleaching |
| LH | Acidified + bleaching |

## Repository layout

| Folder | Contents |
|---|---|
| `code/` | R Markdown scripts (`.Rmd`) and their knitted reports (`.html`), one subfolder per analysis |
| `data/` | Input data, one subfolder per analysis, plus `genotype lookup table.csv` |
| `outputs/` | Everything written by the scripts (statistics tables, figures, intermediate files) |
| `tables/` | Supplementary tables (Tables S1–S11) and datasets (Datasets S1–S4) |
| `figures/` | Manuscript figures (Figures 1–5, S1–S13) |

## Analyses

| Analysis | Script (`code/`) | Data (`data/`) |
|---|---|---|
| Tank conditions (pH, temperature) | `treatments/urban tank treatments.Rmd` | `treatments/` |
| Nutrients | `calcification/urban nutrients.Rmd` | `calcification/` |
| Survivorship | `survivorship/urban survivorship.Rmd` | `survivorship/` |
| Photochemical efficiency (Fv/Fm) | `photosynthesis/urban ipam.Rmd` | `photosynthesis/` |
| Growth (buoyant weight) | `growth/urban buoyant weight.Rmd` | `growth/` |
| Calcification and oxygen flux (incubations) | `calcification/urban carbonate.Rmd` | `calcification/` |
| Differential gene expression | `transcriptomics/urban limma ofav.Rmd`, `urban limma ssid.Rmd` | `transcriptomics/ofav/`, `transcriptomics/ssid/` |
| KOG enrichment | `transcriptomics/urban kogmwu ofav.Rmd`, `urban kogmwu ssid.Rmd` | outputs of the limma scripts |
| Gene co-expression networks (WGCNA) | `transcriptomics/urban wgcna ofav.Rmd`, `urban wgcna ssid.Rmd` | outputs of the limma scripts, plus the physiology data |
| Cross-species ortholog comparison | `transcriptomics/urban orthofinder.Rmd` | outputs of both limma scripts, `transcriptomics/orthofinder/` |

Fragments excluded from analysis are kept in the data files, with the reason recorded there; the exclusions are applied in the scripts.

## Statistical approach

The same model structure is used for the physiology and the gene expression data: site (four levels) and treatment as fixed effects, with parent colony as a random effect. Urban vs reef is tested as a planned contrast of the site model, (Star Island + MacArthur North)/2 − (Emerald + Rainbow)/2, alongside pairwise site comparisons.

- **Physiology:** linear mixed models (`lme4`/`lmerTest`, contrasts with `emmeans`), Cox models for survivorship, and generalised additive mixed models (`mgcv`) for Fv/Fm.
- **Gene expression:** a mixed model for each gene, fitted with `dream` (`variancePartition`) on voom-weighted counts, with Kenward-Roger degrees of freedom. Overall expression differences are tested by PERMANOVA (`vegan`).
- **Sensitivity analyses** in the differential expression scripts compare the primary model with DESeq2, with limma using `duplicateCorrelation`, and with Satterthwaite degrees of freedom. The *O. faveolata* script also checks the effect of pooling second-cohort samples and of individual heat-stress samples.

## Reproducing the analyses

Each script is knitted from its own folder and writes to the matching folder in `outputs/`. The physiology scripts can be knitted in any order. The transcriptomics scripts depend on each other and should be knitted in this order:

1. `urban limma ofav.Rmd` and `urban limma ssid.Rmd`
2. `urban kogmwu ofav.Rmd` and `urban kogmwu ssid.Rmd`
3. `urban wgcna ofav.Rmd` and `urban wgcna ssid.Rmd`
4. `urban orthofinder.Rmd`

Notes:

- **Run time.** The differential expression scripts are slow, because the mixed model is fitted for about 20,000 genes with Kenward-Roger degrees of freedom. Each has a switch in its setup chunk, `run_sensitivity`. Set to `TRUE` (as committed, matching the knitted reports), the sensitivity analyses run as part of the knit, which takes several hours for *S. siderea* and about ten hours for *O. faveolata*. Set it to `FALSE` for a routine knit of the primary analysis only.
- **Annotation files.** The transcriptomics scripts read gene name and KOG annotations from two separate repositories, which must sit beside this one (same parent folder): `Orbicella-faveolata-annotated-transcriptome` and `Siderastrea-siderea-annotated-transcriptome`.
- **Intermediate files.** Large intermediate `.RData` files are not tracked (see `.gitignore`); they are recreated by knitting the scripts in the order above.
- **Folder names.** The `limma` script and output folder names refer to the limma/voom framework that `dream` builds on.

## Main R packages

`tidyverse`, `lme4`, `lmerTest`, `emmeans`, `MuMIn`, `mgcv`, `gratia`, `survival`, `seacarb`, `vegan`, `DESeq2` (count normalisation and variance-stabilising transformation), `edgeR`, `limma`, `variancePartition`, `arrayQualityMetrics`, `WGCNA`, `KOGMWU`, `pheatmap`, `VennDiagram`, `ggpubr`, `cowplot`.
