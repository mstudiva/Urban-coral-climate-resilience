# Urban-coral-climate-resilience

Scripts, data, outputs, and figures/tables associated with a laboratory experiment conducted at UM-CIMAS/NOAA-AOML (https://www.aoml.noaa.gov/experimental-reef-lab/) to determine whether 'urban corals' found on artificial substrates in the Port of Miami are more resilient to thermal and OA stress compared to their natural reef counterparts.

### Michael Studivan -- <studivanms@gmail.com>
### version: October 26, 2026

---

## Experiment

Cores of *Orbicella faveolata* and *Siderastrea siderea* were collected in February 2021 from two urban sites in the Port of Miami (Star Island, MacArthur North) and two reef sites (Emerald Reef, Rainbow Reef): four cores from each of five parent colonies per site and species (six *O. faveolata* colonies at MacArthur North; 84 *O. faveolata* and 80 *S. siderea* cores). The experiment ran from February to June 2021 in the Experimental Reef Lab and had four treatments in a factorial design: contemporary (pH<sub>T</sub> 8.05) or acidified (7.80) pH, at ambient temperature (22.5 °C) or bleaching temperature (33 °C).

- **Acidification:** pH was ramped over March 1–7 and held for the rest of the experiment (from March 8).
- **First heat cohort:** the bleaching-temperature tanks were ramped from April 5 and reached 33 °C on April 19. Tissue for gene expression was sampled on April 23 (four days at 33 °C), while the ambient tanks were still at 22.5 °C.
- **Second heat cohort:** with no mortality at ambient temperature, the ambient tanks were ramped from April 26 and reached 33 °C on May 10, so every fragment was eventually heat-stressed. Ten acidified *O. faveolata* fragments were resampled on May 14; these are used only in a sensitivity analysis.
- **Measurements:** tank pH and temperature (loggers, weekly bottle samples), buoyant-weight growth (February 26–April 3), light/dark incubations for calcification, oxygen flux and nutrients (March 29–April 2), photochemical efficiency (Fv/Fm), survivorship under heat stress, and Tag-Seq gene expression.
- **Exclusions:** 31 *O. faveolata* fragments with partial mortality unrelated to heat stress (an apparent ciliate infection) and one *S. siderea* fragment without health observations are excluded from all analyses. They stay in the data files, with the reason recorded there, and the scripts apply the exclusions.

Several fragments were cut from each parent colony, so parent colony is the unit of replication for site and habitat comparisons throughout.

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
| Gene co-expression networks (WGCNA) | `transcriptomics/urban wgcna ofav.Rmd`, `urban wgcna ssid.Rmd` | outputs of the limma scripts, the physiology data, and `transcriptomics/symbiont alignment counts.csv` |
| Front-loading test (Figure 4c) | `transcriptomics/urban frontloading.Rmd` | outputs of both limma scripts |
| Cross-species ortholog comparison | `transcriptomics/urban orthofinder.Rmd` | outputs of both limma scripts, `transcriptomics/orthofinder/` |

Symbiont genus proportions (Table S2; WGCNA traits) come from Tag-Seq reads aligned to Symbiodiniaceae 28S rRNA references: per-sample counts are in `data/transcriptomics/symbiont alignment counts.csv` and the references in `symbiont references.csv`.

Raw Tag-Seq reads are in the NCBI Sequence Read Archive (BioProject PRJNA899003).

## Statistical approach

The same model structure is used for the physiology and the gene expression data: site (four levels) and treatment as fixed effects, with parent colony as a random effect. Urban vs reef is tested as a planned contrast of the site model, (Star Island + MacArthur North)/2 − (Emerald + Rainbow)/2, alongside pairwise site comparisons.

- **Physiology:** linear mixed models (`lme4`/`lmerTest`, contrasts with `emmeans`), Cox models for survivorship, and generalised additive mixed models (`mgcv`) for Fv/Fm.
- **Gene expression:** a mixed model for each gene, fitted with `dream` (`variancePartition`) on voom-weighted counts, with Kenward-Roger degrees of freedom. Overall expression differences are tested by PERMANOVA (`vegan`).
- **Sensitivity analyses** in the differential expression scripts compare the primary model with DESeq2, with limma using `duplicateCorrelation`, and with Satterthwaite degrees of freedom. The *O. faveolata* script also checks the effect of pooling the second-cohort resamples and of leaving out each acidified + bleaching sample in turn.
- **Front-loading:** a habitat × treatment `dream` model, plus a permutation test (habitat labels shuffled among colonies, 5,000 times) on colony means under control conditions.
- **Functional and network analyses:** KOG-MWU on signed −log10 p, WGCNA (signed networks, colony-level habitat tests), and OrthoFinder one-to-one orthologs compared between species.

## Reproducing the analyses

Each script is knitted from its own folder and writes to the matching folder in `outputs/`. The physiology scripts can be knitted in any order. The transcriptomics scripts depend on each other and should be knitted in this order:

1. `urban limma ofav.Rmd` and `urban limma ssid.Rmd`
2. `urban kogmwu ofav.Rmd`, `urban kogmwu ssid.Rmd` and `urban frontloading.Rmd`
3. `urban wgcna ofav.Rmd` and `urban wgcna ssid.Rmd`
4. `urban orthofinder.Rmd`

Notes:

- **Run time.** The differential expression scripts are slow, because the mixed model is fitted for about 20,000 genes with Kenward-Roger degrees of freedom. Each has a switch in its setup chunk, `run_sensitivity`. Set to `TRUE` (as committed, matching the knitted reports), the sensitivity analyses run as part of the knit, which takes several hours for *S. siderea* and about ten hours for *O. faveolata*. Set it to `FALSE` for a routine knit of the primary analysis only.
- **Annotation files.** The transcriptomics scripts read gene name and KOG annotations from two separate repositories, which must sit beside this one (same parent folder): `Orbicella-faveolata-annotated-transcriptome` and `Siderastrea-siderea-annotated-transcriptome`.
- **Intermediate files.** Large intermediate `.RData` files are not tracked (see `.gitignore`); they are recreated by knitting the scripts in the order above.
- **Folder names.** The `limma` script and output folder names refer to the limma/voom framework that `dream` builds on.

## Main R packages

`tidyverse`, `lme4`, `lmerTest`, `emmeans`, `MuMIn`, `mgcv`, `gratia`, `survival`, `seacarb`, `vegan`, `DESeq2` (count normalisation and variance-stabilising transformation), `edgeR`, `limma`, `variancePartition`, `arrayQualityMetrics`, `WGCNA`, `KOGMWU`, `pheatmap`, `VennDiagram`, `ggpubr`, `cowplot`.
