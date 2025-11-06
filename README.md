
# Plastic stomach contraction and hypertrophy analysis

This repository contains analysis code and outputs for Sable Shearwater
(*Ardenna carneipes*) research, known in te reo Māori as **Toanui**. The
analysis measures stomach contraction following AChE stimulation and
relates stomach weight to plastic ingestion.

This work supports the manuscript:

**From pollution to pathology: biomechanical consequences of plastic
ingestion in a seabird**  
Ingrid L. Pollet¹², William M. Connelly³, Michelle L. Feenstra, Alix M.
de Jersey²³, Alexander L. Bond²⁴, Jennifer L. Lavers²⁴⁵, Jack Rivers
Auty²³

1.  Acadia University, Biology Department, Wolfville, Nova Scotia, B4P
    2R6, Canada  
2.  Adrift Lab, Underwood, Lutruwita/Tasmania, 7268, Australia  
3.  Tasmanian School of Medicine, University of Tasmania, Hobart,
    Tasmania, 7000, Australia  
4.  Bird Group, The Natural History Museum, Tring, Hertfordshire, HP23
    6AP, United Kingdom  
5.  Gulbali Institute, Charles Sturt University, Wagga Wagga, New South
    Wales, 2678, Australia

## What this repo does

The R Markdown document
**`stomach_contraction_analysis_FIGS_v3_1.Rmd`**:

- Loads raw time-series, video mapping, and mass data files.  
- Normalises bird IDs and collapses duplicate rows per bird by summing
  numeric fields.  
- Baseline-corrects traces per video at *t* = 0 (or first non-NA).  
- Groups birds by plastic mass threshold (≤ 0.5 g vs \> 0.5 g).  
- Aggregates contraction and relaxation phases to mean ± SEM.  
- Produces publication figures and a hypertrophy CSV.

## Repository structure

    .
    ├─ stomach_contraction_analysis_FIGS_v3_1.Rmd   # Main analysis
    ├─ final_trim.csv                               # Time-series data (raw and smoothed columns)
    ├─ video_map.csv                                # VIDEO to bird mapping with PHASE
    ├─ masses.csv                                   # Per-bird weights and counts
    ├─ Figure 1.pdf                                 # Scatter + contract mean ± SEM (legend below)
    ├─ Figure 2.pdf                                 # ABCD square with legend in panel D (270 × 180 mm)
    └─ stomach_hypertrophy.csv                      # Export: ID, plastic mass/number, stomach weight

## How to run

1.  Open the `.Rproj` in RStudio (or run `usethis::proj_activate()` in
    the folder).  

2.  Install required packages once per machine:

    ``` r
    install.packages(c("ggplot2","cowplot","viridisLite","sjPlot","grid","usethis"))
    ```

3.  Knit `stomach_contraction_analysis_FIGS_v3_1.Rmd` to HTML. This also
    writes the PDFs and CSV listed above.

### Data file expectations

- `final_trim.csv`: wide table with columns ending in `_raw` and
  `_smoothed`. Time is implicit by row index.  
- `video_map.csv`: has `VIDEO`, `BIRD_ID`, `PHASE` columns. IDs are
  normalised to a `BIRD_KEY`.  
- `masses.csv`: includes `ID`, `Total_mass`, `Total_number`, `Wt`
  (stomach weight), and related fields. Duplicate rows per bird are
  summed.

## Figure notes

- **Figure 1**: Panel A shows stomach weight vs plastic mass with a
  linear fit. Panel B shows mean contraction with SEM by plastic group
  on a minutes scale. Legend sits below panel B.  
- **Figure 2**: Square layout A–D with a **standalone centred legend**
  in panel D. Saved at **270 × 180 mm**.

All figures are saved with embedded fonts suitable for journal
submission when using a PDF device in `ggsave()`.

## Reproducibility

- Baseline correction: per `VIDEO`, baseline is taken at time 0 if
  non-NA, else the first non-NA point. Traces are then expressed as Δ
  from baseline.  
- Aggregation: mean ± SEM per time and plastic group.  
- Time base: one unit equals 30 seconds, so minutes = `time / 2`.  
- Plastic grouping: `Total_mass` ≤ 0.5 g vs \> 0.5 g.

## How to cite

If you use this code or outputs, please cite the manuscript above once
published. You can also cite the software directly from within R:

``` r
citation("ggplot2"); citation("cowplot"); citation("sjPlot")
```

## Contact

For questions about analysis or data, contact **Jack Rivers Auty**.  
For manuscript queries, contact **Ingrid L. Pollet**.
