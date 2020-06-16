# larval_phenology

Data and R code for a manuscript in preparation. The manuscript measures how larval _Ixodes scapularis_ questing phenology changes along an elevation gradient. It then presents a mechanistic, temperature-driven model which explains that variation in questing phenology.

 - `data/drag_sampling.csv` gives the raw data from the study. This gives the number of larval ticks found on each day of sampling.
- `data/phenology_fits.RData` gives the fit phenology curves to `data/drag_sampling.csv`. 
- `data/processed_prism.RData` gives the mean daily temperature at the three elevation classes. 
- `r_code.R` gives all the code for the analysis, phenology model, and to make the figures.
- `flow_diagraim/flow_diagram.tex` is TeX code to make the flow diagram in Figure 1.
