# larval_phenology

Data and R code for a manuscript in submission. The manuscript measures how larval _Ixodes scapularis_ questing phenology changes along an elevation gradient. It then presents a mechanistic, temperature-driven model which explains that variation in questing phenology.

* `data/drag_sampling.csv` gives the raw data from the study. Rows are individual 200 m<sup>2</sup> tick drag cloth samples. The `csv` has the following columns.
  * `site` the name of the location where the drag-cloth sample took place.
  * `elev` the elevation, in meters, of the site.
  * `data` the date when the sample took place. In year-month-day format.
  * `julian` the julian day of the sample.
  * `larva` the number of larval ticks found in the 200 m<sup>2</sup> sample.
* `data/leaf_litter_temp.csv` this gives the estimated daily mean leaf-litter temperature at each site. These are generated from PRISM modeled temperature and an air to leaf-litter temperature difference. See manuscript for details. The columns are:
  * `jday` julian day.
  * `tmean` estimated mean leaf litter temperature (C). 
  * `site` name of site.
  * `elevCat` elevation category of site.
* `r_code.R` gives all the code for the analysis, phenology model, and to make the figures.
* `flow_diagraim/flow_diagram.tex` is TeX code to make the flow diagram in Figure A1.
