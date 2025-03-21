# Data and code for "A mechanistic model explains variation in larval tick questing phenology an elevation gradient"

Data and R code for a manuscript in press. The manuscript measures how larval _Ixodes scapularis_ questing phenology changes along an elevation gradient. It then presents a mechanistic, temperature-driven model which explains that variation in questing phenology.

* `data/drag_sampling.csv` gives the raw data from the study. Rows are individual 200 m<sup>2</sup> tick drag cloth samples. The `csv` has the following columns.
  * `site` the name of the location where the drag-cloth sample took place.
  * `elev` the elevation, in meters, of the site.
  * `data` the date when the sample took place. In year-month-day format.
  * `julian` the julian day of the sample.
  * `larva` the number of larval ticks found in the 200 m<sup>2</sup> sample.
* `data/leaf_litter_temp.csv` this gives the estimated daily mean below-canopy, leaf-litter temperature at each site. These are generated from PRISM modeled temperature and an adjustment based on the difference between above- and below-canopy temperature. See manuscript for details. The columns are:
  * `jday` julian day.
  * `tmean` estimated mean below-canopy, leaf-litter temperature (C). 
  * `site` name of site.
* `r_code.R` gives all the code for the analysis, phenology model, and to make the figures.
