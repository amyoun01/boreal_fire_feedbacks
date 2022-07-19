# boreal_fire_feedbacks

This repository houses a collection of scripts that conducts an analysis aimed at characterizing how negative fire-vegetation feedbacks may moderate future climate-driven fire regime changes in boreal forest ecosystems of North America.

## Organization

The code for this analysis is organized into three different directories:

```bash
sh/
py/
R/
```

While this analysis is primarily conducted using Python 3.8, there are some additional statistical analyses run in R. Furthermore, we use shell scripts as a wrapper to independently organize and run different segments of the analysis, example:

 ```bash
sh/1_organize_fire_data.sh # Convert fire perimiter shapefiles into raster datasets
sh/2_organize_era5_data.sh # Read, extract, and organize ERA5 climate data
sh/3_organize_cmip6_data.sh # Organize CMIP5 GCM output
# etc ...
```