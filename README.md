# boreal_fire_feedbacks

### Contact: Adam Young, ayoung49[at]alaska.edu

---

This repository encompasses a collection of scripts to conduct an analysis aimed at characterizing how negative fire-vegetation feedbacks may moderate future climate-driven fire regime changes in boreal forest ecosystems of North America.

## Organization

The code for this analysis is organized into three different directories:

```bash
sh/
py/
R/
```

While this analysis is primarily conducted using Python 3.8, there are some additional statistical analyses run in R. Furthermore, we use shell scripts as a wrapper to independently organize and run these python and R scripts for  different segments of the analysis, e.g.,:

 ```bash
$ sh/1_organize_fire_data.sh # Convert fire perimiter shapefiles into raster datasets
$ sh/2_organize_era5_data.sh # Read, extract, and organize ERA5 climate data
$ sh/3_organize_cmip6_data.sh # Organize and bias-correct CMIP6 GCM output
# etc ...
```

This README file and parts of the analysis will be regularly updated and added as this project
continues.

