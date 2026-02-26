# Peptiverse_R
This repo holds the code for plotting all of the images for PeptiVerse in R.
## Project Structure

### Directory Overview

```
Peptiverse_R/
├── R/                          # Shared plotting style
├── scripts/                    # Analysis and processing scripts
├── raw/                        # Raw data files
├── data_processed/             # Cleaned and processed datasets
├── data_distribution/          # Data distribution analysis files
├── figures/                    # Generated plots and visualizations
├── iptm_scatter_csv_all/       # pTM score scatter plot data (CSV format)
└── peptideverse.Rproj          # RStudio project file
```

## Folder Descriptions

### `/R`
Contains reusable R functions and utility code for the project. This includes customized font sizes, colors, and font.

**Usage:** Source functions from this directory in analysis scripts.

### `/scripts`
Main analysis and processing scripts. This directory contains R scripts that:
- Process raw data
- Perform statistical analyses
- Generate figures and reports

**Typical workflow:** Scripts read from `/raw`, process data, save outputs to `/data_processed`, and generate results in `/figures`.

### `/raw`
Original, unmodified raw data files taken from the OPTUNA trails from the original repo [PeptiVerse](https://huggingface.co/ChatterjeeLab/PeptiVerse)


### `/data_processed`
Cleaned, processed, and analyzed datasets ready for use in downstream analyses.

### `/data_distribution`
Original training data for distribution visualization.

### `/iptm_scatter_csv_all`
Processed ipTM scatter plot data in CSV format for correlation analysis.

## Required packages
```
install.packages(c(
  "tidyverse","data.table","yardstick","stringr","readr","scales","dplyr",
  "forcats","jsonlite","patchwork","ggplot2","showtext","sysfonts"
))
```
