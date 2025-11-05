# DLEM-Ag-SIF
Integrating satellite SIF with agroecosystem modeling to constrain carbon-water coupling in Midwest U.S. croplands

---

##  Repository Structure
```
DLEM-Ag-SIF/
│
├── Data/                    # Input data files (.mat, .dlem, boundary shapefiles)
├── Fig/                     # Output figures (maps, statistics, correlations)
│
├── DLEM_AG_read.m           # Function to read DLEM-Ag/DLEM-SIF output files
├── evaluate_regression_metrics.m  # Computes R², RMSE, and NRMSE for model evaluation
├── Flux_Site_proc.m         # Site-level analysis and comparison with AmeriFlux observations
├── Midwest_proc.m           # Regional analysis and map generation for GPP, ET, and WUE
├── tight_subplot.m          # Utility for generating compact multi-panel figure layouts
├── brewermap.m          # ColorBrewer colormaps
│
└── README.md                # This file
```

---

## Requirements
- MATLAB **R2023b** or newer  
- **Toolboxes:**
  - Mapping Toolbox  
  - Statistics and Machine Learning Toolbox  

---

## How to Run
1. Place all input data files into the `/Data/` folder.  
2. Open MATLAB and set this folder as your working directory.  
3. Run the following scripts in order:
   - `Flux_Site_proc.m` — processes AmeriFlux site data and evaluates DLEM-SIF vs DLEM-Ag  
   - `Midwest_proc.m` — generates Midwest-wide maps and statistics (GPP, ET, WUE, SPEI)  
4. All figures will be automatically exported to the `/Fig/` folder.

---

## Output
The scripts generate:
- Annual and multi-year maps of **GPP**, **ET**, and **WUE** for both DLEM-SIF and DLEM-Ag  
- Correlation maps with **SPEI** to assess hydroclimatic sensitivity  
- Model–observation comparison plots and performance metrics  
- Distribution and summary statistics of model outputs  

---
