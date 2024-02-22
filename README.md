# Flux-tower-attribute-for-LSM


The code performs further quality control on the existing flux tower reprocessing dataset (PLUMBER2, Ukkola et al., 2022). It also conducts plant functional types (PFTs) classification of the collected fractional vegetation cover (FVC). The final dataset is then created in NetCDF format by reading csv files. This dataset is intended to elevate the proficiency of flux tower data to serve as benchmarking data for LSMs.

The 'flux_vpd_selected.py' is used to screen the PLUMBER2 (Ukkola et al., 2022) dataset for flux variables (latent and sensible heat) and VPD to obtain site and years that met the requirements.

The 'creat_nc_read' folder holds the attribute data and year information for the selected sites. Run program 'Creat_nc.py', reading the csv file in this folder to create the dataset we have provided, and continue to run program 'Add_BD_SOC.py' to add BD and OC information to the dataset.

The 'PFT_postprocess' folder contains two programs: 'Veg_climate.py' and 'Par_C3C4.F90'. 'Veg_climate.py' categorizes vegetation climate types by reading Köppen climate maps (Beck et al., 2018), based on the method proposed by Poulter et al. (2011). 'Par_C3C4.F90' divides C3 and C4 grasslands by temperature and precipitation measured by flux towers, and site LAI, following the approach described by Still et al. (2003).


**Directory:**
```bash
├── Creat_nc.py                : Read siteinfo.csv, select_year.csv, igbp.csv and lai_qc.csv files to create NetCDF files.
├── Add_BD_SOC.py              : Read bd_soc.csv, adding soil bulk density and organic carbon concentration data to the NetCDF files.
├── creat_nc_read              : The files required to create NetCDF files. 
│   ├── bd_soc.csv             : Soil bulk density and organic carbon concentration, and their reference sources.
│   ├── igbp.csv               : Site IGBP (International Geosphere–Biosphere Programme) classifications. 
│   ├── lai_qc.csv             : QC (quality control) information for site LAI. 
│   ├── select_year.csv        : Sites and years screened.
│   └── siteinfo.csv           : Collected site attribute data and their reference sources.
├── flux_vpd_selected.py       : For flux (latent and sensible heat) and VPD (vapor pressure deficit) screening.
├── PFT_postprocess            : Programs for PFTs classification. 
│   ├── Par_C3C4.F90           : Program for C3/C4 grass segmentation.
│   ├── site_latlon.csv        : Site latitude and longitude
│   └── Veg_climate.py         : Program for vegetation climate classification.
├── README.md                  : Repositories description.
├── Excluded_sites.csv         : Sites and years excluded after screening and the reasons. 
└── Selected_sites.csv         : Compliant sites and years. 
```


<br>

**Usage:**
Run Creat_nc.py and Add_BD_SOC.py in turn to generate the final dataset.
Run flux_vpd_selected.py to get the screened sites and years.
Prepare Köppen-Geiger climate classification maps and run the program to get the vegetation climate type.
Prepare site percent plant functional type cover (PCT_PFT), monthly average precipitation, air temperature, and LAI values. Running the program to divide C3/C4 grass ratios.


**Citation:**

Jiahao Shi, Hua Yuan et al., 2024. A flux tower site attribute dataset intended for land surface modeling. To be submitted.


**References:**

Beck, H. E., Zimmermann, N. E., McVicar, T. R., Vergopolan, N., Berg, A., and Wood, E. F.: Present and future Köppen-Geiger climate classification maps at 1-km resolution, Sci Data, 5, 180214, https://doi.org/10.1038/sdata.2018.214, 2018.

Poulter, B., Ciais, P., Hodson, E., Lischke, H., Maignan, F., Plummer, S., and Zimmermann, N. E.: Plant functional type mapping for earth system models, Geosci. Model Dev., 4, 993–1010, https://doi.org/10.5194/gmd-4-993-2011, 2011.

Still, C. J., Berry, J. A., Collatz, G. J., and DeFries, R. S.: Global distribution of C3 and C4 vegetation: carbon cycle implications, Global Biogeochem. Cycles, 17, 6-1-6–14, https://doi.org/10.1029/2001GB001807, 2003.

Ukkola, A. M., Abramowitz, G., and De Kauwe, M. G.: A flux tower dataset tailored for land model evaluation, Earth Syst. Sci. Data, 14, 449–461, https://doi.org/10.5194/essd-14-449-2022, 2022.
