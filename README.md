# Flux-tower-attribute-for-LSM

[说明整个软件包实现的功能]

The 'flux_vpd_selected.py' is used to screen the PLUMBER2 dataset for flux variables (latent and sensible heat) and VPD to obtain site and years that met the requirements.

The 'creat_nc_read' folder holds the attribute data and year information for the selected sites. Run program 'Creat_nc.py', reading the csv file in this folder to create the dataset we have provided, and continue to run program 'Add_BD_SOC.py' to add BD and OC information to the dataset.

The 'PFT_postprocess' folder contains two programs: 'Veg_climate.py' and 'Par_C3C4.F90'. 'Veg_climate.py' categorizes vegetation climate types by reading Köppen climate maps (Beck et al., 2018), based on the method proposed by Poulter et al. (2011). 'Par_C3C4.F90' divides C3 and C4 grasslands by temperature and precipitation measured by flux towers, and site LAI, following the approach described by Still et al. (2003).

├── Add_BD_SOC.py              : 一句话说明
├── Creat_nc.py                : 下同，如果某些文件确定不需要可删除不说明
├── creat_nc_read
│   ├── bd_soc.csv
│   ├── igbp.csv
│   ├── lai_qc.csv
│   ├── select_year.csv
│   └── siteinfo.csv
├── Excluded_sites.csv
├── Excluded_sites.xlsx
├── flux_vpd_selected.py
├── PFT_postprocess
│   ├── Par_C3C4.F90
│   └── Veg_climate.py
├── README.md
├── Selected_sites.csv
├── Selected_sites.xlsx
└── select_year.csv


[添加用法说明]
Usage:


Citation:

Jiahao Shi, Hua Yuan et al., 2024. A flux tower site attribute dataset intended for land surface modeling. To be submitted.


References:

Beck, H. E., Zimmermann, N. E., McVicar, T. R., Vergopolan, N., Berg, A., and Wood, E. F.: Present and future Köppen-Geiger climate classification maps at 1-km resolution, Sci Data, 5, 180214, https://doi.org/10.1038/sdata.2018.214, 2018.

Poulter, B., Ciais, P., Hodson, E., Lischke, H., Maignan, F., Plummer, S., and Zimmermann, N. E.: Plant functional type mapping for earth system models, Geosci. Model Dev., 4, 993–1010, https://doi.org/10.5194/gmd-4-993-2011, 2011.

Still, C. J., Berry, J. A., Collatz, G. J., and DeFries, R. S.: Global distribution of C 3 and C 4 vegetation: Carbon cycle implications: C 4 PLANTS AND CARBON CYCLE, Global Biogeochem. Cycles, 17, 6-1-6–14, https://doi.org/10.1029/2001GB001807, 2003.


