# Flux-tower-attribute-for-LSM


The 'flux_vpd_selected.py' is used to screen the PLUMBER2 dataset for flux variables (latent and sensible heat) and VPD to obtain site and years that met the requirements.

The 'creat_nc_read' folder holds the attribute data and year information for the selected sites. Run program 'Creat_nc.py', reading the csv file in this folder to create the dataset we have provided, and continue to run program 'Add_BD_SOC.py' to add BD and OC information to the dataset.

The 'PFT_postprocess' folder contains two programs: 'Veg_climate.py' and 'Par_C3C4.F90'. 'Veg_climate.py' categorizes vegetation climate types by reading Köppen climate maps (Beck et al., 2018), based on the method proposed by Poulter et al. (2011). 'Par_C3C4.F90' divides C3 and C4 grasslands by temperature and precipitation measured by flux towers, and site LAI, following the approach described by Still et al. (2003).


Reference:

Beck, H. E., Zimmermann, N. E., McVicar, T. R., Vergopolan, N., Berg, A., and Wood, E. F.: Present and future Köppen-Geiger climate classification maps at 1-km resolution, Sci Data, 5, 180214, https://doi.org/10.1038/sdata.2018.214, 2018.

Poulter, B., Ciais, P., Hodson, E., Lischke, H., Maignan, F., Plummer, S., and Zimmermann, N. E.: Plant functional type mapping for earth system models, Geosci. Model Dev., 4, 993–1010, https://doi.org/10.5194/gmd-4-993-2011, 2011.

Still, C. J., Berry, J. A., Collatz, G. J., and DeFries, R. S.: Global distribution of C 3 and C 4 vegetation: Carbon cycle implications: C 4 PLANTS AND CARBON CYCLE, Global Biogeochem. Cycles, 17, 6-1-6–14, https://doi.org/10.1029/2001GB001807, 2003.


