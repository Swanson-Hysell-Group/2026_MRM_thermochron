# Montreal River Monocline Thermochronology

Data analysis and modeling code for study on the thermochronology of the Montreal River monocline in northern Wisconsin associated with the manuscript:

Alemu, T. B., Guenthner, W. R., & Swanson-Hysell, N. L. (2025). Mesoproterozoic, Neoproterozoic, and Mesozoic erosion of the Montreal River monocline in the Laurentia midcontinent. *in prep for Geochemistry, Geophysics, Geosystems*.

## Notebooks

The repository contains the notebooks:

- [📓 MRM_UTh_He_data.ipynb](MRM_UTh_He_data.ipynb) This notebook imports the reduced apatite and zircon (U-Th)/He data developed for the project directly from the https://github.com/wrguenthner/HAL_data repository. Within the https://github.com/wrguenthner/HAL_data/tree/main/projects_2023/TA23_MRM_MCR subdirectory of that repository are the raw datafiles and associated processing notebooks. Plots of these data are developed in the notebook including the main data figure of the paper: 

![Age vs. eU Plots](./figures/markdown_images/MRM_UTh_He_data_image.png)

- [📓 MRM_forward_models.ipynb](MRM_forward_models.ipynb) This notebook implements forward model scenarios using the `pythermo` package that can be compared to the apatite and zircon (U-Th)/He data. The code in this notebook generates the forward model figures of the manuscript.