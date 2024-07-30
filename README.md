# andes_paleoelevation

 Repository containing data and code to reproduce the analysis presented in Liu et al (published online 29 July 2024). The code generates maps of paleoelevation from igneous geochemistry, with specific application to the western margin of South America 0-350 Ma

The published paper can be found here:
https://www.sciencedirect.com/science/article/pii/S1342937X24002016

 ## Contents

- Jupyter notebooks are located in the 'Figures' folder
- Much of the underlying python code is in the 'python' folder
- The 'luffi' folder contains data files required to derive paleoelevations from the Luffi and Ducea (2022) calibration
- The boschman folder contains data required to compare results with independent paleoelevations for 0-80 Ma
- datafiles contains tables of data captured from individual studies and compilations


## Requirements
- The best place to start is the 'Minimal_Example' notebook. This should work without requiring too many exotic modules.
- More complex examples require specific functions for plotting, raster manipulation reconstructing data. See 'module not found' messages for more information. 


## Method Notes

The code implements calibrations according to three published works
- Farner and Lee (2017), EPSL
- Hu et al (2020), GRL
- Luffi and Ducea (2022), Rev. Geophys.

Data filtering implemented within the code is mostly complete, however it is still safer to filter the data (according to stated criteria in above studies) first before loading


