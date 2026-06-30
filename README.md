# Scripts4Paper_Hansetal_OSSE-AMOC11S
Scripts used in: "Evaluating Transport Observations of the Atlantic Meridional Overturning Circulation at 11°S Using an Ocean Model"

These are the scripts used for the analyses and figures in Hans et al. (2026) "Evaluating Transport Observations of the Atlantic Meridional Overturning Circulation at 11°S Using an Ocean Model". The study mimics an observational array at 11°S in a high-resolution ocean model to assess the capability of the array to compute the temporal variability of the Atlantic Meridional Overturning Circulation. Accordingly, an observing system simulation experiment is carried out. Various observational approaches are tested based on subsampled bottom pressure recorders, subsampled moored temperature and salinity sensors, and subsampled inverted echo sounders.


## Overview of the provided scripts:

1. To compute bottom pressure and geostrophic transport from the model output:

  - gen_geostrophy.sh  
  - gen_bp.sh

2. To apply the different observational strategies to the model data:

  - cut_ModelTruth.m  
  - f_cut_BPRSetup.m  
  - run_f_cut_BPRSetup.m  
    requires:  
    f_cut_MooringSetup(2,2,0,1,0,1,1,'')  
    f_cut_MooringSetup(2,3,0,1,0,1,1,'')  
    f_cut_MooringSetup(2,0,1,1,0,1,1,'')  
  - f_cut_MooringSetup.m  
  - run_f_cut_MooringSetup.m  
    requires:  
    f_cut_BPRSetup(2,2,2,1,0,1,1,4,'')  
  - cut_IESTest.m

3. To compute statistics for the tables of the manuscript:

  - f_stats_calc.m  
  - f_stats_BPRSetup.m  
  - run_f_stats_BPRSetup.m  
  - f_stats_MooringSetup.m  
  - run_f_stats_MooringSetup.m

4. To generate figures for the manuscript:

  - plot_Figure_1.m (first run cut_TRACOS_Vsection.m)  
  - plot_Figure_2.m  
  - plot_Figure_3.m  
  - plot_Figure_4.m  
  - plot_Figure_5.m  
  - plot_Figure_6.m  
  - plot_Figure_7.m (first run cut_TRACOS_Vsection.m)  
  - plot_Figures_Appendix.m

## Additional remarks:
For 1), CDFTOOLS is required. For 2) the seawater library is required. For 4), the colour maps cmocean and Scientific colour maps are required. All references are given in the manuscript. Except for the two bash files, the code is written in MATLAB_R2021b.

The scripts require the following data sets ([link to data](https://hdl.handle.net/20.500.12085/f72c4932-6132-4f98-a7c1-2cdd9c1c377f)):  
  - 1_VIKING20X.L46-KFS003_1d_${YYYY}0101_${YYYY}1231_detrendedSSH_11S_box.nc  
  - 1_VIKING20X.L46-KFS003_1d_${YYYY}0101_${YYYY}1231_grid_T_11S_box.nc  
  - 1_VIKING20X.L46-KFS003_1d_${YYYY}0101_${YYYY}1231_grid_U_11S_box.nc  
  - 1_VIKING20X.L46-KFS003_1d_${YYYY}0101_${YYYY}1231_grid_V_11S_box.nc  
  - 1_VIKING20X.L46-KFS003_mesh_mask_11S_box.nc
