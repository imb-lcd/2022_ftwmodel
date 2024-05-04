# Overview
The example code of the ferroptotic trigger wave analysis (i.e., wave speed measurement, vector field analysis, and wave simulation)

## Folder description
1. `mathematical_simulation`

   This folder contains a script (`sim_2D_trigger_wave.m`) to simulate a 2D wave given different levels of Erastin and a colormap for plotting simulation results (`custom_parula.mat`).
   The instructions for the simulation and details of model parameters are given in the script.
   
3. `wave_speed_measurement`

   This folder has three subfolders (`script`, `img`, and `data`).
   The `img` folder is to keep the image for wave speed measurement.
   An example image (~2GB) is provided [here](https://figshare.com/ndownloader/files/46021953 "figshare") and kept in the `img` folder for the following analysis.
   The `data` folder contains the information on the initiation of a wave (`setting_s19.mat`) in the example image and the resulting kymograph (`kymograph_s19.mat`).
   The `script` folder contains the scripts for estimating speed from the kymograph (`speed_estimation.m`) and plotting the kymograph (`make_kymograph.m`).
   
5. `vector_field_analysis`

   This folder has two subfolders (`script` and `data`).
   The `data` folder contains the information on the boundary and initiation of a wave (`mask_s11.mat` and `center_s11.mat`).
   The `script` folder contains the scripts for plotting vector fields and quantifying entropy (`vector_field.m` and `entropy.m`) of the death pattern over time.

The provided code runs in MATLAB_R2024a.

OS: Windows 10
CPU: AMD Ryzen 9 5900X
RAM: 32 GB 1800 MHz DDR4

# Setting up the environment
# Demo
## Wave speed measurement
## Vector field and entropy quantification
## Wave simulation
