# Overview
The example code of the ferroptotic trigger wave analysis (i.e., wave speed measurement, vector field analysis, and wave simulation)
This is a collection of customized scripts used for image analysis, simulation, and plotting.

The provided code runs in MATLAB_R2024a.  
OS: Windows 10  
CPU: AMD Ryzen 9 5900X  
RAM: 32 GB 1800 MHz DDR4

## Folder description
### Figure1
| Figure  | Scripts | Description |
| :------------- | :------------- | :------------- |
| 1a and 1b  | `Figure1\long_distance\script\export.ijm` `Figure1\long_distance\script\outline.m`  | plot the outlines of ferroptotic trigger wave and export microscopic images |
| 1c  | `Figure1\long_distance\script\get_slices.m`  | plot ferroptotic trigger wave with an ROI |
| 1d  | `Figure1\long_distance\script\make_kymograph.m`  | plot the kymograph of ferroptotic trigger wave |
| 1e  | `Figure1\long_distance\script\diffusion_plot.m`  | plot the data of ferroptotic trigger wave with theoretical diffusion |
| 1f  | `Figure1\lipid_dye\script\outline_lipid.m` `Figure1\lipid_dye\script\plot_profile.m`  | plot the outlines and profile of ferroptotic trigger wave (cell death and lipid) |
| 1g  | `Figure1\lipid_dye\script\outline_cy5.m` | plot the outlines of ferroptotic trigger wave (cell death) |
| 1h  | `Figure1\lipid_dye\script\get_slices.m` | plot ferroptotic trigger wave with an ROI (cell death and lipid) |
| 1i  | `Figure1\lipid_dye\script\make_kymograph.m`  | plot the kymograph of ferroptotic trigger wave |

### Figure2
| Figure  | Scripts | Description |
| :------------- | :------------- | :------------- |
| 2a | `Figure2\script\plot_profile.m` | plot the profile of ROS across gaps over time |
| 2b | `Figure2\script\passing_prob.m` | plot a logistic curve fitted to the gap size that can continue a wave |
