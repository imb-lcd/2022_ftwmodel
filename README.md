# Overview
The example code of the ferroptotic trigger wave analysis (i.e., wave speed measurement, vector field analysis, and wave simulation)
This is a collection of customized scripts used for image analysis, simulation, and plotting.

The provided code runs in MATLAB_R2024a.  
OS: Windows 10  
CPU: AMD Ryzen 9 5900X  
RAM: 32 GB 1800 MHz DDR4

## Folder description
| Figure | Folder | Description |
| :------------- | :------------- | :------------- |
| Figure 1a to 1e | `Figure1\long_distance\` | Generate microscopic images and plot outlines, kymograph, intensity profile, and comparison between experimental data and theoretical diffusion. |
| Figure 1f to 1i | `Figure1\lipid_dye\` | Generate microscopic images and plot outlines, kymograph, and intensity profile.  |
| Figure 2 | `Figure2\` | Generate microscopic images, plot intensity profile, and fit a logistic curve. |
| Figure 3 | `Figure3\` | Generate microscopic images, plot outlines and kymograph, and fit a dose-response curve. |

Please see `readme.txt` inside the folder for detailed instructions to reproduce the figures.

| Folder  | Scripts | Description |
| :------------- | :------------- | :------------- |
| Figure1 | 1.`\long_distance\script\export.ijm` 2.`\long_distance\script\outline.m` 3.`\long_distance\script\get_slices.m` 4.`\long_distance\script\make_kymograph.m` 5.`\long_distance\script\diffusion_plot.m` 6.`\lipid_dye\script\outline_lipid.m` 7.`\lipid_dye\script\plot_profile.m` 8.`\lipid_dye\script\outline_cy5.m` 9.`\lipid_dye\script\get_slices.m` 10.`\lipid_dye\script\make_kymograph.m`| 1.export microscopic images; 2.plot the outlines of ferroptotic trigger wave; 3.plot ferroptotic trigger wave with an ROI; 4.plot the kymograph of ferroptotic trigger wave; 5.plot the data of ferroptotic trigger wave with theoretical diffusion; 6.plot the outlines of ferroptotic trigger wave (lipid dye); 7.plot profile of ferroptotic trigger wave (cell death and lipid); 8.plot the outlines of ferroptotic trigger wave (cell death); 9.plot ferroptotic trigger wave with an ROI (cell death and lipid); 10.plot the kymograph of ferroptotic trigger wave| 

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
| 1i  | `Figure1\lipid_dye\script\make_kymograph.m` | plot the kymograph of ferroptotic trigger wave |

### Figure2
| Figure  | Scripts | Description |
| :------------- | :------------- | :------------- |
| 2a | `Figure2\script\export.ijm` `Figure2\script\plot_profile.m` | the cropped images over time and plot the profile of ROS across gaps over time |
| 2b | `Figure2\script\passing_prob.m` | plot a logistic curve fitted to the gap size that can continue a wave |

### Figure3
| Figure  | Scripts | Description |
| :------------- | :------------- | :------------- |
| 3c and d | `Figure3\script\export_DFO.ijm` `Figure3\script\outline_DFO.m` `Figure3\script\make_kymograph_DFO.m`| export microscopic images and plot the outlines and kymograph of ferroptotic trigger wave of DFO |
| 3e and f | `Figure3\script\export_FAC.ijm` `Figure3\script\outline_FAC.m` `Figure3\script\make_kymograph_FAC.m` | export microscopic images and plot the outlines and kymograph of ferroptotic trigger wave of FC |
| 3g, 3h and 3i | `Figure3\script\make_kymograph.m` | plot the kymograph of ferroptotic trigger wave under the chemical perturbations (GKY, LY294002, and Dasatinib) |
| 3j to 3n | `Figure3\script\summarize_DFO_96.m` `Figure3\script\summarize_FAC_96.m` `Figure3\script\summarize_GKT_96.m` `Figure3\script\summarize_LY294002_96.m` `Figure3\script\summarize_Dasatinib_96.m`  | the dose-response curves of 5 chemical perturbations |
