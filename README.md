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
| 1i  | `Figure1\lipid_dye\script\make_kymograph.m` | plot the kymograph of ferroptotic trigger wave |

### Figure2
| Figure  | Scripts | Description |
| :------------- | :------------- | :------------- |
| 2a | `Figure2\script\export.ijm` `Figure2\script\plot_profile.m` | the cropped images over time and plot the profile of ROS across gaps over time |
| 2b | `Figure2\script\passing_prob.m` | plot a logistic curve fitted to the gap size that can continue a wave |

### Figure3
| Figure  | Scripts | Description |
| :------------- | :------------- | :------------- |
| 3c and d | `Figure3\script\export_DFO.ijm` `Figure3\script\generate_outline_DFO.m`| export microscopic images and plot the outlines of ferroptotic trigger wave of DFO |
| 3e and f | `Figure3\script\export_FAC.ijm` `Figure3\script\generate_outline_FAC.m` | export microscopic images and plot the outlines of ferroptotic trigger wave of FC |
| 3d, 3f, 3g, 3h and 3i | `Figure3\script\make_kymograph.m` | plot the kymograph of ferroptotic trigger wave under the 5 chemical perturbations (DFO, FC, GKT137831, LY294002, and Dasatinib) |
| 3j to 3n | `Figure3\script\summarize_DFO_96.m` `Figure3\script\summarize_FAC_96.m` `Figure3\script\summarize_GKT_96.m` `Figure3\script\summarize_LY294002_96.m` `Figure3\script\summarize_Dasatinib_96.m`  | the dose-response curves of 5 chemical perturbations |

### Figure4
| Figure  | Scripts | Description |
| :------------- | :------------- | :------------- |
| 4a | `Figure4\priming\script\bifurcation_diagram.m`| Bifurcation diagram |
| 4b | `Figure4\priming\script\export_ros_steadystate.ijm`| export microscopic images (ROS) after photoinduction |
| 4c | `Figure4\priming\script\quantify_ros_staedystate.m`| plot ROS steady state before and after photoinduction |
| 4d | `Figure4\erastin6\script\sim_trigger_waves.m` `Figure4\erastin6\script\export.ijm` `Figure4\erastin6\script\outline.m` | plot simulated trigger waves, export microscopic images and plot the outlines |
| 4e to 4g | `Figure4\erastin6\script\summarize_speed.m` `Figure4\erastin6\script\measure_wavefront_width_amp.m` | plot wave speed, wavefront width, and ROS amplitude under different erastin levels|

### Figure5
| Figure  | Scripts | Description |
| :------------- | :------------- | :------------- |
| 5a and 5b | `Figure5\muscle_4hne\script\export_fig5a.ijm` `Figure5\muscle_4hne\script\export_fig5b.ijm` | export microscopic images |
| 5c and 5d | `Figure5\tunel_4hne\script\export_fig.ijm` | export microscopic images |
| 5e | `Figure5\wave_in_limb\script\processing.ijm` `Figure5\wave_in_limb\script\outline_limb.m` `Figure5\wave_in_limb\script\outline_wave.m` | export microscopic images and plot outlines |
| 5f | `Figure5\death_area\script\summarize_area.m` | plot the ratio of death area to limb area under different conditions |
| 5g to 5i | `Figure5\photoinduction\script\export_v2.ijm` `Figure5\photoinduction\script\plot_ROS_quantification_v2.m` | export microscopic images and plot quantification result |
| 5j | `Figure5\muscle_fiber\script\generate_maxprojection_subset_fiber.ijm` `Figure5\muscle_fiber\script\fiber_orientation_v2.ijm` | export fiber images with colormap representing fiber orientation|
