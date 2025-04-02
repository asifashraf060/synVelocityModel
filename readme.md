# 3D Velocity Model Construction for StingRay & TomoLab

## Project Overview
This repository contains MATLAB scripts (tested on **MATLAB 2022a**) for constructing a 3D initial seismic velocity model, intended for use in local earthquake tomography. The scripts integrate various geophysical datasets (e.g. shallow velocity profiles, crustal interfaces, station/event info) to build a continuous P-wave velocity model. The resulting model and associated files are formatted for **StingRay** (a 3D ray-tracing toolkit) and **TomoLab** (a MATLAB-based travel-time tomography inversion library). Researchers in geophysics can use this toolkit to prepare a starting model and input files before running inversion with StingRay/TomoLab.

## Script Descriptions
Below is a brief description of each script in the repository and its purpose:

- **`make_upCrustModel.m`** – Generates a smooth 2D upper-crust velocity model from sparse data points. This script reads an input table (e.g. Excel file) of P-wave velocities at given depths and longitude positions, then performs interpolation and smoothing (Gaussian filtering) to create a refined velocity profile for the upper crust. It saves the result as a MATLAB structure (e.g. `upCrust_*.mat`) for later use.

- **`make_interface.m`** – Processes geological interface data (e.g. Moho depth, slab geometry, topography) to create a merged 3D interface structure. It reads input files such as a slab model `.mat`, longitude/latitude grids, elevation data, and mapping geometry. The script outputs a consolidated interface model (`int_3D_*.mat`) defining key discontinuities (e.g. basement and Moho surfaces) across the region. Diagnostic plots of the interfaces are also produced for verification.

- **`make_srModel.m`** – Assembles the full 3D velocity model (the “srModel”) by integrating the upper crust model and interface structure with additional data. This script loads the upper crust velocity matrix, the 3D interfaces, and supporting datasets (station locations, event locations, elevation grids, and StingRay control parameters). It fills in velocity values across the entire model grid, assigning velocities in sediments, crust, and mantle, and interpolating across gaps. Optional steps include smoothing sharp boundaries, patching or adjusting velocities, applying a Gaussian filter, or superimposing a checkerboard perturbation for resolution testing. The final **srModel** (a structured MATLAB variable containing grids of slowness/velocity and other info) is saved to `outputs/entireModel/models/` (e.g. `srModel_<date>_<region>.mat`) for use in StingRay.

- **`apply_srModel_checkerboard.m`** – Utility function to apply a 3D checkerboard velocity perturbation to an existing `srModel`. Given a base model, it perturbs the P-wave velocities by a specified percentage in a checkerboard pattern (controlled by wavelength parameters). This can be applied to the whole model or restricted to certain depth ranges (e.g. only the upper mantle). The output is a modified velocity matrix (or updated srModel) useful for synthetic resolution tests.

- **`apply_srModel_mAvg.m`** – Utility function to apply a moving-average smoothing filter to an `srModel`. It convolves the velocity grid with a defined window size in a specified direction (vertical or horizontal), optionally followed by Gaussian filtering. This helps to create a more gradational model (e.g. “softening” layer boundaries) and can be used to prepare a smoothed initial model variant.

- **`make_tlPert.m`** – Prepares the **TomoLab perturbation grid** structure (`tlPert`). This script defines the inversion grid for tomography by down-sampling or coarsening the `srModel` domain. It creates a set of nodes in X, Y, and Z (depth) at specified intervals (e.g. 2 km spacing) covering the model area, and assigns initial uncertainty values to each node. It also extracts the main interface (e.g. Moho) from the model to include as an interface parameter in the inversion, with its own node spacing and uncertainty. The result is saved as `tlPert_*.mat`, which TomoLab uses to parameterize velocity perturbations during inversion.

- **`make_tlArrival.m`** – Converts observed travel-time picks into the **TomoLab arrival** structure (`tlArrival`). Using the station and event information and the prepared model, this script calls StingRay/TomoLab utility functions (e.g. `tlPick2tlArrival`) to read travel-time pick files (from a specified directory) and bundle them into a unified structure. It supports multiple seismic phases (e.g. direct P, Pg, Pn, etc.) as specified by the user. The output `tlArrival_*.mat` contains the travel times and associated metadata, ready for input to the tomography inversion in TomoLab.

- **`prep_tomoLab_files.m`** – A driver script that ties everything together. It sets up file paths and parameters, then runs the above scripts in the correct sequence to generate the final model and tomography inputs. Specifically, it will call `make_srModel.m` to build and save the 3D velocity model, then `make_tlPert.m` to produce the perturbation grid, and `make_tlArrival.m` to produce the arrival-time data structure. This script is convenient for automating the entire preparation workflow once the input data and configuration are set.

## Dependencies
- **MATLAB 2022a** – The code is developed and tested with MATLAB R2022a and requires this (or a compatible) version.
- **StingRay Toolbox** – The **StingRay** ray-tracing and TomoLab inversion utilities must be installed and on the MATLAB path. These scripts rely on StingRay’s functions such as `load_srModel`, `load_srStation`, `load_srEvent`, `tlPick2tlArrival`, etc. You can obtain StingRay/TomoLab from the University of Oregon (see the [StingRay project page](https://pages.uoregon.edu/drt/Stingray/)). Ensure the StingRay toolbox is properly configured in MATLAB.
- **Input Data Files** – Actual geophysical input data are not included in this repository. The user is expected to have:
  - An upper-crust velocity profile dataset (e.g. an Excel or `.mat` table of depth vs velocity at various locations).
  - A crustal interface model (e.g. slab/Moho depth grid in a `.mat` file, plus supporting lat/lon grids).
  - Station coordinates and earthquake event lists (as StingRay `.mat` structures).
  - Elevation/topography data for the region (`.mat` file for StingRay).
  - A StingRay control parameter file (e.g. defining model grid extents, etc.).
  - Travel-time pick files for the seismic arrivals to be inverted (for use in `tlPick2tlArrival`).
  
  **Note:** Paths to these input files are set inside the scripts (in the “Input” sections of `make_interface.m`, `make_upCrustModel.m`, etc.). Before running, update those file names and ensure the files are placed in the expected `inputs/` subdirectories.

## How to Use
**Execution Order:** The typical workflow is to run the scripts in the following sequence to build the model and prepare inversion inputs:
1. **Prepare Input Data:** Ensure all required input files are available and update the script parameters (file names, region limits, etc.) as needed. For example, set the correct Excel file name in `make_upCrustModel.m` and the correct interface file and region bounds in `make_interface.m`.
2. **Generate Upper Crust Model:** Run `make_upCrustModel.m` first. This will create a refined upper crust velocity profile and save a file (in `outputs/upCrust/structures/`) containing the 2D velocity matrix for the upper crust.
3. **Generate Interface Structure:** Next, run `make_interface.m`. This processes the slab/Moho and elevation data to produce a merged interface structure file (saved in `outputs/slab/structures/`). Check the output plots to ensure the interface surfaces cover your region of interest.
4. **Assemble 3D Model:** Update `make_srModel.m` (or `prep_tomoLab_files.m`) to point to the newly created upper crust and interface files (the script uses `which()` to find files by name). Then run **either** `make_srModel.m` standalone or use the `prep_tomoLab_files.m` driver:
   - Running `make_srModel.m` will construct the full 3D velocity model (`srModel`). You may interactively inspect the model or plots it generates for quality control (the script checks for NaNs and shows model coverage).
   - If you use `prep_tomoLab_files.m`, it will internally call `make_srModel.m` followed by the next steps automatically.
5. **(Optional) Modify Model:** If you wish to test alternative initial models, you can use the utility functions:
   - Apply a checkerboard perturbation by calling `apply_srModel_checkerboard` on the saved `srModel` structure (or by enabling the checkerboard option within `make_srModel.m` if configured). This produces a variant of the model useful for resolution tests.
   - Smooth the model by calling `apply_srModel_mAvg` with a chosen window size and direction to remove sharp velocity contrasts.
   These steps are optional and mainly for experimentation or sensitivity testing.
6. **Prepare Tomography Inputs:** After obtaining the final velocity model, run `make_tlPert.m` to generate the inversion grid (`tlPert`) and then `make_tlArrival.m` to convert pick data to the arrival structure. (If using `prep_tomoLab_files.m`, these are done in order automatically.) Ensure that `tlPickDir` (the directory of pick files) and `phaseIn` (the list of seismic phase codes) are correctly set in the script before running. The scripts will save `tlPert` and `tlArrival` files in the `outputs/tomolab/` directory.
7. **Ready for Inversion:** With `srModel`, `tlPert`, and `tlArrival` prepared, you can proceed to run the TomoLab inversion routines (outside this repository’s scope) to invert the travel times for velocity updates. Refer to StingRay/TomoLab documentation for running the inversion using these input structures.

**Notes:**  
- Each script is designed to be run independently if needed (ensuring the input paths are set). However, the `prep_tomoLab_files.m` provides a convenient way to run the critical model-building steps sequentially without manual intervention.  
- Be mindful of coordinate systems: StingRay works in a local Cartesian coordinate system. The scripts use a provided `srGeometry` (mapping from lat/lon to X/Y) to ensure consistency. Always use the same geometry reference for all inputs (station locations, interfaces, etc.).
- Before inversion, double-check the output model and picks for consistency (e.g., station and event coordinates should lie within the model bounds, no NaN velocities, etc.).

## Output Files
All results are saved under the `outputs/` directory in structured subfolders. Key output files include:

- **3D Velocity Model (`srModel_*.mat`):** The primary output is a MATLAB `.mat` file containing the `srModel` structure. This includes the 3D grid coordinates (local X, Y, Z), P-wave slowness (or velocity) values, topography, and interface definitions. This file serves as the initial velocity model for StingRay and as the starting model for the tomography inversion. *(Location: `outputs/entireModel/models/`)*

- **TomoLab Perturbation Grid (`tlPert_*.mat`):** A MATLAB structure defining the inversion grid nodes and their uncertainties. It contains the subsampled X, Y, Z grid used by TomoLab to solve for velocity perturbations, along with uncertainty values (or damping factors) for each node. It also includes a representation of the main velocity interface (e.g. Moho) with its own grid and uncertainties. *(Location: `outputs/tomolab/`)*

- **TomoLab Arrival Data (`tlArrival_*.mat`):** A structure holding all picked travel times and associated information (phase, station, event, etc.) formatted for TomoLab. This is derived from the user-provided pick files and is required by the inversion algorithm to compute misfits and update the model. *(Location: `outputs/tomolab/`)*

- **Upper Crust Model (`upCrust_*.mat`):** The intermediate upper-crust velocity model file produced by `make_upCrustModel.m`. This 2D (depth vs. distance) velocity profile is used as input for building the 3D model. *(Location: `outputs/upCrust/structures/`)*

- **Interface Structure (`int_3D_*.mat`):** The merged interface surfaces file from `make_interface.m`. It contains grids (in X–Y) for the depth of key interfaces (e.g. Moho, slab top, etc.) across the region. This is incorporated into the 3D model assembly. *(Location: `outputs/slab/structures/`)*

- **Diagnostic Plots:** Several JPEG/PNG plots are generated to visualize the model and interfaces (e.g. `upCrust_model.jpg`, interface maps, cross-sections). These help in verifying that the interpolation and model assembly are correct. *(Location: under corresponding `outputs/*/plots/` directories)*

All output files are meant for **research use**. They provide the necessary starting model and parameters for running a StingRay/TomoLab inversion. Once these files are generated, users can proceed with StingRay (for ray tracing or synthetic travel times) and TomoLab (for iterative tomographic inversion) to refine the velocity model using observed data. This toolkit streamlines the preparation stage, ensuring that the inversion has a geologically reasonable starting model and properly formatted input data.
