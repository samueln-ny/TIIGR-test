# computeSystemMatrix_01.py

This script computes **system matrices** from GEANT4 simulation `.root` files (data containing information such as energies and cordinates of interaction of photons in a compton camera based detector). These matrice characterises the detector response when a photon interact with it. This component is essential in the image reconstruction of the source of photons in the detector. The script applies physics-based filtering,  and then a 2D Kernel Density Estimation (KDE) is used to construct a system matrice.

The final output is saved in both **ROOT** and **HDF5** formats.

---

##  Overview

- **Input**: `.root` files from GEANT4 simulations  
- **Parameters**: YAML file defining detector and processing configuration  
- **Output**: System matrices saved as `.root` and `.h5` files  
- **Main techniques**:
  - Physics-based filtering (e.g., position and energy cuts)
  - 2D KDE 


  # FWHM Analysis of 3D Photon Source Reconstructions

This script computes the **Full Width at Half Maximum (FWHM)** for 1D slices of reconstructed 3D photon source distributions. It is designed to evaluate and compare the quality of images produced using different **system matrices** in a simulated **Compton camera-based detector**. The FWHM metric is used as a quantitative indicator of spatial resolution.

## 📘 Description

In simulated photon imaging using Compton cameras, image reconstruction quality depends on the system matrix used. This script analyzes the reconstructions by:

- Extracting 1D slices (X, Y, and Z axes) through the center of the 3D volume.
- Calculating the FWHM using `scipy.signal.find_peaks`.
- Tracking FWHM evolution over multiple training epochs.
- Comparing final image sharpness across different system matrices via plots and metrics.

## 📂 Input

The script expects:

- A directory structure containing subdirectories for each system matrix (e.g., `OUTPUT_0_0_3`, `ORIG_OUTPUT`, etc.).
- Inside each subdirectory: 3D activity reconstructions in **HDF5 format**, named by epoch (e.g., `0001.hdf5`, `0060.hdf5`).
- Each `.hdf5` file must contain:
  - A dataset named `'activities'`
  - An attribute `'voxels'` describing voxel dimensions and grid shape



  



