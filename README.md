This script computes **system matrices** (KDE-based system matrix and Histogram based system matrix) from GEANT4 simulation `.root` files (data containing information such as energies and coordinates of interaction of photons from our simulation(the simulation setup was designed to model the interaction of gamma rays within a two layer liquid xenon (LXe) detector composed of two concentric cylinders). These system matrix encodes the probability that a photon emitted from a particular point in the detector is detected with a specific scattering angle and energy. First of all the data obtained from the simulation is filtered using physics-based filtering to make it more realistic like it was from a real detector and then we construct the system matrices from the preprocessed data using the energies deposited and the angle at which photons were scattered.

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



  



