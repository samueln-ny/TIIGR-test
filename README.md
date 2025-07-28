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
  



