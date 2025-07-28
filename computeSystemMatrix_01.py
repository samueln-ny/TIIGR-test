### Bandwidth: 0.1
import ROOT as R
import numpy as np
import uproot as ur
import argparse
import h5py
import os
from sklearn.neighbors import KernelDensity
import matplotlib.pyplot as plt
import time
from tiigr import ProcessingParams
from tiigr import EventProcessor
from tiigr import RunInfos
import sys
import pandas as pd

dirPathOfThisFile = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    script_start_time = time.time()
    R.EnableImplicitMT()

    nBinsEnergy = 10000
    nBinsAngle = 1800

    parser = argparse.ArgumentParser(
        prog="computeSystemMatrix_01.py",
        description="Compute a system matrix from one or multiple root files from GEANT4",
    )

    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument("-i", "--input", nargs="+", required=True)
    required.add_argument("-o", "--output", required=True)
    required.add_argument(
        "-p",
        "--params",
        required=True,
        help="parameters file name",
    )

    optional.add_argument(
        "--nBinsEnergy",
        type=int,
        default=10000,
        required=False,
        help="number of bins for the energy dimension (default: %(default)i)",
    )
    optional.add_argument(
        "--nBinsAngle",
        type=int,
        default=1800,
        required=False,
        help="number of bins for the angle dimension (default: %(default)i)",
    )

    args = parser.parse_args()

    inputFileNames = args.input
    outputFileName = args.output

    print(f"{len(inputFileNames)} input files:")
    print(f"inputFileNames : {inputFileNames}")
    print("")

    paramsFileName = args.params

    processingParams = ProcessingParams.fromYaml(paramsFileName)
    processingParams.print()

    print(f"Reading run infos...")
    runInfos = RunInfos.fromRootFiles(inputFileNames)

    runInfos.print()

    nDetectors = len(runInfos.detectorInfos)

    eResol = processingParams.eResol
    pResol = processingParams.pResol

    if len(eResol) != nDetectors:
        raise ValueError(
            f"ERROR: there are {nDetectors} detectors but {len(eResol)} energy resolutions"
        )

    # system matrix bins
    if args.nBinsEnergy:
        nBinsEnergy = args.nBinsEnergy
    if args.nBinsAngle:
        nBinsAngle = args.nBinsAngle

    histoDegOutputFileName = f"histo_{outputFileName}_DEG.root"
    histoCosOutputFileName = f"histo_{outputFileName}_COS.root"
    histoDegSmoothOutputFileName = f"histo_{outputFileName}_smooth_DEG.root"
    histoCosSmoothOutputFileName = f"histo_{outputFileName}_smooth_COS.root"
    kdeHistoCosOutputFileName = f"kde_histo_{outputFileName}_COS.root"
    

    sysMatDegOutputFileName = f"sysMat_{outputFileName}_DEG"
    sysMatCosOutputFileName = f"sysMat_{outputFileName}_COS"
    sysMatDegSmoothOutputFileName = f"sysMat_{outputFileName}_smooth_DEG"
    sysMatCosSmoothOutputFileName = f"sysMat_{outputFileName}_smooth_COS"
    kdeSysMatCosOutputFileName = f"kde_sysMat_{outputFileName}_COS"

    print("10 files will be created:")

    os.system(f"mkdir -p {dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos")

    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/{sysMatDegOutputFileName}.hdf5")
    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/{sysMatCosOutputFileName}.hdf5")
    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/{sysMatDegSmoothOutputFileName}.hdf5")
    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/{sysMatCosSmoothOutputFileName}.hdf5")
    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/{kdeSysMatCosOutputFileName}.hdf5")
    print("***************************************************************************")
    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{histoDegOutputFileName}")
    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{histoCosOutputFileName}")
    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{histoDegSmoothOutputFileName}")
    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{histoCosSmoothOutputFileName}")
    print(f"\t{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{kdeHistoCosOutputFileName}")
    print("*****************************************************************************")

    dataFrame = R.RDataFrame("tree", inputFileNames)

    R.RDF.Experimental.AddProgressBar(dataFrame)

    countInitFuture = dataFrame.Count()

    eventProcessor = EventProcessor(processingParams)

    initEnergy = runInfos.initEnergy
    maxEnergy = 0.91 * initEnergy

    dataFrame = eventProcessor.basicSelection(dataFrame)
    dataFrame = eventProcessor.topoCuts(dataFrame)
    dataFrame = eventProcessor.computePoints(dataFrame, runInfos.initPos)
    dataFrame = eventProcessor.distanceSelection(dataFrame)
    dataFrame = eventProcessor.computeEnergies(dataFrame, runInfos.initEnergy)

    dataFrame = dataFrame.Alias("eScat", "eScatCombined")
    dataFrame = dataFrame.Alias("eAbso", "eAbsoCombined")

    dataFrame = dataFrame.Filter(f"eScat<{maxEnergy}")

    dataFrame = eventProcessor.computePositions(dataFrame)
    dataFrame = eventProcessor.computeAngles(dataFrame)

    dataFrame = dataFrame.Alias("cosAngle", "cosAngleSmeared")
    dataFrame = dataFrame.Define("angleDeg", "angleSmeared*180/3.14159265")

    nEventsFinalFuture = dataFrame.Count()

    print("SYSTEM MATRIX: ")
    print(f"      Energy: {nBinsEnergy} bins from 0 to {maxEnergy:.2f} keV")
    print(f"       Angle: {nBinsAngle} bins\n")

    systemMatrixDeg = dataFrame.Histo2D(
        (
            "sysMatDeg",
            f";Scatter energy (keV); Scatter angle (deg)",
            nBinsEnergy,
            0,
            maxEnergy,
            nBinsAngle,
            0,
            180,
        ),
        "eScat",
        "angleDeg",
    )

    systemMatrixCos = dataFrame.Histo2D(
        (
            "sysMatCos",
            f";Scatter energy (keV); Scatter cos angle",
            nBinsEnergy,
            0,
            maxEnergy,
            nBinsAngle,
            -1,
            1,
        ),
        "eScat",
        "cosAngle",
    )

    print("Extracting data for KDE...")
    kde_start_time = time.time()
    # Get columns directly from RDataFrame as numpy arrays
    numpy_dict = dataFrame.AsNumpy(["eScat", "cosAngle"])
    kde_data = np.column_stack((numpy_dict["eScat"], numpy_dict["cosAngle"]))


    KDE_nBinsEnergy = nBinsEnergy
    KDE_nBinsAngle = nBinsAngle
    
    # Create KDE grid matching original histogram bins
    energy_grid = np.linspace(0, maxEnergy, KDE_nBinsEnergy)
    cosAngle_grid = np.linspace(-1, 1,KDE_nBinsAngle)
    X, Y = np.meshgrid(energy_grid, cosAngle_grid)
    grid_points = np.vstack([X.ravel(), Y.ravel()]).T


    print("Fitting KDE...")
    kde = KernelDensity(kernel="epanechnikov", bandwidth=0.1,algorithm="ball_tree", metric="euclidean")
    kde.fit(kde_data)
    print("Computing log density...")
    log_density = kde.score_samples(grid_points)
    Z = np.exp(log_density).reshape(X.shape)
    kde_end_time = time.time()
    print(f"KDE computation took {kde_end_time - kde_start_time:.2f} seconds\n")
    
    


    systemMatrixDeg.SetContour(99)
    systemMatrixCos.SetContour(99)

    systemMatrixDegSmooth = systemMatrixDeg.Clone("sysMatDegSmooth")
    systemMatrixCosSmooth = systemMatrixCos.Clone("sysMatCosSmooth")

    systemMatrixDegSmooth.Smooth(1)
    systemMatrixCosSmooth.Smooth(1)

    nEventsFinal = nEventsFinalFuture.GetValue()

    print(f"N events simulated      = {runInfos.nEvents}")
    print(f"N events kept in GEANT4 = {countInitFuture.GetValue()}")
    print(f"N events final          = {nEventsFinal}")

    sensitivity = nEventsFinal / runInfos.nEvents

    print(f"Sensitivity: {sensitivity:.2%}")


    print("\nProcessing KDE output into ROOT TH2D and HDF5...")
        
        # Calculate bin widths
    bin_width_eScat = maxEnergy / KDE_nBinsEnergy
    bin_width_cosAngle = 2.0 / KDE_nBinsAngle 

        
        # Z_kde has shape (nBinsAngle, nBinsEnergy)
    kde_counts_per_bin = Z * nEventsFinal * bin_width_eScat * bin_width_cosAngle #(Z_kde) to represent expected counts per bin
        
        # Create a ROOT TH2D for the KDE  data
       
    kdeRootHistCos = R.TH2D("kdeSysMatCos", "KDE System Matrix (cosAngle);Scatter energy (keV); Scatter cos angle",
                                KDE_nBinsEnergy, 0, maxEnergy,
                                KDE_nBinsAngle, -1, 1)

        # Fill the ROOT histogram bin by bin
        # Z_kde[j, i] corresponds to cosAngle_grid_kde[j] and energy_grid_kde[i]
    for i_energy in range(KDE_nBinsEnergy):      
        for j_cosAngle in range(KDE_nBinsAngle): 
            

            bin_content = kde_counts_per_bin[j_cosAngle, i_energy]
            kdeRootHistCos.SetBinContent(i_energy + 1, j_cosAngle + 1, bin_content)
        
    print(f"KDE-derived ROOT TH2D '{kdeRootHistCos.GetName()}' created and filled.")
        
    

     # Save the KDE into ROOT histogram
    kdeRootFile = R.TFile.Open(
            f"{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{kdeHistoCosOutputFileName}",
            "RECREATE" )
    kdeRootHistCos.Write("sysMat") 
    kdeRootFile.Close()
    print(f"KDE-derived ROOT histogram saved to: {dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{kdeHistoCosOutputFileName}")

    sysMatDegFile = R.TFile.Open(
        f"{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{histoDegOutputFileName}",
        "RECREATE",
    )
    systemMatrixDeg.Write("sysMat")
    sysMatDegFile.Close()

    sysMatCosFile = R.TFile.Open(
        f"{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{histoCosOutputFileName}",
        "RECREATE",
    )
    systemMatrixCos.Write("sysMat")
    sysMatCosFile.Close()

    sysMatDegFileSmooth = R.TFile.Open(
        f"{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{histoDegSmoothOutputFileName}",
        "RECREATE",
    )
    systemMatrixDegSmooth.Write("sysMat")
    sysMatDegFileSmooth.Close()

    sysMatCosFileSmooth = R.TFile.Open(
        f"{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{histoCosSmoothOutputFileName}",
        "RECREATE",
    )
    systemMatrixCosSmooth.Write("sysMat")
    sysMatCosFileSmooth.Close()

    def writeSysMatHDF5(histoFileName, outputFileName, angleMode):
        start_time = time.time()
        urFile = ur.open(f"{dirPathOfThisFile}/SYSTEM_MATRIXES_10/histos/{histoFileName}")
        urHist = urFile["sysMat"].to_numpy()

        data = urHist[0] / nEventsFinal

        outputFileName = f"{dirPathOfThisFile}/SYSTEM_MATRIXES_10/{outputFileName}.hdf5"
        with h5py.File(outputFileName, "w") as file:
            matrixData = file.create_dataset(
                "matrix", data=data, dtype="f", compression="gzip"
            )
            matrixData.attrs["angleMode"] = angleMode
            matrixData.attrs["maxEnergy"] = maxEnergy
            matrixData.attrs["normFactor"] = nEventsFinal

            file.attrs["sensitivity"] = sensitivity

            file.attrs["processingParams"] = str(processingParams.toDict())
            file.attrs["runInfos"] = str(runInfos.toDict())
        end_time = time.time()
        print(f"HDF5 file creation for {outputFileName} took {end_time - start_time:.2f} seconds")

    writeSysMatHDF5(histoDegOutputFileName, f"{sysMatDegOutputFileName}", "DEG")
    writeSysMatHDF5(histoCosOutputFileName, f"{sysMatCosOutputFileName}", "COS")
    writeSysMatHDF5(
        kdeHistoCosOutputFileName, f"{kdeSysMatCosOutputFileName}", "COS"
    )

    writeSysMatHDF5(
        histoDegSmoothOutputFileName, f"{sysMatDegSmoothOutputFileName}", "DEG"
    )
    writeSysMatHDF5(
        histoCosSmoothOutputFileName, f"{sysMatCosSmoothOutputFileName}", "COS"
    )
    
    script_end_time = time.time()
    print(f"\nTotal script execution time: {script_end_time - script_start_time:.2f} seconds")
