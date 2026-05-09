# Radon-DIT: Discrete Solutions of the Inverse Radon Problem

This repository contains the source code and datasets for the following paper:

> **"The Direct Integration Theorem: A Rigorous Framework for Consistent Discrete Solutions of the Inverse Radon Problem"** (M. G. Mozerov). 
> *Currently under review in IEEE Transactions on Pattern Analysis and Machine Intelligence (Submitted May 2026).*


## Overview
This paper presents a novel Direct Integration Theorem (DIT), derived as a non-trivial corollary of the classical Central Slice Theorem (CST). The DIT provides a mathematically consistent transition from the continuous to the discrete domain—a fundamental challenge in computed tomography—thereby eliminating the need for frequency-domain interpolation without resorting to conventional ramp-filtering. 

The proposed approach circumvents two principal limitations inherent in traditional methods: 
1. The zero-frequency singularity and spectral distortions introduced by the mandatory ramp-filtering step.
2. Discretization inaccuracies associated with frequency-domain interpolation.

Mathematical modeling demonstrates that this approach achieves quasi-exact reconstruction, with errors constrained solely by sampling parameters and grid geometry. Comparative simulations confirm that the proposed method eliminates common artifacts, such as intensity cupping, and consistently outperforms FBP in terms of PSNR, SSIM, and reprojection fidelity.

## Repository Structure
* `DIT_VERSUS_FBP/` — C++ source code and MS Visual Studio 2022 project files.
* `DIT_BIN/` — Dataset of 10 test images (512x512) used in the study.
## Instructions
1. **Data:** Place the `DIT_BIN` folder in the root of any local drive (e.g., `D:\DIT_BIN\`).
2. **Build:** Open `Direct__Integration__Mapping.sln` in Visual Studio 2022.
3. **Run:** Build in **Release / x64** mode. 
    * **Note:** Upon execution, the interactive console will ask for the drive letter (C, D, E, or F) where the `DIT_BIN` folder is located. Ensure the folder is in the root directory for the program to find the test images.
## Citation
If you find this work useful for your research, please cite the corresponding paper (refer to the manuscript or arXiv for the latest citation format).
