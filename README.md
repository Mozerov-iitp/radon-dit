# radon-dit
Discrete Inverse Radon Transform based on the Direct Integration Theorem
markdown
# Radon-DIT: Discrete Inverse Radon Transform

This repository contains the implementation of the **Discrete Inverse Radon Transform** based on the **Direct Integration Theorem**, as described in the paper:

> **"Название вашей статьи"** (Authors, Year) — *Under review / Published in Journal Name*

## Overview
In this work, we present a numerical method for the inverse Radon transform that utilizes direct integration instead of classical frequency-domain approaches (like FBP). 

### Key Features
* **Direct Integration Theorem:** Accurate analytical reconstruction.
* **Verification:** Scripts and images to verify the results presented in the paper.
* **Accuracy:** High-precision results for discrete data sets.

## Repository Structure
* **DIT_VERSUS_FBP/** — Main C++ project directory.
* **DIT_VERSUS_FBP/Direct__Integration__Mapping.sln** — Visual Studio 2022 Solution file.
* **DIT_BIN/** — Test images and input data (phantom and sinograms).

## Getting Started
1. Install **Visual Studio 2022**.
2. Open `DIT_VERSUS_FBP/Direct__Integration__Mapping.sln`.
3. Set the configuration to **Release** / **x64**.
4. Build and Run (F5).

## Citation
If you use this code in your research, please cite:
