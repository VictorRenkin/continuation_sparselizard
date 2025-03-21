# Continuation Sparselizard
![ViewCount](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https://github.com/VictorRenkin/continuation_sparselizard&count_bg=%2379C83D&title_bg=%23555555&icon=github.svg&icon_color=%23E7E7E7&title=views&edge_flat=false)
### Academic Year 2024 – 2025

#### Author: Victor Renkin s2306326

## Table of Contents
- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [AI](#ai)

## Introduction

## Requirements  
This repository works with **Sparselizard**. Follow the documentation available at [this link](https://github.com/halbux/sparselizard-users/tree/main/api/python) to install it properly. This project requires the packages listed in the `requirements.txt` file. To install these dependencies, ensure you have [Python](https://www.python.org/) and [pip](https://pip.pypa.io/en/stable/) installed, then run the following command:

```bash
pip install -r requirements.txt
```

## Usage

### Eigenfrequencies and Mode Shapes
To compute the eigenfrequencies and mode shapes, run:
```bash
python3 src/frequence_propre&Mode.py
```
The frequence will be find on the treminal and the mode shapes will be saved in the `data/` directory.

### Nonlinear Frequency Responses
To compute the nonlinear frequency response (NLFR), run:
```bash
python3 src/main.py
```
### Compute the FRF
To compute the FRF, run:
```bash
python3 src/linear_FRF.py
```

## Project Structure

### **`data/`**
Directory for storing generated data.

### **`geo_GMSH/`**
Contains mesh files used in the analysis.

### **`Quanscient/`**
Scripts for Quanscient computations related to NLFR.

### **`backbone-curve-nonlinear-harmonic-balance/`**
C++ implementation for:
- Clamped-clamped membrane with geometric nonlinearity.
- Frequency sweep around the resonance peak using the multiharmonic (harmonic balance) FEM.
- Example kindly provided and validated by S. Saxena (based on the paper: *Parallel Harmonic Balance Method for Analysis of Nonlinear Dynamical Systems*).

### **`figures/`**
Contains visual representations and figures related to the project.

### **`src/`**
Source code directory, organized as follows:
- **`frequence_propre&Mode.py`**: Computes eigenfrequencies and mode shapes.
- **`import.py`**: Handles all necessary imports for the project.
- **`requirements.txt`**: Lists required dependencies for the project.
- **`sparselizard_extension.py`**: Contains utility functions built using SparseLizard.
- **`VizData.py`**: Generates visualizations and figures.
- **`main.py`**: Computes nonlinear frequency responses (NLFR).

## AI
The AI is used to occasionally correct the code and to reformulate sentences from the report.