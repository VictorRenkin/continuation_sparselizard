# Continuation Sparselizard

### Academic Year 2024 – 2025

#### Author: Victor Renkin s2306326

## Table of Contents

* [Introduction](#introduction)

* [Requirements](#requirements)

* [Usage](#usage)

* [Project Structure](#project-structure)

* [AI Integration](#ai-integration)

## Introduction

This project focuses on implementing a **parallel Harmonic Balance Method** within the *Sparselizard* framework for the analysis of a nonlinear mechanical system. The primary objective is to accurately trace the system's **nonlinear frequency responses (NLFR)**. This is achieved through the application of a **continuation method** based on a robust predictor-corrector scheme.

## Requirements

This repository relies on **Sparselizard**. Please refer to the official documentation at [this link](https://github.com/halbux/sparselizard-users/tree/main/api/python) for proper installation instructions. Additionally, this project requires specific Python packages, which are listed in the `requirements.txt` file. To install these dependencies, ensure that both [Python](https://www.python.org/) and [pip](https://pip.pypa.io/en/stable/) are installed on your system, then execute the following command:

```bash
pip install -r requirements.txt
```

## Usage

For all computations, it is essential to precisely define the **physical regions** of the system. This includes the entire volume, any clamped surfaces, and the specific point (or surface) of excitation.

### Nonlinear Frequency Response (NLFR) Computation

To compute the nonlinear frequency response of the system, execute the following command:

```bash
python3 src/main_NLFR.py
```
### Nonlinear Normal Mode (NNM) / Frequency Response Function (FRF) Computation

To compute the Nonlinear Normal Mode (NNM) run:
```bash
python3 src/main_NNM.py
```

## Project Structure

The codebase is predominantly organized around a **class-based architecture**, distinguishing between the **NLFR** and **NNM** functionalities. Each of these main functionalities is encapsulated within a class that incorporates **abstract methods**. These abstract methods serve as guidelines for implementing custom predictor and corrector schemes, allowing for flexible extension of the framework.

The core of the code follows the **continuation loop**, which is represented by the `continuation_loop_NNM` and `continuation_loop_NLFR` functions. The adaptive adjustment of the predictor step size is managed by the `StepSizeRule` function, enabling efficient convergence. Furthermore, a dedicated function is implemented to **store previous solutions**, typically maintaining a history equivalent to the predictor's order plus one, which is crucial for the continuation process.

## AI Integration

Artificial intelligence tools are utilized periodically throughout this project. Their primary roles include **code correction** and **rephrasing sentences** within the accompanying report to enhance clarity and conciseness.






