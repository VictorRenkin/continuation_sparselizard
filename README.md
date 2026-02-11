# Continuation Sparselizard

## Table of Contents

* [Introduction](#introduction)

* [Requirements](#requirements)

* [Usage](#usage)

* [Project Structure](#project-structure)

* [AI Integration](#ai-integration)

## Introduction

This project is dedicated to the implementation of a **parallel Harmonic Balance Method** within the *Sparselizard* framework for the analysis of a **nonlinear mechanical system**. The main objective is the accurate computation and tracing of the system’s **nonlinear frequency responses (NLFRs)** and **nonlinear normal modes (NNMs)**. This is achieved through the use of a **continuation strategy** relying on a robust **predictor–corrector scheme**, ensuring numerical stability and efficiency when following solution branches.

The *Sparselizard* framework is used as a **development and prototyping environment**, enabling rapid implementation and validation of the proposed numerical formulations. Once validated, all developed algorithms have been **integrated into Quanscient Allsolve** (https://quanscient.com/), in order to benefit from advanced **high-performance computing capabilities**, such as **domain decomposition methods** and large-scale parallelisation. For clarity and reproducibility, the parts of the code transcribed from the Quanscient implementation are explicitly provided in the present project under the **`transcription_quanscient`** directory.

In addition, all numerical results presented and analysed in this work are fully reproducible using the **test cases** included in the **`testcase`** directory. These test cases cover a representative set of benchmark configurations, including a **clamped–clamped beam**, a **cantilever beam**, a **fan blade**, and a **cantilever beam in rotation**, thereby ensuring transparency and consistency between the numerical implementations and the reported results.

Finally, the governing formulation is expressed in a **rotating reference frame**, allowing **centrifugal effects** to be consistently taken into account in the system dynamics, which is essential for applications involving rotating mechanical structures.

## Requirements

This repository relies on *Sparselizard*. Please refer to the official documentation at [this link](https://github.com/halbux/sparselizard-users/tree/main/api/python) for proper installation instructions. Additionally, this project requires specific Python packages, which are listed in the `requirements.txt` file. To install these dependencies, ensure that both [Python](https://www.python.org/) and [pip](https://pip.pypa.io/en/stable/) are installed on your system, then execute the following command:

```bash
pip install -r requirements.txt
```

## Usage

The physical regions of the system must be clearly defined, including the full volume, clamped surfaces, excitation and measurement regions, together with the relevant material properties. For rotational effects, the distance to the axis of rotation must also be specified.

### Nonlinear Frequency Response (NLFR) Computation

To compute the nonlinear frequency response of the system, execute the following command:

```bash
python3 src/main_NLFR.py
```
### Nonlinear Normal Mode (NNM) Computation

To compute the Nonlinear Normal Mode (NNM) run:
```bash
python3 src/main_NNM.py
```

## Project Structure

The codebase is predominantly organized around a **class-based architecture**, distinguishing between the **NLFR** and **NNM** functionalities. Each of these main functionalities is encapsulated within a class that incorporates **abstract methods**. These abstract methods serve as guidelines for implementing custom predictor and corrector schemes, allowing for flexible extension of the framework.

The core of the code follows the **continuation loop**, which is represented by the `continuation_loop_NNM` and `continuation_loop_NLFR` functions. The adaptive adjustment of the predictor step size is managed by the `StepSizeRule` function, enabling efficient convergence. Furthermore, a dedicated function is implemented to **store previous solutions**, typically maintaining a history equivalent to the predictor's order plus one, which is crucial for the continuation process.

## AI Integration

Artificial intelligence tools are utilized periodically throughout this project. Their primary roles include **code correction** and **rephrasing sentences** within the accompanying report to enhance clarity and conciseness.






