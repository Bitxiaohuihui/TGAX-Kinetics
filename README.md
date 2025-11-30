# TGAX-Kinetics

**A Comprehensive Python Framework for Non-Isothermal Kinetic Analysis & Lifetime Prediction**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## 📖 Overview

**TGAX-Kinetics** is an open-source, GUI-based software designed for the rigorous kinetic analysis of thermogravimetric (TGA) data. 

While originally optimized for metastable hydrogen storage materials (e.g., $\alpha$-AlH$_3$), it is a **general-purpose tool** suitable for pharmaceuticals, polymers, energetic materials, and other solid-state decomposition processes. The software bridges the gap between raw thermal data and practical engineering decisions by offering transparent algorithms for mechanism identification and shelf-life prediction.

This repository contains the source code associated with the paper:
> **TGAX Kinetics: An Open-Source Python Framework for Non-Isothermal Kinetic Analysis of Hydrogen Storage Materials**

## ✨ Key Features

* **User-Friendly GUI:** Built with Python's Tkinter, offering a responsive interface without complex installation requirements.
* **Advanced Isoconversional Methods:**
    * **Differential:** Friedman method.
    * **Integral:** KAS (Kissinger-Akahira-Sunose) and OFW (Ozawa-Flynn-Wall).
    * **Non-Linear:** Advanced Vyazovkin method (minimizes error without integral approximations).
* **Robust Global Fitting:** * Implements **Multi-variate Non-Linear Least Squares (NLLS)** optimization.
    * Includes an **Interquartile Range (IQR) filter** to exclude statistical outliers during fitting.
    * Supports complex mechanisms: Sestak-Berggren (SB), Kamal-Sourour (Autocatalytic), and Nucleation-Growth models.
* **Stochastic Lifetime Prediction:** * Calculates shelf-life at arbitrary temperatures.
    * Uniquely incorporates **Jacobian-based uncertainty propagation** to provide 95% Confidence Intervals (CI) for predictions.
* **Comprehensive Reporting:**
    * One-click export of analysis reports to **Microsoft Word (.docx)** and **LaTeX (.tex)**.
    * High-quality plots suitable for publication (Angewandte/ACS style).

## ⚙️ Installation

### Prerequisites
* Python 3.8 or higher.

### Step 1: Clone the Repository
```bash
git clone [https://github.com/](https://github.com/)[YourUsername]/TGAX-Kinetics.git
cd TGAX-Kinetics
