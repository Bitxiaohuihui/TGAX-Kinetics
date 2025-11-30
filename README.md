# TGAX-Kinetics

**A Comprehensive Python Framework for Non-Isothermal Kinetic Analysis & Lifetime Prediction**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![DOI](https://zenodo.org/badge/1107020370.svg)](https://doi.org/10.5281/zenodo.17769287)

## 📖 Overview

**TGAX-Kinetics** is an open-source, GUI-based software designed for the rigorous kinetic analysis of thermogravimetric (TGA) data.

While originally developed for the safety assessment of metastable hydrogen storage materials (e.g., $\alpha$-AlH$_3$), TGAX-Kinetics is a **general-purpose tool** widely applicable to pharmaceuticals, polymers, energetic materials, and other solid-state decomposition processes. It bridges the gap between raw thermal data and practical engineering decisions by offering transparent algorithms for mechanism identification, robust global fitting, and stochastic shelf-life prediction.

This repository contains the source code associated with the research paper:
> **TGAX Kinetics: An Open-Source Python Framework for Non-Isothermal Kinetic Analysis of Hydrogen Storage Materials**

## ✨ Key Features

### 1. User-Friendly Interface
* Built with Python's **Tkinter**, offering a responsive graphical user interface (GUI) that requires no programming knowledge to operate.
* Interactive plots with "Pop-out" capabilities for high-resolution export.

### 2. Advanced Isoconversional Methods (Model-Free)
Determines the activation energy ($E_a$) as a function of conversion ($\alpha$) using ICTAC-recommended algorithms:
* **Differential:** Friedman method.
* **Integral:** KAS (Kissinger-Akahira-Sunose) and OFW (Ozawa-Flynn-Wall).
* **Non-Linear:** **Advanced Vyazovkin method**, which minimizes the error function without relying on approximate integral solutions.

### 3. Robust Multi-Variate Global Fitting
* **Complex Mechanism Support:** Fits experimental data to advanced models including:
    * Sestak-Berggren (SB)
    * Kamal-Sourour (Autocatalytic/Curing)
    * General Autocatalytic Initiation (GAI)
    * Parallel Reaction Pathways
* **IQR Filter:** Implements an Interquartile Range (IQR) filter to statistically exclude outliers in $E_a$ during global fitting, ensuring parameters are derived from the most physically relevant data segments.
* **Cross-Validation:** Includes Leave-One-Out Cross-Validation (LOOCV) to prevent overfitting.

### 4. Stochastic Lifetime Prediction
* Predicts conversion time ($t_{\alpha}$) at arbitrary isothermal storage temperatures.
* **Uncertainty Propagation:** Uniquely incorporates **Jacobian-based uncertainty propagation** to provide rigorous **95% Confidence Intervals (CI)** for all lifetime predictions, essential for safety assessments.

### 5. Automated Reporting
* **One-Click Export:** Generates comprehensive analysis reports in **Microsoft Word (.docx)** and **LaTeX (.tex)** formats.
* **Publication-Ready Plots:** Exports charts in high-DPI formats (PNG, TIFF, SVG) conforming to journal standards (e.g., *Angewandte*, *ACS*).

## ⚙️ Installation

### Prerequisites
* Python 3.8 or higher.

### Option 1: Quick Install (Recommended)
1.  Clone this repository:
    ```bash
    git clone [https://github.com/Bitxiaohuihui/TGAX-Kinetics.git](https://github.com/Bitxiaohuihui/TGAX-Kinetics.git)
    cd TGAX-Kinetics
    ```
2.  Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

### Option 2: Manual Install
Ensure you have the following libraries installed:
* `numpy`
* `pandas`
* `matplotlib`
* `scipy`
* `python-docx`
* `openpyxl`

*(Note: `tkinter` is included with standard Python installations. Linux users may need to run `sudo apt-get install python3-tk`)*

## 🚀 Usage Guide

1.  **Launch the Software:**
    Run the main script from your terminal or IDE:
    ```bash
    python TGAX_Kinetics.py
    ```

2.  **Import Data:**
    * Click **"Import Data Files"**.
    * Select multiple TGA files (supports `.csv`, `.xlsx`, `.txt`).
    * *Format:* Files should contain columns for Time, Temperature, and Mass (or TG%). The software includes an intelligent column parser to auto-detect headers.

3.  **Kinetic Analysis Workflow:**
    * **Isoconversional Analysis:** Choose a method (e.g., Vyazovkin) and set the conversion range (e.g., 0.05 - 0.95).
    * **Model Fitting:** Navigate to the `Analysis` menu to perform CKA (Combined Kinetic Analysis) or Autocatalytic Model Fitting.
    * **Prediction:** Use the results to predict shelf-life at ambient or specific temperatures.

4.  **Export Results:**
    * Use the `Export` menu to save data tables, plots, or generate a full textual report.

## 📂 Repository Structure

* `TGAX_Kinetics.py`: The main application source code.
* `requirements.txt`: List of Python dependencies.
* `BIT_Kinetics_Icon_Tight.ico`: Application icon.
* `splash.png`: Software splash screen image.

## 📄 Citation

If you use TGAX-Kinetics in your research, please cite the following paper:

> **Xing, X., Li, S., & Liu, J.** (2025). TGAX Kinetics: An Open-Source Python Framework for Non-Isothermal Kinetic Analysis of Hydrogen Storage Materials. *[Submitted/Journal Name]*.

[![DOI](https://zenodo.org/badge/1107020370.svg)](https://doi.org/10.5281/zenodo.17769287)

## 📜 License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details. This permissive license allows for reuse, modification, and distribution for both academic and commercial purposes, provided the original copyright notice is retained.

---

*Developed at the School of Materials Science and Engineering, Beijing Institute of Technology.*

