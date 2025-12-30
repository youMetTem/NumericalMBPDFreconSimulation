# Numerical Reconstruction of the Maxwell-Boltzmann Distribution via 3D Monatomic Hard-Sphere Simulation
This project explores the numerical reconstruction of the Maxwell-Boltzmann probability density function based on data generated from a 3D hard-sphere system. Additionally, the study evaluates the accuracy of reconstruction by comparing cubic spline interpolation and standard non-linear curve fitting (scipy.optimize) against the theoretical curve.
<table>
  <tr>
    <td width="40%" align="center" valign="middle">
      <img src="assets/forProjectOverview/03_dynamicPDFKDESimulation.gif" width="100%" />
    </td>
    <td width="20%" align="center" valign="middle">
      <img src="assets/forProjectOverview/06_3DParticlesSimulation.gif" width="100%" />
    </td>
    <td width="40%" align="center" valign="middle">
      <img src="assets/forProjectOverview/04_GenericFunctionCurveFitted.png" width="100%" />
    </td>
  </tr>
</table>

## Tabel of Contents
1. [Project Overview](#numerical-reconstruction-of-the-maxwell-boltzmann-distribution-via-3D-monatomic-hard-sphere-simulation)
2. [Theoretical Background](#theoretical-background)
3. [Methodology](#methodology)
4. [Results](#results)
5. [Installation & Usage](#installation--usage)
6. [AI Use Declaration](#ai-use-declaration)

## Theoretical Background
in process ...

## Methodology
in process ...

## Results
in process ...

## Installation & Usage

### 1. Prerequisites & Setup
First, ensure you have Python 3.x installed. For Window OS replace `python3` with `python`.

Clone repository:
```sh
git clone https://github.com/youMetTem/NumericalMBPDFreconSimulation.git
cd NumericalMBPDFreconSimulation
```
Install requirements:
```sh
python3 -m venv env
. env/bin/activate
pip install -r requirements.txt
```

### 2. Running the Simulation and Analysis (`main.py`)
By default, the simulation models 2000 Helium molecules at a temperature of 373 K confined within a cubic box of side length $10^{-9}\\,\text{m}$. These parameters can be modified in the `main.py` file.

Run the script:
```sh
python3 src/main.py
```

An interactive menu will appear in your terminal. Select a mode by typing the corresponding number to view visualization of the analysis.
The following is an example of the output:
```text
$ python3 src/main.py

=====Helium Molecules Simulation=====
------------------------------------------
Please Select Mode of Computation
------------------------------------------
[1] Visualizing Initial PDF Compared to Theoretical PDF
[2] Visualize Dynamic Simulated PDF Histogram
[3] Visualize Dynamic Simulated PDF (KDE curve)
[4] Visualize Estimated Function From Hypothesized Generic Function
[5] Visualize Histogram's Bin Midpoint Spline Interpolated Curve
[6] Exit
<*> Default mode is set to Helium, change this in main.py, parameters block
<*> Mode [4], [5] might take a few second
------------------------------------------
Selected Mode: 
```
After selecting a mode, wait a few seconds for the Matplotlib visualization window will appear after a brief delay. Once you are satisfied with viewing the graphs, close the window to select another mode within the same simulation batch. To exit the program, enter `6`

### 3. Running the 3D Simulation Visualization (`vis3D.py`)
The script lauches a real-time hard-sphere particle simulation of particles colliding within a box of fixed dimensions using `VPython`.

Run the script:
```bash
python3 src/vis3D.py
```

A local web server will start, and your default we browser will automatically open a new tab to display the 3D animation.

#### 3D Simulation Controls (MacOS/Trackpad)
Once the browser window opens, you can navigate the 3D space using these trackpad gestures:

| Action | Trackpad Gesture | Keyboard Alternative |
| :--- | :--- | :--- |
| **Rotate View** | **Two-finger Click** & Drag | Hold **Ctrl** + Click & Drag |
| **Zoom In/Out** | **Two-finger Swipe** (Up/Down) | Hold **Option (Alt)** + Click & Drag |
| **Pan (Move)** | Hold **Shift** + Click & Drag | -- |

Note: For more details on camera control, refer to the official [VPython Documentation](https://www.glowscript.org/docs/VPythonDocs/index.html#).

## AI Use Declaration
This project utilized AI tools to assits in the development process and text refinement.
* **GitHub Copilot:** Used for code autocompletion, debugging scripts, and snipped generation.
* **Google Gemini:** Used for drafting the initial structure of this README file, text refinement, and documentation formatting.

All scientific logic, mathematical derivations, and final code implementation were verified and reviewed by me.
