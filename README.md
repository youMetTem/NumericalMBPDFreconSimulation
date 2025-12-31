# Numerical Reconstruction of the Maxwell-Boltzmann Distribution via 3D Monatomic Hard-Sphere Simulation
This project explores the numerical reconstruction of the Maxwell-Boltzmann probability density function based on data generated from a 3D hard-sphere simulation system. Additionally, the study evaluates the accuracy of reconstruction by comparing cubic spline interpolation and standard non-linear curve fitting (`scipy.optimize`) against the theoretical curve.
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
6. [Acknowledgements & Resources](#acknowledgements--resources)
7. [AI Use Declaration](#ai-use-declaration)



## Theoretical Background

### The Maxwell-Boltzmann Distribution
For an ideal gas at thermodynamic equilibrium, the speeds of particles are not unitform. Instead, they follow a specific probability distribution known as *Maxwell-Boltzmann PDF*.

The probability density function $f(v)$ for a particle of mass $m$ at temperature $T$ is given by:

$$f(v)=4\pi \left( \frac{m}{2\pi k_B T} \right)^{3/2} v^2 \exp\left(-\frac{mv^2}{2k_B T}\right)$$

Conventional analytical approaches to derive the aforementioned function rely on frameworks from statistical mechanics or entropy maximization. The primary objective of this project is to numerically reconstruct the *Maxwell-Boltzmann probability density function (PDF)* from the data collected from a simulated 3D thermodynamic system using a hard-sphere particle model.

Beyond simple verification, this study aims to evaluate the accuracy of statistical reconstruction techniques. Specifically, *Cubic Spline Interpolation* against *Non-Linear Curve Fitting* (via `scipy.optimize`) to determine which method better reconstruct the theoretical distribution from discrete, stochastic simulation data.



## Methodology
This project consists of 2 phases: the Kinematic Simulation, which generate raw physical data, and the Numerical Reconstruction, which analyzes the statistical properties and implement *Cubic Spline Interpolation* and *Non-Linear Curve Fitting* of the system (`main.py`). Additionally, 3D particle simulation visualizer is also provided (`vis3D.py`)

### Physical Assumptions
I treat the system as an *Ideal Gas*. The require specific constraints on how particles behave:
1. *No Intermolecular Forces*: Particles do not attract or repel each other at a distance. They only interact when they physically collide (*Hard-Sphere Model*).
2. *Random Motion*: Particles move in straight lines in random directions until they collide with each other or the container wall.
3. *Elastic Collisions*: All collisions are perfectly elastic, no energy is lost to heat or deformation.
4. *Small Atomic Radius*: Particles posses very small atomic radius (matching the real size of an Element in nanometers).

Note: As the simulation needs to check for particle colliding, I can not assume of *Point Mass*.

### Simulation Environment
The simulation models an ideal gas within a bouned 3D cubic container of length $L$. The simulation is initilized with $N$ spherical particles, each having mass $m$ and radius $R$ at $T$ temperature in Kelvin. ($L$, $N$, $m$, $T$, and particle element can be modified in the `main.py` file.)
#### Initialization Conditions:
1. *Positions* ($\vec{r}$): Initialized uniformly apart from each other within the domain $(R, L, -R)$ for all dimension $(x, y, z)$. Ensuring no two particles overlapping when spawned.
2. *Velocityes* ($\vec{v}$): Initialized with random components in $(x, y, z)$ but all components add up to magnitude of root mean square speed ($v_{rms}$). Directly correlate Temperature ($T$) with the simulation's environment.
3. *Discrete Time Steps* ($dt$): Set as $0.2$ factor of the time interval particle takes to move with displacement equals to its radius. Preventing unexpected particle tunnelling from excessive initial $dt$.

### Kinematics Simulation
The simulation evolves over discrete time steps $dt$. With each step $dt$, particle position are updated using Forward Euler integration:

$$\vec{r}(t+dt) = \vec{r}(t) + \vec{v}(t)dt$$

### Collision Handling Mechanisms
Collision observed in this simulation consists of 2 distinct type: *Wall-Particle Collisions* and *Particle-Particle Collisions*. The simulation handles these interaction differently.
#### Wall-particle Collisions
The container walls are treated as infinite mass barriers with perfect elasticity. If particle component $r_{i}$ reaches or exceeds the boundary limits $L$, the velocity component perpendicular to the wall will be inverted:

$$\vec{v}_{\perp, new} = -\vec{v}_{\perp, old}$$

#### Particle-Particle Collisions
At every step of evolving $dt$, the scripts checks for overlapping particle pairs. By finding Euclidean distance between each particles pair $i$, $j$ and check if they are below particle's diameter:

$$\lVert \vec{r}_{i}-\vec{r}_{j} \rVert \le 2R$$

Once the collision event is triggered, the velocities are updated based on the convervation of linear momentum and kinetic energy. For two particles of equal mass, the post-collision state $v_{1, f}$, $v_{2, f}$ are calculated using vector projection along the line of impact: 

$$\vec{v}_{1, f} = \vec{v}_{1, i} - \frac{(\vec{v}_{1, i}-\vec{v}_{2, i}) \cdot (\vec{r}_{1, i}-\vec{r}_{2, i})}{\lVert \vec{r}_{1, i}-\vec{r}_{2, i} \rVert^2} (\vec{r}_{1, i}-\vec{r}_{2, i})$$
$$\vec{v}_{2, f} = \vec{v}_{2, i} - \frac{(\vec{v}_{2, i}-\vec{v}_{1, i}) \cdot (\vec{r}_{2, i}-\vec{r}_{1, i})}{\lVert \vec{r}_{2, i}-\vec{r}_{1, i} \rVert^2} (\vec{r}_{2, i}-\vec{r}_{1, i})$$

The vector form of velocites ensures that the changes in momentum is directly and only applied along the normal vector connecting the particles pair center.


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



## Acknowledgements & Resources
This project was inspired by serveral excellent resources. Special thanks to the following creators, authors for their high quality educational content:
* **Physics for Scientists and Engineers (Serway & Jewett)** - Primary reference for the Kinetic Theory of Gases and Maxwell-Boltzmann derivation.
* **[Numerical Methods for Engineers](https://youtube.com/playlist?list=PLkZjai-2Jcxn35XnijUtqqEg0Wi5Sn8ab&si=NUFmcl7sFmRZM9Wr)** by *Prof. Jeffrey Chasnov* - Provided the concept of numerical matrix operation and cubic spline interpolation.
* **[EGME206 Numerical Methods for Engineers, Spring 2021](https://youtube.com/playlist?list=PLLM1AZpDbYI1HPhHWmA-_5YqF26Rnt9ZM&si=xaMQ15EgJC063XTu)** by *Prof. Ittichote Chuckpaiwong* - In-dept explanation in the topic of Regressions.
* **[VPython for Beginners](https://youtube.com/playlist?list=PLdCdV2GBGyXOnMaPS1BgO7IOU_00ApuMo&si=8MC25WVnnu3wRcem)** by *Let's Code Physics* - Introduce the usage of `VPython` for 3D simulation.



## AI Use Declaration
This project utilized AI tools to assits in the development process and text refinement.
* **GitHub Copilot:** Used for code autocompletion, debugging scripts, and snipped generation.
* **Google Gemini:** Used for drafting the initial structure of this README file, text refinement, and documentation formatting.

All scientific logic, mathematical derivations, and final code implementation were verified and reviewed by me.
