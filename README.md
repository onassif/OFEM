# OFEM: A Finite Elemnent Analysis Software Utilizing Object-Oriented Programming

**Author: Omar Nassif, Ph.D**

**For Inquiries, Feel Free To Contact Me At: omar@onassif.com**

# Table of Contents:
* [About OFEM](#About-OFEM)
    * [Current Functionality](#Current-Functionality)
    * [Current Material Models Supported](#Current-Material-Models-Supported) 
* [Installation](#Installation)
    * [Quick Start](#Quick-Start)
* [Elements Configuration](#Elements-Configuration) 
* [Writing input files](#Writing-input-files)
* [To do list](#To-do-list)


# About OFEM

OFEM is a solid mechanics software that started as a PhD project that took a form of itself.

The software was developed with an emphasis on readiness and modularity, something I feel is not common in FEA software.

It's a great tool if you're learning FEM or you want to check your results.


## Current Functionality
Currently All elements are standard FEA Elements (No reduced integration, incompatible modes, etc...)

**All 2D Elements are Plane-Strain**

### Current Material Models Supported
* Linear Elasticity (Infinitesimal Strain Theory).
  * Material-ID: 1
  * Supported Elements: 
    * 2D: T3,  T6, Q4, Q9
    * 3D: T4, T10, Q8
* Mixed Formulation (Displacement+Pressure) Linear Elasticity
  * Material-ID: 6
  * Supported Elements:
    * 2D: T3,  T6, Q4, Q9
* Hypo-Elasticity
  * Material-ID: 2
  * Supported Elements:
    * 2D: T3,  T6, Q4, Q9 
    * 3D: T4, T10, Q8
* Hyper-Elasticity (NeoHookean)
  * Material-ID: 3
  * Supported Elements:
    * 2D: T3,  T6, Q4, Q9 
    * 3D: T4, T10, Q8
* Rate-Independent Plasticity
  * Material-ID: 5
  * Supported Elements:
    * 2D: T3,  T6, Q4, Q9 
    * 3D: T4, T10, Q8
* Visco-Plasticity (Power-Law Isotropic Hardening)
  * Material-ID: 4
  * Supported Elements:
    * 2D: T3,  T6, Q4, Q9 
    * 3D: T4, T10, Q8
* Crystal Plasticity (MTS Hardening)
  * Material-ID: 11
  * Supported Elements:
    * 2D: Q4
    * 3D: Q8
* Linear Elastic Discontinuous Galerkin
  * Material-ID: 8
  * Supported Elements:
    * 2D: T3,  T6, Q4, Q9 
    * 3D: T4, T10, Q8
* Hyper Elastic (NeoHookean) Discontinuous Galerkin
  * Material-ID: 9
  * Supported Elements:
    * 2D: Q4
* Crystal Plasticity Discontinuous Galerkin
  * Material-ID: 10
  * Supported Elements:
    * 3D: Q8

# Installation

## Quick Start
Here are the very simple steps for running the bot on Windows, however most of these instructions should be followed
regardless of your OS:
1. [Turn on your computer](https://www.google.com/search?q=how+do+I+turn+on+my+computer)
2. [Install Matlab](https://www.mathworks.com/products/matlab.html). 
3. Download GitHub Desktop and Open the OFEM Repository with GitHub Desktop (or download the zip file). 
4. Open the OFEM folder in File Explorer. Double click main.m and run it.
5. The UI will show all files currently in the input folder, you can either run them or make your own input file and run it.
6. Wait for OFEM to finish running the simulation..
7. Another dialog would appear with results, the results are also saved in the history folder.


# Elements Configuration
TBW

# Writing input files
Please try it to figure it yourself for now. I believe in you!!

# To do list
* Fix documentation
* Convert to python
* Add shell elements
* Add reduced integration, incompatible modes, hourglass control, etc...
* Add parallization
* Make it compatible with abaqus input files
* Include Paraview
* C++ Computations
* Add more comments to code
* Web developement?
* X-FEM
