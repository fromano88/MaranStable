# MaranStable 3.1

MaranStable was developed at TU Wien in the Institute of Fluid Mechanics and Heat Transfer under the guidance of Professor Hendrik C. Kuhlmann. Michael Lukasser developed MaranStable 1.0 in the framework of the project Engineering Marangoni Flows (EMA), which was supported by FFG. Mario Stojanovic took over in the framework of the project SAJE, which is the acronym of the FFG Project Stability Analysis for the JEREMI Experiment. He developed MaranStable 2.x and MaranStable 3.x with the help of Francesco Romanò, where they removed several bugs, extended the source code by many features, and created a graphical user interface (GUI) for MaranStable.
_________________

How to use MaranStable 3.1?

There are two possibilities:
   1) Open Matlab and run the "main.m" file in the main directory. Note that MaranStable 3.1 is only supported on Windows and Linux.
   2) You can install the stand-alone application by running the executable "MS3.1_win.exe" on a Windows machine or "MS3.1_linux.install" on a Linux machine. You find both files in the "bin" folder. There, you also find additional readme_*.txt files that guide through the installation.

In both cases, the software requires Matlab Version R2022a or newer. All source files required to run MaranStable 3.1 are stored in the "src" folder.
_________________

Load our tutorial cases from the folder "tutorials" and go through "tutorials.pdf" in the "docs" folder to get started with MaranStable 3.1. They are intended to guide you all the way to the results provided in the paper "MaranStable: A linear stability solver for multiphase flows in canonical geometries" by M. Stojanovic et al. (2023).
_________________

Supplementary material:

The supplementary material is stored in the "docs" folder and contains:

   1) "benchmarks.pdf":   
      In this file several test cases for verification and validation are explained, as well as a grid convergence study that demonstrating the accuracy of the software.

   2) "discretization.pdf":   
      Explains discretization of the differential equations and the numerical methods and used to solve the discretized equations.

   3) "flowchart.pdf":   
      Provides an overview on the main elements of MaranStable.
      
   4) "tutorials.pdf" (see above):   
      Is a guide on how to arrive at the results presented in the main paper.

   5) "VDI.pdf", "Shin-Etsu04.pdf":   
      These files provide the temperature dependence of all material parameters implemented in MaranStable 3.1. For the silicone oils, the manufacturers do not specify the temperature dependence of the surface tension coefficient and of the kinematic viscosity. For these quantities, please see
      
      [1] I. Ueno et. al., Phys. Fluids 15, 408-416 (2003)
      
      [2] T. Yano et al., J. Phys.: Conf. Ser. 327, 012029 (2011)
      
      [3] F. Romanò et. al., Phys. Fluids 29, 093303 (2017)
