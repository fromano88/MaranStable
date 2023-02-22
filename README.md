# MaranStable 3.1
 
How to use MaranStable 3.1?

There are two possibilities:
   1) Open Matlab and run the "main.m" file in the main directory. Note that
      MaranStable 3.1 is only supported by Windows and Linux.
   2) You can install the standalone application by running the executable
      "MS3.1_win.exe" on a windows machine or "MS3.1_linux.install" on a
      linux machine. You find both files in the "bin" folder. There, You also
      find additional readme_*.txt files that guide through the installation.

In both cases, the software requires Matlab Version R2022a or newer. All
source files required to run MaranStable 3.1 are stored in the "src" folder.


Load our tutorial cases from the folder "tutorials" and go through
"tutorials.pdf" from the "docs" folder to get started with MaranStable 3.1.
They are intended to guide you to the results shown in our paper
"MaranStable: A linear stability solver for multiphase flows in canonical
geometries".


Supplementary material:

The supplementary material is stored in the "docs" folder and contains:

   1) "benchmarks.pdf":
      Extensive verification and validation cases as well as a grid
      convergence study that proves the accuracy of our software.

   2) "discretization.pdf":
      Explains the numerical method and discretization used in MaranStable.

   3) "flowchart.pdf"

   4) "tutorials.pdf" (see above)

   5) "VDI.pdf", "Shin-Etsu04.pdf"
      All fluids implemented in MaranStable 3.1. are taken from those files.
      For the silicone oils, the manufacturers do not provide informations
      about the surface tension coefficient and the exponential decay of the
      kinematic viscosity. In this case, please take a look at:
      [1] I. Ueno et. al., Phys. Fluids 15, 408-416 (2003)
      [2] T. Yano, J. Phys.: Conf. Ser. 327 012029 (2011)
      [3] F. Romano et. al., Phys. Fluids 29, 093303 (2017)
