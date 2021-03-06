# BnSSA-eN

> This package is a modification of the ELSEPA code by F. Salvat, A. Jablonski, and C. J. Powell  [_Comp. Phys. Comm. 165, 157 (2005)_](https://www.sciencedirect.com/science/article/abs/pii/S0010465504004795?via%3Dihub) to make predictions for the beam-normal single-spin asymmetry (_B<sub>n</sub>_) in elastic electron-nucleus scattering at GeV beam energies. Inelastic intermediate-state excitations of the target nucleus were included into the Coulomb problem in the form of an optical potential. Theoretical details are discussed in the paper by O. Koshchii, M. Gorchtein, X. Roca-Maza, and H. Spiesberger, which is available on arXiv: [arXiv:2102.11809](https://arxiv.org/abs/2102.11809). The ELSEPA source package is available on [Mendeley](https://mendeley.figshare.com/articles/dataset/elsepa_Dirac_partial-wave_calculation_of_elastic_scattering_of_electrons_and_positrons_by_atoms_positive_ions_and_molecules/11332520) and is  distributed under a [CC BY 4](https://creativecommons.org/licenses/by/4.0/) license. 


## Table of Contents

* [Setup](#setup)
* [Usage](#usage)
* [General Info](#general-information)
* [Project Status](#project-status)


## Setup

The distribution package consists of Fortran 77 programs (source code files) and numerical data files. There are 9 files in the distribution package:
1. ELSEPA-related files (5 files):
    1. 'readme.txt' -  the present file.
    2. 'Makefile' -  makefile.
    3. 'elsepa.f' - calculation of elastic scattering of electrons and positrons by nuclei. Dirac partial wave analysis for real and complex central potentials; high-energy factorizations.
    4. 'bnssa.f' -  the main program file for predicting B<sub>n</sub> in elastic electron-nucleus scattering at GeV beam energies.
    5. 'mycase.in' --  example of an input data file for 'bnssa'. Relevant input parameters that can be changed for obtaining various predictions for _B<sub>n</sub>_ are: IZ, MNUCL, MABS, ISOT, IFR, UNC, and EV. The meaning of these parameters is described inside the input file. Other parameters should not be changed, and they are marked as NC (no change) inside the input file.

2. Files (4 files) related to the parametrization of structure functions that are extracted from inelastic virtual photoabsorption cross sections on the proton:
    1. 'f1f2wd.f' 
    2. 'f2allm.f' 
    3. 'f2glob.f'
    4. 'r1998.f'

The last 4 files are external to the ELSEPA code. They contain the Bosted-Christy parametrization [Phys. Rev. C 81, 055213 (2010)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.81.055213) (arXiv: https://arxiv.org/abs/0712.3731) for the proton structure functions that are needed for the evaluation of the absorptive potential _V<sub>abs</sub>_. The absorptive potential is computed according to the framework outlined in our paper https://arxiv.org/abs/2102.11809.
           
The file 'f1f2wd.f' includes an implementation of the first derivatives of the photoabsorption cross sections _??<sub>L</sub>_ and _??<sub>T</sub>_ at _Q<sup>2</sup>_=0. These derivatives are needed for the calculation of the absorptive potential within our framework.


## Usage

1. Generate the executable binary file 'bnssa.exe' using your Fortran compiler. For example, when I use the gfortran compiler in Linux, this executable is generated by means of the following command:

       gfortran -fno-align-commons -fno-automatic -ffixed-line-length-none -std=legacy bnssa.f -o bnssa.exe `cernlib -safe -G Motif mathlib`
      Alternatively, one can generate the executable file using the `make` command. If your compiler is different from gfortran, the Makefile needs to be modified appropriately.

2. The  <code>\`cernlib -safe -G Motif mathlib\`</code> part is needed in order to execute subroutines that evaluate nucleon structure functions. One may need to install some additional libraries for this command to work properly. In my case (Ubuntu 16.04), I needed to run the following lines prior to the first execution of the command provided in part 1:
    * sudo apt-get install libmathlib2-dev
    * sudo apt-get install libx11-dev
    * sudo apt-get install libmotif-dev

3. Use/modify the input datafile 'mycase.in' for your problem. Note that the following nuclei are implemented in this package:
    * <sup>208</sup>Pb (_Z_ = 82, _IS_ = 0)
    * <sup>90</sup>Zr (_Z_ = 40, _IS_ = 0)
    * <sup>48</sup>Ca (_Z_ = 20, _IS_ = 1)
    * <sup>40</sup>Zr (_Z_ = 20, _IS_ = 0)
    * <sup>28</sup>Si (_Z_ = 14, _IS_ = 0)
    * <sup>27</sup>Al (_Z_ = 13, _IS_ = 0)
    * <sup>12</sup>C (_Z_ = 6, _IS_ = 0)
    * <sup>4</sup>He (_Z_ = 2, _IS_ = 0)

4. Execute the code by entering the command:
            
    `./bnssa.exe < mycase.in` 
      
   which redirects the standard input unit (=5) to your input data file.


## General Information
The program 'bnssa.f' reads data from the input file 'mycase.in', whose structure is described in the heading comments of the Fortran 77 'bnssa.f' source file. A short description of the input parameters is provided in the file 'mycase.in'. The program uses external subroutines and functions contained in modules 'elsepa.f', 'f1f2wd.f', 'f2allm.f', 'f2glob.f', 'r1998.f'; these modules have been inserted into the source files of the main programs by means of an 'Include' statement.

For electron scattering from a target nucleus with with the atomic number _Z_ and isotop indicator _IS_ at a fixed beam energy, the calculated differential cross section (DCS) and beam-normal single-spin asymmetry _B<sub>n</sub>_ (also known as the Sherman function) are written to an output file in a format ready for visualization with a plotting program.
                    
When the absorptive potential _V<sub>abs</sub>_ is turned on (MABS=1), the output filename is given in the following format:\
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   'dcs_BPx1_Zx2_ISx3_UNCx4_x5_mx6.dat', where\
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    x1 - value of the Compton slope parameter used for the calculation of V<sub>abs</sub>,\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    x2 - input value for the parameter IZ (atomic number)\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    x3 - input value for the parameter ISOT (isotope indicator number)\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    x4 - input value for the parameter UNC (model for calculation of V<sub>abs</sub>)\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    x5 - input value for the parameter EV (beam energy in eV, in the E form format: x.yyyezz)\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    x6 - input value for the parameter MNUCL (nuclear charge density distribution model)
    
When the absorptive potential is turned off, the output filename is given in the following format:\
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    'dcs_NoVabs_Zx2_ISx3_x5_mx6.dat'\
\
Several subroutines write partial information (phase shifts, absorptive potential, ...) on the standard output unit (=6) as the calculation progresses; they may also generate files with the extension '*.dat' with self-explanatory contents. This does not impair the calculation speed and helps the user to get a  feeling of the time that a planned calculation may take and to track partial results of the calculation. 
   
Typical running times with the parameter MABS = 0 (V<sub>abs</sub> = 0) are between few seconds and several minutes, depending on the atomic number of the target and the beam energy. When MABS is set to 1 in the input file, the absorptive potential is evaluated using the formulas provided in our paper https://arxiv.org/abs/2102.11809 (choose UNC = 0 for V<sub>abs</sub> = V<sup>(0)</sup><sub>abs</sub> or choose UNC = 3 for V<sub>abs</sub> = V<sup>(0)</sup><sub>abs</sub> + V<sup>(1)</sup><sub>abs</sub>). Typical running times in this case are 5-30 hours depending on the atomic number of the target and the kinetic energy of the projectile.


## Project Status
Project is: _Completed_. 


<!---## Room for Improvement
Room for improvement:
- Improvement to be done 1
- Improvement to be done 2
To do:
- Feature to be added 1
- Feature to be added 2--->

