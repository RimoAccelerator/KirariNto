# KirariNto
**KirariNto** is a program that traces the electronic state during TD-DFT geometry optimization, by comparing the NTO overlap between each state with the target state of the previous geometry. It is based on a modified algorithm of Marco Campetella and Juan Sanz Garcia (10.1002/jcc.26162). Gaussian and ORCA are currently supported.

![image](https://github.com/RimoAccelerator/KirariNto/blob/main/kirarinto_logo.png)

# Requirements
Linux operation system

**Python modules:**

  PySCF
  
  Geometri
  
  PyBerny
  
  Numpy

**Programs:**

Multiwfn

Gaussian or ORCA

# Usage

1. Set the following enviromental variants:
   KIRARINTO_MULTIWFN
   
    KIRARINTO_GAUSSIAN
   
    KIRARINTO_FORMCHK
   
    KIRARINTO_ORCA
   
    KIRARINTO_ORCA2MKL

   Among them, KIRARINTO_MULTIWFN is required, and set how to invoke Gaussian or ORCA, along with the formchk or orca_2mkl tools according to your task.

   Example:

```export KIRARINTO_MULTIWFN=~/Multiwfn_3.8_dev_bin_Linux_noGUI/Multiwfn

export KIRARINTO_GAUSSIAN=~/g16/g16

export KIRARINTO_FORMCHK=~/g16/formchk

export KIRARINTO_MULTIWFN=~/Multiwfn_3.8_dev_bin_Linux_noGUI/Multiwfn_noGUI

export KIRARINTO_ORCA=/opt/orca5/orca

export KIRARINTO_ORCA2MKL=/opt/orca5/orca_2mkl
```

2. Prepare an input file for a TD-DFT single point calculation of your initial geometry. Number of states and the interested state should be set explicitly. Example:

(Gaussian)

```
%mem=8GB
%nprocshared=8
%chk=u_td_1.chk
# td=(nstates=3,root=1) b3lyp sto-3g em=gd3bj iop(9/40=4)

TC

0 1
 C                  1.49653099   -0.03970163    0.25798558
 C                  0.96661422    1.10990407   -0.32719382
 C                 -0.29560541    1.04203217   -0.91629061
 H                  2.46280152   -0.01555637    0.71693905
 H                  1.51812645    2.02681494   -0.32424024
 C                 -0.42910744   -1.21052831   -0.33465934
 O                 -1.07773001   -2.28888230   -0.33812694
 N                 -0.96370843   -0.11703076   -0.90600943
 N                  0.78578354   -1.17302439    0.24032710
 H                  1.16155785   -1.99994675    0.65864886
 N                 -0.88665177    2.23323444   -1.54286970
 H                 -1.88298815    2.18859514   -1.46992353
 H                 -0.55480592    3.05391114   -1.07771587
```

(ORCA)

```! blyp def2-sv(p) def2/j
%pal nprocs 8 end
%tddft nroots 5 iroot 1 end
* xyz 0 1
 C                  1.49653098   -0.03970163    0.25798558
 C                  0.96661422    1.10990407   -0.32719382
 C                 -0.29560541    1.04203217   -0.91629061
 H                  2.46280151   -0.01555637    0.71693905
 H                  1.51812644    2.02681493   -0.32424024
 C                 -0.42910744   -1.21052830   -0.33465934
 O                 -1.07773001   -2.28888229   -0.33812694
 N                 -0.96370843   -0.11703076   -0.90600943
 N                  0.78578354   -1.17302438    0.24032710
 H                  1.16155784   -1.99994674    0.65864886
 N                 -0.88665177    2.23323443   -1.54286969
 H                 -1.88298814    2.18859513   -1.46992352
 H                 -0.55480592    3.05391113   -1.07771587
*
```

3. Run the script:

```
usage: kirariNto.py [-h] [-c CONSTRAINTS] [-p PROGRAM] [-d DELETE] [-g GRID] input

positional arguments:
  input                 An input file for the TD-DFT calculation of your initial
                        geometry.

optional arguments:
  -h, --help            show this help message and exit
  -c CONSTRAINTS, --constraints CONSTRAINTS
                        (optional) An external geometry constraint file that follows the
                        format requirement of PySCF.
  -p PROGRAM, --program PROGRAM
                        The program you are using (Gaussian or ORCA; case insensitive).
                        Default: Gaussian if input file contains .gjf, and ORCA if input
                        file contains .inp
  -d DELETE, --delete DELETE
                        (optional) If the .chk or .gbw files to be deleted after each
                        step. yes (default) or no.
  -g GRID, --grid GRID  (optional) 1, 2 (default) or 3. The accuracy of integral grid for
                        overlap evaluation. Denser for a large value.
```

# How does it work?

**KirariNto** acts as an interface between Gaussian/ORCA and PySCF. It deals with the energies and gradients outputted by Gaussian/ORCA, and inputs them into the **geomeTRIC** optimizer of PySCF.

In each cycle unless the first one, a TD-DFT single point calculation is first performed. Then KirariNto invokes **Multiwfn**, to calculate the **Natural Transition Orbital** (NTO) of each excited state. The overlap integrals between each state and the previous target state for the electron and hole are calculated separately, summed up, and ordered. The state with the largest overlap integral is determined to be the state of interest for the next step.

The calculation of overlap integral between state A and B is achieved by the following process. This process is designed to fit in the available functions of Multiwfn.

1. Get the geometries of A and B aligned.
2. Perform NTO analysis for A and B separately, to obtain the top 10 NTO pairs and their eigenvalues. Save the NTO orbitals to a temporary wavefunction file.
3. Set the occupation number of each NTO orbitals to be the 4th root of its NTO eigenvalue, and store the orbitals of the hole and electron separately into 4 new wavefunction files.
4. Calculate the electron densitiy based on the 4 wavefunction files: 2 holes, and 2 electrons.
5. Calculate the overlap integral between the two hole densities, and the two electron densities. Sum them up to give the final overlap integral.
