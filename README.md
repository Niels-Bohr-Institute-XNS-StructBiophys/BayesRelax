# BayesRelax

## Introduction

Software for estimating a relaxation spectrum from rheological measurements usually obtained by a frequency sweep analysis.

The mathematical underpinnings are presented in:

> Hansen SL (2008). Estimation of the relaxation spectrum from dynamic experiments using Bayesian analysis and a new regularization constraint. Rheologica Acta 47: 169-178

The first example of use in a biological context was published here:

> Hansen SL, Ray PM, Karlsson AO, Jørgensen B, Borkhardt B, Petersen BL, Ulvskov P (2011). Mechanical properties of plant cell walls probed by relaxation spectra. Plant Physiol 155: 246-258

The latter paper employed a web-server hosted implementation of the algorithm. This web-server no longer exists so here we provide a method for running the analysis locally from the command line.

BayesRelax is a Perl script that ties together three executables written in Fortran by Steen Laugesen Hansen.

## Systems requirements

The present implementation will compile and run on Linux and MacOS systems. MacOS 10.6.8 (with .eps for graphics output rather than pdf) and 10.13.6 have been tested as have Linux Mint version 19.3 both yielding plots in pdf.

Perl 5 is required and the script make use of the File::Basename and File::Copy modules

Gnuplot 5 is required (in PATH). This implementation was developed using version 5.4.

The gfortran compiler is required. On Linux it is installed via synaptic. Fortran is not included by default in the developer tools on newer Mac systems. It may conveniently be installed using homebrew but some libraries do not end up where the complier looks for them. This is solved as follows:

Create a folder named gfortran in /usr/local/ and a folder called lib inside the gfortran folder. Then create two aliases as follows:

    ln /usr/local/Cellar/gcc/10.2.0/lib/gcc/10/libgfortran.5.dylib /usr/local/gfortran/lib/libgfortran.5.dylib
    ln /usr/local/Cellar/gcc/10.2.0/lib/gcc/10/libquadmath.0.dylib /usr/local/gfortran/lib/libquadmath.0.dylib
    
## Compiling

    gfortran error2_cl.f -o error2
    gfortran conv_cl.f -o conv
    gfortran bayes_cl.f -O1 -o bayes

error2_cl.f and bayes_cl.f have beed edited relative to the original versions destined for the web-server. conv_cl.f is identical to the original conv.f.

## Running

The perl script and the three Fortran executables must reside in the same folder, so create a folder and place BayesRelax.pl, error2, bayes and conv there.

Navigate your terminal to that folder and:

    ./BayesRelax.pl <pathTo/inputDataFile.txt> <desired number of points in spectrum>

inputDataFile.txt is a tab-delimited text file with three columns: Frequency, storage modulus (G') and loss modulus (G''). Replicates need not be sorted nor balanced.

“desired number of points in spectrum” is an integer between 30 and 100.

A folder called Results will be generated and if a solution is found two plots will be saved plus the corresponding tab-delimited data files.

## Command line options

The present implementation only allows for setting the ‘number of points optional parameter’ as the defaults of the other ones are usually adequate. However, the script is easily modified to accommodate the remaining ones:

| Parameter | Description                                                                                       | Default           |
| --------- | ------------------------------------------------------------------------------------------------- | ----------------- |
| $qmin     | First point used in data file.                                                                    | All points used   |
| $qmax     | Last point used in data file.                                                                     | All points used   |
| $alphamin | Lower limit for the Lagrange multiplier, minimum is -20.                                          | -10               |
| $alphamax | Upper limit for the Lagrange multiplier, maximum is 20.                                           | 10                |
| $ncalc    | Number of steps between alphamin and alphamax (more steps increase the cpu-time). Maximum is 100. | 30                |
| $npoints  | Number of points in the relaxation spectrum (more points increase the cpu-time). Maximum is 100.  | 30                |
| $nfree    | Y' for free end points.                                                                           | H(End points) = 0 |

The bayes executable reads a 7 line text file with empty lines meaning ‘use default’. This file is written in line 31 of the script with ‘desired number of points in spectrum’ in line 6 of the output and newlines otherwise. Line 31 and line 11 of the the script should be edited to allow for other parameters to be set.
