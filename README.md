# GalactICS galaxy model

This is a C++ reimplementation and extension of the models developed by Kuijken
& Dubinski (1995) and further modified by Widrow, Pym and Dubinski (2008). Those
models included a DF-based halo, bulge and single disk. The current model 
extends those models to include multiple collisionless disks, as well as 
optionally allowing for an Einasto halo profile instead of a cuspy NFW-type 
profile. The code has been cleaned up (no more GOTOs), comments added and cruft
removed, with a few bug fixes and optimisations. The optimisations include 
switching to a log grid, removing unneeded file IO (which was a big choke point
in the original code) and removing redundant calculations. In addition, nearly 
all the arrays have now been switched to STL::vector types, which should improve
memory management and eliminate bugs caused by lack of memory (or at least allow
the program to die gracefully).

You can either run the galaxy generator as a standalone program to generate a 
galaxy by running GenerateGalaxy, or you can run a MCMC chain to fit the 
parameters to a set of observations by running MCMC.

## Prerequisites

The GSL is required in the standard directories. The Makefile indicates G++, so
modify as needed.

## GenerateGalaxy:

The input for the code is several files that define the required parameters. The
file in.dbh has been reworked to look like this:
```
<halo_flag>                                                    
<halo cutoff> <halo velocity> <halo radius> <halo truncation width> <cusp>
    <halo stream fraction>
<disk_flag>                                                    
<bulge_flag>                                                   
<sersic n> <p> <bulge velocity> <bulge radius> <bulge stream fraction>
<blackhole_flag>                                               
<black hole mass>
<dr> <nr> <lmax> 
<do_file_io> <do_chisq_file_io> <nbody_flag> <chisq_flag>  
<n_spline> <n_iter>
```
if `<halo_flag>` is 0 there is no halo, otherwise we have a halo
if `<disk_flag>` 0 no disks, otherwise disk parameters are defined in in.diskpars
if `<bulge_pars>` 0 no bulge, otherwise we have a bulge
if `<blackhole_flag>` is 0 no black hole, otherwise... you guessed it
`<dr>` is the innermost bin - dr is a fraction of the smallest scale length above or
   in in.diskpars
`<nr>` is the number of radial bins 
`<lmax>` is the maximum harmonic
`<do_file_io> 0 if you don't want extraneous file IO done. If not 0, the code will 
   output dbh.dat, cordbh.dat, and all the other fun data files that the 
   original output.
`<do_chisq_file_io>` 0 if you don't want the chi square output files (with the model
   observations)to be generated
`<dochisq>` 0 if you don't want the chi square calculations done
`<nbody_flag>` to 1 if you want the N-body realisations to be generated.

in.gendenspsi:
```
<n_psi>
<n_int>
<max_iter> (for root bisection)
<tolerance_factor> (for root bisection)
<n_los>
<n_vel>
```
in.diskpars: 

This file has the following structure: 
```
<whitespace delineated parameters> 0
<whitespace delineated parameters> 0
<whitespace delineated parameters> 0
etc
```
where each line corresponds to a different disk component and the 0s indicate 
the end of the set of parameters for that component. The program counts the 
number of disks based on the number of 0s in this file. If the code encounters a
-1, it will stop reading disk parameters completely. This is useful for testing
say, 1 versus 2 disks; simply add a -1 to the end of the line for the last disk 
instead of deleting lines and then having to retype them. Currently, the code 
allows for an exponential disk; if you want extra components they need to be 
added in the DiskProfile functions of diskprofiles.cpp. This function takes an 
integer i as an argument, which corresponds to the i-th component of the disks 
(where i = 0,1... (# of disks-1)). For extra disk components, modify these 
functions to read 
```
if (i==0)
{
    <get f(r),f',f" here for first component>
}
else if (i==1)
{
    <get f(r),f',f" here for second component>
}
else if (i==2)
ETC.
```
Note that there are three DiskProfile functions, one for the density, and two 
for the density derivatives. This is faster than calculating all derivatives 
with each call of a single DiskProfile function, because not all three values 
are always needed when DiskProfile is called.

The current parameters for in.diskpars are the following:
```
<Exponential disk mass M_d> <disk scale radius R_d> <truncation radius> 
<sech scale height z_d> <truncation width> <Central dispersion sigma_0> 
<Dispersion scale length R_sigma> <inner truncation radius R_h> 
<inner truncation shape parameter alpha>
```
Other parameters will require modification of the GetDiskParameters function in 
diskprofiles.cpp to read in these variables. You'll want to declare any such 
additional parameters in galaxy.cpp and galaxy.h (in the galaxy struct) as 
global variables so that you can actually access them.

## MCMC:

For the chi square calculations you need an input file called 'mcmc_input'. This
file determines the starting parameters, what components you want, and anything 
else about the models for a MCMC run (this information is not read from the 
input files for GenerateGalaxy as described above). The structure of this file 
is as follows (where Halo, Bulge etc. are the actual strings "Halo", etc):
```
Halo <halo_flag>
<number> <name> <start_value> <sigma> <minimum> <maximum>
<etc for all halo params>

Bulge <bulge_flag>
<number> <name> <start_value> <sigma> <minimum> <maximum>
<etc for all bulge params>

BlackHole <smbh_flag>
<number> <name> <start_value> <sigma> <minimum> <maximum>
<etc for all smbh params>

Disk <disk_flag>
<number> <name> <start_value> <sigma> <minimum> <maximum>
<etc for all disk params for this component>

Disk <disk_flag>
<number> <name> <start_value> <sigma> <minimum> <maximum>
<etc for all disk params for this component>

<etc for all desired disk components>

Astro
<number> <name> <start_value> <sigma> <minimum> <maximum>
<etc for all astronomical params>

Error
<number> <name> <start_value> <sigma> <minimum> <maximum>
<etc for all error params>
```
It's your responsibility to ensure that the number of disk ML ratios in the 
astro parameters matches the number of disks. The code doesn't check this 
because doing so would have made it harder to add other astro parameters if 
required. The code does check to make sure you have as many error parameters as 
the number of data sets read in earlier, otherwise it spits out an error and 
exits.

The file 'mcmc.in' requires the following parameters:
```
<continue_flag> <chainlength> <scale_factor> <temperature> <scale_length>
<dr> <nr> <lmax> <dofileio> <n_psi> <n_int> <max_iter> <tolerance_factor>
<n_los> <n_vel>
```
where

`<continue_flag>` here is a flag telling the program whether to start from scratch
or read in the last line of parameters.out and start from there. If zero it's 
the former, otherwise it's the latter. The program will exit if mcmc_input and 
parameters.out do not match.

`<chainlength>` is the number of steps to take in the MCMC chain, and is added to 
whatever was already in parameters.out if <continue_flag> is true. 

`<scale_factor>` is the amount by which to multiply all the sigmas in mcmc_input, 
and is an easy way to scale them all by a constant if they're generally too 
large or too small. 

`<temperature>` is the starting temperature if simulated annealing is used. As of 
writing this has not been implemented.

`<scale_length>` is the scale length for the decline of the temperature if 
simulated annealing is used.

The remaining parameters are as described above.

The file parameters.out is the output for the MCMC program. Its default
structure is as follows
```
<step> <accepted steps> <chi square> <chisquare> <input parameters> <extra
parameters>
```
The parameters found in parameters.out include ALL the parameters found in
mcmc_input, even those that do not change. This is for two reasons: first, it
makes easier to write a program to take the parameters in any line and write
them to data files to get a specific model (otherwise we'd have to determine
which parameters were left out using mcmc_input, and if mcmc_input has changed
we probably can't do this at all). Second, if a run is interrupted and
restarted, it allows you to continue with previously fixed parameters now
variable and vice-versa without a change in the structure of the output file.

## Notes

If you want to change how MCMC handles the file IO (like whether or not it
writes the parameter files at each step), you'll want to edit the 
WriteDataFiles() and GenGalaxyAndGetChiSquare() functions in mcmc.cpp. These
functions write the input files for the next step of the chain and call the 
relevant galaxy generation functions respectively, but these would need to be 
modified to ignore file IO as well. This can be done by reading the appropriate 
element of the parameter array to each galaxy parameter in WriteDataFiles() and
by commenting out the sections involving file input in the functions in
GenGalaxyAndGetChiSquare().

## New features

This version of the code includes adiabatic contraction and an optional gas disk.
These features are not fully realised, so use at your own risk.
