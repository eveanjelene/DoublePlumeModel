{\rtf1\ansi\ansicpg1252\cocoartf2513
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\csgray\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww14920\viewh12340\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 DoublePM description\
\
07/04/2020\
\
\
1) run_doublePM.m\
\
MATLAB script to run the double plume integral model. Within this script, set parameters to decide whether the mass from the inner plume is added to the single plume at the top of the double plume region or not, and which grain size distribution. The initial velocity, temperature,  power law coefficient and gas content for the inner plume is defined as the input matrix. double_PM function is called within. \
\
2) double_PM.m\
\
MATLAB function of the double plume model. \
Inputs: initial mass flow rate, eruption source parameters of inner plume, eruption source parameters of outer plume, mode of inner plume grain-size distribution, addMass\
Output: output struct that contains the dynamic profiles and the ode output (mass, momentum, energy) profiles. The structure of the data is that they are vectors containing this data for the inner plume, outer plume and the single plume region. See schematic in model description for further detail of what/where these regions are in a plume.\
\
This function define the constants, atmospheric constants and the initial source parameters for the plume model. The ode solvers are called  and are the following.\
a i) doublePlume_return2one_ODE_ED.m- ode45 solver function of the double plume region of the model\
a ii) doublePlume_return2one_stopinner.m  - ode45 stop function. Stops the ode45 when the inner plume velocity reaches 0.01 m/s.\
\
b i) onePume_ode.m - ode45 solver function of the single plume region of the model\
b ii) onePlume_stop.m  - ode45 stop function. Stops the ode45 when the outer plume velocity reaches 0.01 m/s.\
\
The output of the ode solver functions are used to define the dynamic profiles of the inner, outer and single plume parts of the model. These, along with he ode solver function output, are assigned to the output struct to output from the model.\
\
If the double plume model doesn\'92t stop because of running out of mass or velocity not reaching 0 (i.e. stopped because a particular size is no longer supported by the plume), model reruns until it does. Each output is stored in a cell. \
\
Within doublePlume_return2one_ODE_ED.m is sedInner.m, which calculates the amount of fallout of each grain-size from the inner plume at each step within the ODE.\
\
c) sedInner.m\
MATLAB function which calculates the amount of fallout of each grain-size from the inner plume at each step within the ODE.\
Inputs: velocity of the inner plume at last height step, radius of the inner plume at last height, vector of the grain-sizes in m, gravitational acceleration constant, density of the inner plume at last height, viscosity, structure of constants, vector of MFR of each grain-size at last height step\
Outputs: vector of MFR of each size the falls out of the inner plume\
\
d) addFallout.m\
MATLAB function that\cf2 \cb3  determines the what grain-sizes of the inner is supported by the surrounding tephra plume\
In\cf0 \cb1 put: velocity of inner plume from last height, vector of grain-sizes in m, gravitational acceleration constant, bulk density of inner plume from last height, viscosity, struct containing constants, MFR distribution of inner plume from last height\cf2 \cb3 \
\cf0 \cb1 Output : vector of MFR related to corresponding grain-size that can be supported by the \cf2 \cb3 surrounding tephra plume\cf0 \cb1 \
\
e) add2Surrounding.m\
\
MATLAB function to determine a how much of the sedimenting material from the inner plume goes into the surrounding tephra plume\
Inputs:  output of last  run of double plume model, struct containing constants\
Output: new source parameters for the next run of the double plume\
\
Output from re-runs are put together to form one vector for each output for the whole double plume region in the following function\'85\
\
f) putTogether.m - loop through ode outputs stored in the cell storage. Add each variable to a vector. Vector is then reduced to be just the size of the compiled ode outputs. Outputs the coupled vectors of the ode outputs from all of the reruns.\
\
e) addFallout.m\
\
MATLAB function that\cf2 \cb3  determines the what grain-sizes of the inner is supported by the start of the single plume\
In\cf0 \cb1 put: velocity of inner plume from last height, vector of grain-sizes in m, gravitational acceleration constant, bulk density of inner plume from last height, viscosity, struct containing constants, MFR distribution of inner plume from last height\cf2 \cb3 \
\cf0 \cb1 Output : vector of MFR related to corresponding grain-size that can be supported by the \cf2 \cb3 surrounding tephra plume\
\cf0 \cb1 \
3) gaussDistribution.m\
\
MATLAB function to determine a gaussian grain-size distribution \
Inputs:  mode grain-size, standard deviation, grain-sizes\
Output: mass fraction\
\
Called within double_PM.m to define the initial grain size distribution of the inner and outer plume regions.}