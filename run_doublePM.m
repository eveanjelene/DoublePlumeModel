clear
close all

set(groot,'defaultLineLineWidth',1.5)

%% Set inputs to doublePM
 
    addMass = 2;                                                           % Variable to define if to add remaining MFR of inner, once it has stop, to start of the single plume 0 = don't add, 1 = add all, 2 = add if supported based on velocities
    ratio   = 0.7;                                                         % Partition ration of the initial mass flow rate (MFR) 
    gs_size = -5;                                                          % Median grain-size for the gaussian distribution
    
    MER = 1e5;                                                             % Initial MFR
    
    %create empty vecrtors for the inputs for inner and out plum
    %format = velocity, temperature, gass mass fraction
    input_i = [150 1500 0.01];  
    input_o = [150 1200 0.03];
    
%% run double PM 
        
    [output] = double_PM(MER,input_i,input_o,addMass,gs_size,ratio);  

