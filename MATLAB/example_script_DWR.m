% This is an example of a script that calls postprocessing_DWR_ll.m
close all
clear, clc

% All parameter are in SI units
% Geometry parameters
geom = struct();
geom.medium = struct();
geom.medium.H = 3/1000;% [m]
geom.medium.R6 = 24.5/1000;% [m]
geom.medium.R5 = 23.5/1000;% [m]
geom.medium.Rr = 24;% [m]
geom.medium.R1 = 20/1000;% [m]
geom.medium.R3 = 28.7875/1000;% [m]
geom.medium.ringW = geom.medium.R6 - geom.medium.R5;% [m]
geom.medium.stepW = geom.medium.ringW;% [m]
geom.medium.R4 = geom.medium.R3 + geom.medium.stepW;% [m]
geom.medium.R2 = geom.medium.R1 - geom.medium.stepW;% [m]
geom.medium.G1 = geom.medium.R5 - geom.medium.R1;% [m]
geom.medium.G2 = geom.medium.R3 - geom.medium.R6;% [m]
geom.medium.inertia = (0.1 + 0.002)/1000;% system (rotor + DWR) inertia [Kg·m^2]
geom.medium.ICorrected = true;% True if torque data are inertia corrected

% Mesh parameters
mesh = struct();
mesh.ringSubs = 40;% Must be pair!!
if rem(mesh.ringSubs, 2) ~= 0
    warning('Ring subs must be pair!!')
end
mesh.upperBC = 'fb';% 'fb' for free interface (open contour) and 'ns' for no-slip (closed contour) boundary conditions
mesh.DOrder = 2;% 1 for first order and 2 for second order approximation in drags calculation

% Upper phase physical parameters
Bulk = struct();
Bulk.rho_bulk1 = 1000;% kg/m^3 water
Bulk.eta_bulk1 = 1e-03;% Ns/m water
Bulk.rho_bulk2 = 1.204;% kg/m^3 air
Bulk.eta_bulk2 = 1.813e-05;% Ns/m air

% Iterative scheme parameters
iteParams = struct();
iteParams.iteMax = 100;% maximum number of iterations
iteParams.tolMin = 1e-5;% threshold tolerance

% Input/output data
IO = struct();
IO.colIndexFreq = 1;% ordinal number of the data of the column that contains the modulus of the frequency
IO.colIndexAR = 2;% ordinal number of the data of the column that contains the modulus of the amplitude ratio
IO.colIndexDelta = 3;% ordinal number of the data of the column that contains the modulus of the phase shift
IO.colIndexAmp = 6;% ordinal number of the data of the column that contains the angular amplitude
IO.colIndexTorq = 4;% ordinal number of the data of the column that contains the Torque amplitude
IO.inputFilepath = pwd;% input filepath
IO.outputFilepath = pwd;% output filepath

% Execute postprocessingDWR_ll.m with the specified input data
[resultStruct, timeElapsedTotal] = postprocessingDWR_ll(geom.medium, mesh, Bulk, iteParams, IO);
save('MATLLAB_results', 'resultStruct', 'timeElapsedTotal')
