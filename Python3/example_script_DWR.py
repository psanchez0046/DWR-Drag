# %% This is an example of a script that calls postprocessing_BiconeDWR_ll
import os
import DWRfuns
from ypstruct import struct

# Geometry parameters
geom = struct()
geom.H = 3/1000# distance between interface and cup bottom [m]
geom.R6 = 24.5/1000
geom.R5 = 23.5/1000
geom.Rr = 24
geom.R1 = 20/1000
geom.R3 = 28.7875/1000
geom.ringW = geom.R6 - geom.R5
geom.stepW = geom.ringW
geom.R4 = geom.R3 + geom.stepW
geom.R2 = geom.R1 - geom.stepW
geom.G1 = geom.R5 - geom.R1 
geom.G2 = geom.R3 - geom.R6
geom.inertia = (0.1 + 0.002)/1000 # system (rotor + DWR) inertia [Kg·m^2]
geom.ICorrected = True

# Mesh parameters
mesh = struct()
mesh.ringSubs = 40# Must be pair!!
mesh.upperBC = 'fb'
mesh.DOrder = 2

# Bulk phases physical parameters
Bulk = struct()
Bulk.rho_bulk1 = 1000# kg/m^3 water
Bulk.eta_bulk1 = 1e-3# Ns/m water
Bulk.rho_bulk2 = 1.204# kg/m^3 air
Bulk.eta_bulk2 = 1.813e-5# Ns/m air

# Iterative scheme parameters
iteParams = struct()
iteParams.iteMax = 100# maximum number of iterations
iteParams.tolMin = 1e-5# threshold tolerance

# Input/output data
IO = struct()
IO.colIndexFreq = 0# ordinal number of the data of the column that contains the modulus of the frequency
IO.colIndexAR = 1# ordinal number of the data of the column that contains the modulus of the amplitude ratio
IO.colIndexDelta = 2# ordinal number of the data of the column that contains the modulus of the phase shift
IO.colIndexAmp = 4# ordinal number of the data of the column that contains the angular amplitude
IO.colIndexTorq = 3# ordinal number of the data of the column that contains the Torque amplitude
IO.inputFilepath = os.getcwd()# input filepath
IO.outputFilepath = os.getcwd()# output filepath

## Execute postprocessingBiconeCPC.m with the specified input data
resultStruct = DWRfuns.postprocessingDWR_ll(geom, mesh, Bulk, iteParams, IO)
