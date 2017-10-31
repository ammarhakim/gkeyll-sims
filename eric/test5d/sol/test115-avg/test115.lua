-- lua file to read in total fields and compute background and fluctuating fields
-- mpiexec -n 4 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i test115.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist

-- phase-space decomposition (5d)
phaseDecomp = DecompRegionCalc5D.CartProd { cuts = {4, 1, 1, 1, 1} }
-- configuration space decomposition (3d)
confDecomp = DecompRegionCalc3D.SubCartProd5D {
   decomposition = phaseDecomp,
   collectDirections = {0, 1, 2},
}
-- logical sheath decomposition for solve on 2d planes
sheathDecomp = DecompRegionCalc2D.SubCartProd3D {
   decomposition = confDecomp,
   collectDirections = {0, 1},
}

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.1
-- parameters to control time-stepping
tStart = 0.0
tEnd = 3158e-6
nFrames = 3158 
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

-- physical parameters
eV        = Lucee.ElementaryCharge
elcCharge = -eV
ionCharge = eV
eps0      = Lucee.Epsilon0 -- permittivity of free space
ionMass   = 3.973*Lucee.ProtonMass 
elcMass   = ionMass/400 -- REDUCED ELECTRON MASS Lucee.ElectronMass
elcTemp   = 6 -- [eV]
ionTemp   = 1 -- [eV]
B0 = 0.0398  -- [T]
n0 = 2*10^18 -- [1/m^3]
-- derived parameters
vtElc   = math.sqrt(0.5*elcTemp*eV/elcMass) -- only used for grid
vtIon   = math.sqrt(ionTemp*eV/ionMass)
omega_i = math.abs(ionCharge*B0/ionMass)

c_s = math.sqrt(elcTemp*eV/ionMass)
rho_s = c_s/omega_i

deltaX  = 100*rho_s
deltaY  = 100*rho_s
R          = 40*rho_s -- [m]
L_parallel = 36*R -- [m]
-- grid parameters: number of cells
N_X = 36
N_Y = 36
N_Z = 10
N_VPARA = 10
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = -deltaX/2
X_UPPER = deltaX/2
Y_LOWER = -deltaY/2
Y_UPPER = deltaY/2
Z_LOWER = -L_parallel/2
Z_UPPER = L_parallel/2
-- source parameters
L_s = 0.5*rho_s
r_s = 20*rho_s
sourceAmplitude = 0.03*n0*c_s/R
-- total volume
totalVol = (X_UPPER-X_LOWER)*(Y_UPPER-Y_LOWER)*(Z_UPPER-Z_LOWER)

VPARA_UPPER_ELC = 4*vtElc
VPARA_LOWER_ELC = -VPARA_UPPER_ELC
MU_LOWER_ELC = 0
MU_UPPER_ELC = 0.75*elcMass*(VPARA_UPPER_ELC*VPARA_UPPER_ELC)/(2*B0)

VPARA_UPPER_ION = 4*vtIon
VPARA_LOWER_ION = -VPARA_UPPER_ION
MU_LOWER_ION = 0
MU_UPPER_ION = 0.75*ionMass*(VPARA_UPPER_ION*VPARA_UPPER_ION)/(2*B0)

-- A generic function to run an updater.
function runUpdater(updater, currTime, timeStep, inpFlds, outFlds)
   updater:setCurrTime(currTime)
   if inpFlds then
      updater:setIn(inpFlds)
   end
   if outFlds then
      updater:setOut(outFlds)
   end
   return updater:advance(currTime+timeStep)
end

-- full 5d phase space grids
grid_elc_5d = Grid.RectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER_ELC, MU_LOWER_ELC},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER_ELC, MU_UPPER_ELC},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
   decomposition = phaseDecomp,
}
grid_ion_5d = Grid.RectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER_ION, MU_LOWER_ION},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER_ION, MU_UPPER_ION},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
   decomposition = phaseDecomp,
}

-- 3d spatial grid
grid_3d = Grid.RectCart3D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER},
   cells = {N_X, N_Y, N_Z},
   decomposition = confDecomp,
}

-- 2d spatial grid for sheath boundary conditions
grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
   decomposition = sheathDecomp,
}

-- create 5d basis functions
basis_elc_5d = NodalFiniteElement5D.SerendipityElement {
   onGrid = grid_elc_5d,
   polyOrder = polyOrder,
}
basis_ion_5d = NodalFiniteElement5D.SerendipityElement {
   onGrid = grid_ion_5d,
   polyOrder = polyOrder,
}
-- create 3d basis functions
basis_3d = NodalFiniteElement3D.SerendipityElement {
   onGrid = grid_3d,
   polyOrder = polyOrder,
}

-- stores result of potential solve
phi = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- stores result of background potential diagnostic
bgPhi = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- variance
varPhi = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- potential fluctuations
fluctPhi = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron guiding center number density
numDensityElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron total tempearture
temperatureElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron total pressure
pressureElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron background pressure
bgPressureElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- variance
varPressureElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron fluctuating pressure
fluctPressureElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

bgPressureCalc = Updater.RunningAverageOfFieldCalc3D {
  onGrid = grid_3d,
  sampleSize = 100,
}

bgPhiCalc = Updater.RunningAverageOfFieldCalc3D {
  onGrid = grid_3d,
  sampleSize = 100,
}

bgDensCalc = Updater.RunningAverageOfFieldCalc3D {
  onGrid = grid_3d,
  sampleSize = 100,
}
-- electron background density
bgDensityElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- variance
varDensityElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron fluctuating density
fluctDensityElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

bgTempCalc = Updater.RunningAverageOfFieldCalc3D {
  onGrid = grid_3d,
  sampleSize = 100,
}
-- electron background temperature
bgTempElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- variance
varTempElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron fluctuating temperature
fluctTempElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- total current = sum(<v>_s*q_s)
totalCurrent = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- <v>
mom1dir3Elc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
mom1dir3Ion = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}


-- start computing bg info at frame 100
tCurr = tCurr + 2000*tFrame

for bgFrame = 2000, nFrames do
  numDensityElc:read("nElc_" .. bgFrame .. ".h5")
  bgDensityElc:accumulate(1/(nFrames-2000+1), numDensityElc)
end

for frame = 2000, nFrames do
  Lucee.logInfo (string.format("-- Computing average at frame %g", frame))
  Lucee.logInfo ("")

  --bgPhi:clear(0.0)
  --bgTempElc:clear(0.0)
  --bgPressureElc:clear(0.0)
  --bgDensityElc:clear(0.0)
  -- first compute average of potential, temperature, pressure, and density fields
  -- using the past 101 frames
  --for bgFrame = frame-50, frame do
    --phi:read("phi_" .. bgFrame .. ".h5")
    --temperatureElc:read("tElc_" .. bgFrame .. ".h5")
    --pressureElc:read("pressureElc_" .. bgFrame .. ".h5")
  --  numDensityElc:read("nElc_" .. bgFrame .. ".h5")
    --mom1dir3Elc:read("mom1dir3Elc_" .. bgFrame .. ".h5")
    --mom1dir3Ion:read("mom1dir3Ion_" .. bgFrame .. ".h5")

    -- accccumulate to bg structures
    --bgPhi:accumulate(1/(100+1), phi)
    --bgTempElc:accumulate(1/(100+1), temperatureElc)
    --bgPressureElc:accumulate(1/(100+1), pressureElc)
  --  bgDensityElc:accumulate(1/(50+1), numDensityElc)
  --end

  -- write out background fields used at this time
  --bgPhi:write( string.format("bgPhi_%d.h5", frame), tCurr)
  --bgTempElc:write( string.format("bgTempElc_%d.h5", frame), tCurr)
  bgDensityElc:write( string.format("bgDensityElc_%d.h5", frame), tCurr)
  --bgPressureElc:write( string.format("bgPressureElc_%d.h5", frame), tCurr)

  -- read in total fields at this time
  --phi:read("phi_" .. frame .. ".h5")
  --temperatureElc:read("tElc_" .. frame .. ".h5")
  --pressureElc:read("pressureElc_" .. frame .. ".h5")
  numDensityElc:read("nElc_" .. frame .. ".h5")

  -- compute fluctuating fields at this time (total - bg)
  --fluctPhi:combine(1,phi,-1,bgPhi)
  --fluctTempElc:combine(1,temperatureElc,-1,bgTempElc)
  fluctDensityElc:combine(1,numDensityElc,-1,bgDensityElc)
  --fluctPressureElc:combine(1,pressureElc,-1,bgPressureElc)

  -- run updater for average
  --runUpdater(bgDensCalc, 0, 0, {numDensityElc}, {bgDensityElc, varDensityElc})
  -- write out fluctuating data data
  --fluctPhi:write( string.format("fluctPhi_%d.h5", frame), tCurr)
  --fluctTempElc:write( string.format("fluctTempElc_%d.h5", frame), tCurr)
  fluctDensityElc:write( string.format("fluctDensityElc_%d.h5", frame), tCurr)
  --fluctPressureElc:write( string.format("fluctPressureElc_%d.h5", frame), tCurr)

  -- advance time
  tCurr = tCurr + tFrame
end
