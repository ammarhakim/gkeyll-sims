-- test of RunningAverageOfFieldCalc
-- mpiexec -n 4 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i test66.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist

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
tEnd = 1000e-6
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 1000
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

-- physical parameters
eV        = Lucee.ElementaryCharge
elcCharge = -eV
ionCharge = eV
eps0      = Lucee.Epsilon0 -- permittivity of free space
elcMass   = Lucee.ElectronMass
ionMass   = 2.014*Lucee.ProtonMass 
elcTemp   = 5 -- [eV]
ionTemp   = 5 -- [eV]
B0 = 0.076  -- [T]
n0 = 10^16 -- [1/m^3]
N = 4 -- field line turns
kPerpTimesRhoS = 0.2
sourceTemp     = 2 -- [eV] -- temperature to base grid off of
sourceTempInit = 2 -- [eV] -- temperature to init source with
-- derived parameters
vtElc   = math.sqrt(sourceTemp*eV/elcMass)
omega_e = math.abs(elcCharge*B0/elcMass)
rho_e   = vtElc/omega_e

vtIon   = math.sqrt(sourceTemp*eV/ionMass)
omega_i = math.abs(ionCharge*B0/ionMass)
rho_i   = vtIon/omega_i

rho_s = math.sqrt(elcTemp*eV/ionMass)/omega_i
deltaX  = 64*rho_s
deltaY  = 64*rho_s/N
R          = 200*rho_s -- [m]
L_parallel = 2*math.pi*R*N -- [m]
-- grid parameters: number of cells
N_X = 16
N_Y = 4
N_Z = 4
N_VPARA = 6
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = R
X_UPPER = R + deltaX
Y_LOWER = -deltaY/2
Y_UPPER = deltaY/2
Z_LOWER = 0
Z_UPPER = L_parallel
-- source parameters
sourceLambda_UH    = 5*rho_s
sourceX_UH         = X_LOWER + 15*rho_s
sourceAmplitude_UH = 1.5e20 -- [1/(m^3*sec)] 
sourceLambda_EC    = 2.5*rho_s
sourceX_EC         = X_LOWER + 35*rho_s
sourceAmplitude_EC = 1e20 -- [1/(m^3*sec)] 
-- total volume
totalVol = (X_UPPER-X_LOWER)*(Y_UPPER-Y_LOWER)*(Z_UPPER-Z_LOWER)

VPARA_UPPER_ELC = 4*vtElc
VPARA_LOWER_ELC = -VPARA_UPPER_ELC
MU_LOWER_ELC = 0
MU_UPPER_ELC = elcMass*(VPARA_UPPER_ELC*VPARA_UPPER_ELC)/(2*B0)

VPARA_UPPER_ION = 4*vtIon
VPARA_LOWER_ION = -VPARA_UPPER_ION
MU_LOWER_ION = 0
MU_UPPER_ION = ionMass*(VPARA_UPPER_ION*VPARA_UPPER_ION)/(2*B0)

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
   periodicDirs = {1}, -- periodic in y only
   decomposition = phaseDecomp,
}
grid_ion_5d = Grid.RectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER_ION, MU_LOWER_ION},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER_ION, MU_UPPER_ION},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
   periodicDirs = {1}, -- periodic in y only
   decomposition = phaseDecomp,
}

-- 3d spatial grid
grid_3d = Grid.RectCart3D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER},
   cells = {N_X, N_Y, N_Z},
   --periodicDirs = {0, 1, 2},
   periodicDirs = {1},
   decomposition = confDecomp,
}

-- 2d spatial grid for sheath boundary conditions
grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
   periodicDirs = {1},
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
-- ion guiding center number density
numDensityIon = DataStruct.Field3D {
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
-- ion background density
bgDensityIon = DataStruct.Field3D {
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
-- ion fluctuating density
fluctDensityIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

-- electron background temperature
bgTempElc = DataStruct.Field3D {
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
-- electron temperature
tempParaElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron background temperature
bgTempParaElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron fluctuating temperature
fluctTempParaElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- ion temperature
tempParaIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- ion background temperature
bgTempParaIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- ion fluctuating temperature
fluctTempParaIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- parallel velocity moment
mom1dir3Ion = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
bgMom1dir3Ion = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
fluctMom1dir3Ion = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- parallel velocity moment
mom1dir3Elc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
bgMom1dir3Elc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
fluctMom1dir3Elc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

-- to divide two 3d fields together on nodes
divide3dCalc = Updater.FieldArithmeticUpdater3D {
  onGrid = grid_3d,
  basis = basis_3d,
  evaluate = function(fld1, fld2)
    return fld1/fld2
  end
}
divisionResult = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

-- start computing bg info at frame 100
tCurr = tCurr + 100*tFrame

for frame = 100, nFrames do
  Lucee.logInfo (string.format("-- Computing average at frame %g", frame))
  Lucee.logInfo ("")

  bgPhi:clear(0.0)
  bgTempElc:clear(0.0)
  bgPressureElc:clear(0.0)
  bgDensityElc:clear(0.0)
  bgDensityIon:clear(0.0)
  bgTempParaIon:clear(0.0)
  bgTempParaElc:clear(0.0)
  bgMom1dir3Ion:clear(0.0)
  bgMom1dir3Elc:clear(0.0)
  -- first compute average of potential, temperature, pressure, and density fields
  -- using the past 101 frames
  for bgFrame = frame-100, frame do
    phi:read("phi_" .. bgFrame .. ".h5")
    temperatureElc:read("tElc_" .. bgFrame .. ".h5")
    pressureElc:read("pressureElc_" .. bgFrame .. ".h5")
    numDensityElc:read("nElc_" .. bgFrame .. ".h5")
    numDensityIon:read("nIon_" .. bgFrame .. ".h5")
    tempParaElc:read("tParaElc_" .. bgFrame .. ".h5")
    tempParaIon:read("tParaIon_" .. bgFrame .. ".h5")
    mom1dir3Elc:read("mom1dir3Elc_" .. bgFrame ..".h5")
    mom1dir3Ion:read("mom1dir3Ion_" .. bgFrame ..".h5")

    -- divide mom1dir3Ion by ion density
    runUpdater(divide3dCalc, 0, 0, {mom1dir3Ion, numDensityIon}, {divisionResult})

    -- accccumulate to bg structures
    bgPhi:accumulate(1/(100+1), phi)
    bgTempElc:accumulate(1/(100+1), temperatureElc)
    bgPressureElc:accumulate(1/(100+1), pressureElc)
    bgDensityElc:accumulate(1/(100+1), numDensityElc)
    bgDensityIon:accumulate(1/(100+1), numDensityIon)
    bgTempParaElc:accumulate(1/(100+1), tempParaElc)
    bgTempParaIon:accumulate(1/(100+1), tempParaIon)
    bgMom1dir3Elc:accumulate(1/(100+1), mom1dir3Elc)
    bgMom1dir3Ion:accumulate(1/(100+1), divisionResult)
  end

  -- write out background fields used at this time
  bgPhi:write( string.format("bgPhi_%d.h5", frame), tCurr)
  bgTempElc:write( string.format("bgTempElc_%d.h5", frame), tCurr)
  bgDensityElc:write( string.format("bgDensityElc_%d.h5", frame), tCurr)
  bgDensityIon:write( string.format("bgDensityIon_%d.h5", frame), tCurr)
  bgPressureElc:write( string.format("bgPressureElc_%d.h5", frame), tCurr)
  bgTempParaElc:write( string.format("bgTempParaElc_%d.h5", frame), tCurr)
  bgTempParaIon:write( string.format("bgTempParaIon_%d.h5", frame), tCurr)
  bgMom1dir3Elc:write( string.format("bgMom1dir3Elc_%d.h5", frame), tCurr)
  bgMom1dir3Ion:write( string.format("bgMom1dir3Ion_%d.h5", frame), tCurr)

  -- read in total fields at this time
  phi:read("phi_" .. frame .. ".h5")
  temperatureElc:read("tElc_" .. frame .. ".h5")
  pressureElc:read("pressureElc_" .. frame .. ".h5")
  numDensityElc:read("nElc_" .. frame .. ".h5")
  numDensityIon:read("nIon_" .. frame .. ".h5")
  tempParaElc:read("tParaElc_" .. frame .. ".h5")
  tempParaIon:read("tParaIon_" .. frame .. ".h5")
  mom1dir3Elc:read("mom1dir3Elc_" .. frame .. ".h5")
  mom1dir3Ion:read("mom1dir3Ion_" .. frame .. ".h5")
  -- divide mom1dir3Ion by ion density
  runUpdater(divide3dCalc, 0, 0, {mom1dir3Ion, numDensityIon}, {divisionResult})

  -- compute fluctuating fields at this time (total - bg)
  fluctPhi:combine(1,phi,-1,bgPhi)
  fluctTempElc:combine(1,temperatureElc,-1,bgTempElc)
  fluctDensityElc:combine(1,numDensityElc,-1,bgDensityElc)
  fluctDensityIon:combine(1,numDensityIon,-1,bgDensityIon)
  fluctPressureElc:combine(1,pressureElc,-1,bgPressureElc)
  fluctTempParaElc:combine(1,tempParaElc,-1,bgTempParaElc)
  fluctTempParaIon:combine(1,tempParaIon,-1,bgTempParaIon)
  fluctMom1dir3Elc:combine(1,mom1dir3Elc,-1,bgMom1dir3Elc)
  fluctMom1dir3Ion:combine(1,divisionResult,-1,bgMom1dir3Ion)

  -- run updater for average
  --runUpdater(bgDensCalc, 0, 0, {numDensityElc}, {bgDensityElc, varDensityElc})
  -- write out fluctuating data data
  fluctPhi:write( string.format("fluctPhi_%d.h5", frame), tCurr)
  fluctTempElc:write( string.format("fluctTempElc_%d.h5", frame), tCurr)
  fluctDensityElc:write( string.format("fluctDensityElc_%d.h5", frame), tCurr)
  fluctDensityIon:write( string.format("fluctDensityIon_%d.h5", frame), tCurr)
  fluctPressureElc:write( string.format("fluctPressureElc_%d.h5", frame), tCurr)
  fluctTempParaElc:write( string.format("fluctTempParaElc_%d.h5", frame), tCurr)
  fluctTempParaIon:write( string.format("fluctTempParaIon_%d.h5", frame), tCurr)
  fluctMom1dir3Elc:write( string.format("fluctMom1dir3Elc_%d.h5", frame), tCurr)
  fluctMom1dir3Ion:write( string.format("fluctMom1dir3Ion_%d.h5", frame), tCurr)

  -- advance time
  tCurr = tCurr + tFrame
end
