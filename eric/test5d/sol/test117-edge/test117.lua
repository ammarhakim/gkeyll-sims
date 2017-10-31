-- lua file to read in total fields and compute background and fluctuating fields
-- mpiexec -n 4 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i test117.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist -r

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
tEnd = 4157e-6
nFrames = 4157 
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

-- outward electron flux
edgeFluxElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- outward ion flux
edgeFluxIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- field of ones
oneField = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
oneField:clear(1)

productAtEdgeCalc = Updater.SOLFluxAcrossEdgeCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
}

elcFluxVec = DataStruct.DynVector { numComponents = 2, }
ionFluxVec = DataStruct.DynVector { numComponents = 2, }

-- Calculate elc and ion source rates
elcSource = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
ionSource = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
fieldInt3d = Updater.IntegrateGeneralField3D {
  onGrid = grid_3d,
  basis = basis_3d,
}
tempVec = DataStruct.DynVector { numComponents = 1, }
elcSource:read("elcSource_0.h5")
runUpdater(fieldInt3d, 0, 0, {elcSource}, {tempVec})
Lucee.logInfo( string.format("Elc source rate = %g",tempVec:lastInsertedData()) )
ionSource:read("ionSource_0.h5")
runUpdater(fieldInt3d, 0, 0, {ionSource}, {tempVec})
Lucee.logInfo( string.format("Ion source rate = %g",tempVec:lastInsertedData()) )

--for frame = Lucee.RestartFrame, nFrames do
--  Lucee.logInfo (string.format("-- Reading data at frame %g", frame))
--  Lucee.logInfo ("")

--  edgeFluxElc:read("edgeFluxElc_" .. frame .. ".h5")
--  edgeFluxIon:read("edgeFluxIon_" .. frame .. ".h5")

-- tCurr = frame*tFrame
--  runUpdater(productAtEdgeCalc, tCurr, 0, {edgeFluxElc,oneField}, {elcFluxVec})
--  runUpdater(productAtEdgeCalc, tCurr, 0, {edgeFluxIon,oneField}, {ionFluxVec})

--end

--elcFluxVec:write(string.format("elcFluxVec_%d.h5", 1))
--ionFluxVec:write(string.format("ionFluxVec_%d.h5", 1))
