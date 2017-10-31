-- Input file for a 3X2V Simulation of NSTX
-- mpiexec -n 4 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i test136.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist
filename = "test136"

-- phase-space decomposition (5d)
phaseDecomp = DecompRegionCalc5D.CartProd { cuts = {1, 2, 2, 1, 1} }
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
nFrames = 4000
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

-- physical parameters
eV        = Lucee.ElementaryCharge
elcCharge = -eV
ionCharge = eV
eps0      = Lucee.Epsilon0 -- permittivity of free space
ionMass   = 2.014*Lucee.ProtonMass -- (Deuterium)
elcMass   = ionMass/400 -- REDUCED ELECTRON MASS Lucee.ElectronMass
elcTemp   = 40 -- [eV]
ionTemp   = 40 -- [eV]
B0 = 0.5  -- [T]
R0 = 0.85 -- [m], Major Radius
a0 = 0.5
n0 = 7*10^18 -- [1/m^3]
V_bias = 0 -- [V]
P_SOL = 5.4e6/10 -- [W]
xSource = -0.04 + R0 + a0 -- [m], source start coordinate
lambdaSource = 0.02 -- [m], characteristic length scale of density and temperature
bFieldReduction = 1 -- factor by which to reduce gradients in B
-- derived parameters
vtElc   = math.sqrt(elcTemp*eV/elcMass) -- only used for grid
vtIon   = math.sqrt(ionTemp*eV/ionMass)
omega_i = math.abs(ionCharge*B0*(R0/(R0+a0))/ionMass)
c_s = math.sqrt(elcTemp*eV/ionMass)
rho_s = c_s/omega_i
epsilonDensity = 1e-7 -- density floor to enforce

deltaX  = 50*rho_s
deltaY  = 200*rho_s
qR = 8 -- [m, from Myra 2016]
L_parallel = 4 -- [m]
-- grid parameters: number of cells
N_X = 18
N_Y = 72
N_Z = 10
N_VPARA = 10
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = (R0+a0) - deltaX/2
X_UPPER = (R0+a0) + deltaX/2
Y_LOWER = -deltaY/2
Y_UPPER = deltaY/2
Z_LOWER = -L_parallel/2
Z_UPPER = L_parallel/2
-- source parameters
L_s = 0.5*rho_s
r_s = 20*rho_s
sourceAmplitude = 1--n0*c_s/L_parallel
-- total volume
totalVol = (X_UPPER-X_LOWER)*(Y_UPPER-Y_LOWER)*(Z_UPPER-Z_LOWER)

VPARA_UPPER_ELC = 4*vtElc
VPARA_LOWER_ELC = -VPARA_UPPER_ELC
MU_LOWER_ELC = 0
MU_UPPER_ELC = 0.75*elcMass*(VPARA_UPPER_ELC*VPARA_UPPER_ELC)/(2*B0*(R0/(R0+a0)))

VPARA_UPPER_ION = 4*vtIon
VPARA_LOWER_ION = -VPARA_UPPER_ION
MU_LOWER_ION = 0
MU_UPPER_ION = 0.75*ionMass*(VPARA_UPPER_ION*VPARA_UPPER_ION)/(2*B0*(R0/(R0+a0)))

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
   periodicDirs = {1}, -- periodic in y only
   decomposition = confDecomp,
}

-- 2d spatial grid for sheath boundary conditions
grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
   periodicDirs = {1}, -- periodic in y only
   decomposition = sheathDecomp,
}

-- create 5d basis functions
basis_elc_5d = NodalFiniteElement5D.SerendipityElement {
   onGrid = grid_elc_5d,
   polyOrder = polyOrder,
   num1DGaussPoints = 3,
}
basis_ion_5d = NodalFiniteElement5D.SerendipityElement {
   onGrid = grid_ion_5d,
   polyOrder = polyOrder,
   num1DGaussPoints = 3,
}
-- create 3d basis functions
basis_3d = NodalFiniteElement3D.SerendipityElement {
   onGrid = grid_3d,
   polyOrder = polyOrder,
   num1DGaussPoints = 3,
}
-- create 2d basis functions for (x,y) grid
basis_2d = NodalFiniteElement2D.SerendipityElement {
   onGrid = grid_2d,
   polyOrder = polyOrder,
   num1DGaussPoints = 3,
}

-- distribution function for electrons
fElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
   --writeGhost = {1,1},
}
-- distribution function for ions
fIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
   --writeGhost = {1,1},
}
-- Electron Hamiltonians
hamilElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
hamilKeElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
hamilDerivKeElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
-- Ion Hamiltonians
hamilIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
hamilKeIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
hamilDerivKeIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
-- Magnetic field
function bFieldProfile(x)
  return B0*R0/(R0+a0+ (x - R0 - a0 )/bFieldReduction )
end
-- Updater to initialize KE part of hamil for electrons
initHamilKeElc = Updater.EvalOnNodes5D {
   onGrid = grid_elc_5d,
   basis = basis_elc_5d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,vPara,mu,t)
      return 0.5*elcMass*vPara*vPara + math.abs(mu)*bFieldProfile(x)
   end
}
runUpdater(initHamilKeElc, 0.0, 0.0, {}, {hamilKeElc})
-- Compute parallel velocity derivative of hamiltonian
calcHamilKeDerivElc = Updater.SOLDerivativeCalc5D {
  onGrid = grid_elc_5d,
  basis = basis_elc_5d,
  scaleFactor = 1/elcMass,
  dir = 3,
}
runUpdater(calcHamilKeDerivElc, 0.0, 0.0, {hamilKeElc}, {hamilDerivKeElc})

-- Updater to initialize KE part of hamil for ions
initHamilKeIon = Updater.EvalOnNodes5D {
   onGrid = grid_ion_5d,
   basis = basis_ion_5d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,vPara,mu,t)
      return 0.5*ionMass*vPara*vPara + math.abs(mu)*bFieldProfile(x)
   end
}
runUpdater(initHamilKeIon, 0.0, 0.0, {}, {hamilKeIon})
-- Compute parallel velocity derivative of hamiltonian
calcHamilKeDerivIon = Updater.SOLDerivativeCalc5D {
  onGrid = grid_ion_5d,
  basis = basis_ion_5d,
  scaleFactor = 1/ionMass,
  dir = 3,
}
runUpdater(calcHamilKeDerivIon, 0.0, 0.0, {hamilKeIon}, {hamilDerivKeIon})

-- Magnetic Field (3D)
bField3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bField3d
initb = Updater.EvalOnNodes3D {
  onGrid = grid_3d,
  basis = basis_3d,
  shareCommonNodes = false,
  evaluate = function(x,y,z,t)
    return bFieldProfile(x)
  end
}
runUpdater(initb, 0.0, 0.0, {}, {bField3d})
bField3d:sync()

-- heat flux at edge diagnostic
heatFluxAtEdgeElcCalc = Updater.SOLHeatFluxCalc {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/elcMass, -- B will be an input
}
heatFluxLowerElc = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
heatFluxUpperElc = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}

heatFluxAtEdgeIonCalc = Updater.SOLHeatFluxCalc {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/ionMass, -- B will be an input
}
heatFluxLowerIon = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
heatFluxUpperIon = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}

fElc:read("fElc_663.h5")
fElc:sync()
fIon:read("fIon_663.h5")
fIon:sync()

  -- compute heat flux at edges
runUpdater(heatFluxAtEdgeElcCalc, 0, 0, {fElc,hamilDerivKeElc,hamilKeElc,bField3d},
  {heatFluxLowerElc, heatFluxUpperElc})
runUpdater(heatFluxAtEdgeIonCalc, 0, 0, {fIon,hamilDerivKeIon,hamilKeIon,bField3d},
  {heatFluxLowerIon, heatFluxUpperIon})

-- write out heat fluxes
heatFluxUpperElc:write( string.format("heatFluxUpperElc_%d.h5", 0), 0)
heatFluxUpperIon:write( string.format("heatFluxUpperIon_%d.h5", 0), 0)
heatFluxLowerElc:write( string.format("heatFluxLowerElc_%d.h5", 0), 0)
heatFluxLowerIon:write( string.format("heatFluxLowerIon_%d.h5", 0), 0)
