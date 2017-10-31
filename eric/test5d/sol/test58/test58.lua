-- Input file for a 3X2V SOL Simulation
-- 01-28-2016 Creation
-- 03-15-2016 First attempt at input file with kinetic elections and ions
-- Command line parameters
-- mpiexec -n 1 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i test58.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist
-- added positivity back in
-- number of parallel velocity space cells changed to 6, max v = 6*vt
-- changed boundary condition on lhs to be the same as that on rhs
-- 4-10-16 testing fixes to sheath boundary condition
-- 4-10-16 testing improvement of using a precomputed dH/dv field
-- 4-14-16: Load with averaged profiles and restart simulation
-- 4-15-16: Added code to look at field in a cell
-- 4-16-16: Nothing changed, but running with a modified potential solve
-- 4-18-16: studying flux imbalance in test23
-- 4-20-16: running with phi = 30 eV on the side walls
-- 4-20-16: changed source init from eval to projection. do not restart.
-- 4-22-16: a test for the positivity updater with energy preservation
-- 4-27-16: testing positivity updater in standard SMT simulation
-- 5-1-16: running with walls set to 15 ev
-- 5-3-16: modified to have biased sheath boundary conditions
-- 5-5-16: testing new density initialization. reduced vMax to 4*vt
-- 5-14-16: testing <u> calculation
-- 5-18-16: added LB operator for e-e collisions

-- phase-space decomposition (5d)
phaseDecomp = DecompRegionCalc5D.CartProd { cuts = {1, 1, 1, 1, 1} }
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
tEnd = 10e-6
dtSuggested = 1e-8 --0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 5
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart
-- relative size of finite difference for Newton's method jacobian
epsilonU = 1
epsilonT = 1e-8
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
N_X = 2
N_Y = 2
N_Z = 2
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
-- create 2d basis functions
basis_2d = NodalFiniteElement2D.SerendipityElement {
   onGrid = grid_2d,
   polyOrder = polyOrder,
}

-- distribution function for electrons
fElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
   --writeGhost = {1,1},
}
-- distribution function after positivity drag
fElcPositivityDrag = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
   --writeGhost = {1,1},
}
-- delta positivity term
fElcPositivityDelta = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
   --writeGhost = {1,1},
}
-- for reverting distribution function when time step fails
fDupElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
-- for rk3 stages
f1Elc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
fNewElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
-- distribution function for ions
fIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
   --writeGhost = {1,1},
}
-- distribution function after positivity drag
fIonPositivityDrag = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
   --writeGhost = {1,1},
}
-- delta positivity term
fIonPositivityDelta = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
   --writeGhost = {1,1},
}
-- for reverting distribution function when time step fails
fDupIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
-- for rk3 stages
f1Ion = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
fNewIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
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
hamilParaElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
hamilDerivKeElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
hamilDupElc = DataStruct.Field5D {
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
hamilParaIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
hamilDerivKeIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
hamilDupIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
-- Magnetic field
function bFieldProfile(x)
  return B0/x
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
-- compute parallel energy
initHamilParaElc = Updater.EvalOnNodes5D {
   onGrid = grid_elc_5d,
   basis = basis_elc_5d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,vPara,mu,t)
      return 0.5*elcMass*vPara*vPara
   end
}
runUpdater(initHamilParaElc, 0.0, 0.0, {}, {hamilParaElc})

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
-- compute parallel energy
initHamilParaIon = Updater.EvalOnNodes5D {
   onGrid = grid_ion_5d,
   basis = basis_ion_5d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,vPara,mu,t)
      return 0.5*ionMass*vPara*vPara
   end
}
runUpdater(initHamilParaIon, 0.0, 0.0, {}, {hamilParaIon})


function elcTempProfile(x)
  return elcTemp
end

function fElcProfile(x,y,z,v,mu)
  return (2*math.pi*elcTempProfile(x)*eV/elcMass)^(-3/2)*
    math.exp(-elcMass*v^2/(2*elcTempProfile(x)*eV))*
    math.exp(-math.abs(mu)*bFieldProfile(x)/(elcTempProfile(x)*eV))
end

function ionTempProfile(x)
  return ionTemp
end

function fIonProfile(x,y,z,v,mu)
  return (2*math.pi*ionTempProfile(x)*eV/ionMass)^(-3/2)*
    math.exp(-ionMass*v^2/(2*ionTempProfile(x)*eV))*
    math.exp(-math.abs(mu)*bFieldProfile(x)/(ionTempProfile(x)*eV))
end

-- set random seed based on processor number to avoid imprinting
-- domain decomp on ICs.
r = Lucee.getRank()
math.randomseed(100000*r+os.time())
function perturbDensityProfile(x,y,z)
  return 1e-3*math.random(-1,1)
end

-- Initial density perturbation
nInitialPerturbCG = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numExclusiveNodes(),
   ghost = {1, 1},
}
nInitialPerturbDG = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- initialize perturbation to electron distribution function
-- does not contain f_0, instead this will be multiplied with f_0
initDensityPerturb = Updater.EvalOnNodes3D {
   onGrid = grid_3d,
   basis = basis_3d,
   shareCommonNodes = true,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 return 1 + perturbDensityProfile(x,y,z)
	 end
}
runUpdater(initDensityPerturb, 0.0, 0.0, {}, {nInitialPerturbCG})
nInitialPerturbCG:sync()
-- copy nInitialPerturb to a DG field
-- create updater to copy a continuous 3D field to a 3D discontinuous field
copyCToD3d = Updater.CopyContToDisCont3D {
   onGrid = grid_3d,
   basis = basis_3d,
}
runUpdater(copyCToD3d, 0.0, 0.0, {nInitialPerturbCG}, {nInitialPerturbDG})
-- Copy 3d dg density perturbation to a 5d field
fElcInitialPerturb = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
fIonInitialPerturb = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
-- to copy phi to the 5d space the electron hamiltonian is stored on
copy3dTo5dElc = Updater.CopyNodalFields3D_5D {
  -- 5D phase-space grid 
  onGrid = grid_elc_5d,
  -- 5D phase-space basis functions
  targetBasis = basis_elc_5d,
  -- 3D spatial basis functions
  sourceBasis = basis_3d,
}
-- to copy phi to the 5d space the ion hamiltonian is stored on
copy3dTo5dIon = Updater.CopyNodalFields3D_5D {
  -- 5D phase-space grid 
  onGrid = grid_ion_5d,
  -- 5D phase-space basis functions
  targetBasis = basis_ion_5d,
  -- 3D spatial basis functions
  sourceBasis = basis_3d,
}
runUpdater(copy3dTo5dElc, 0.0, 0.0, {nInitialPerturbDG}, {fElcInitialPerturb})
-- do the same initial perturbation calculation for ions
--runUpdater(initDensityPerturb, 0.0, 0.0, {}, {nInitialPerturbCG})
--runUpdater(copyCToD3d, 0.0, 0.0, {nInitialPerturbCG}, {nInitialPerturbDG})
runUpdater(copy3dTo5dIon, 0.0, 0.0, {nInitialPerturbDG}, {fIonInitialPerturb})
-- to multiply two 5d fields together
multiply5dCalcElc = Updater.FieldArithmeticUpdater5D {
  onGrid = grid_elc_5d,
  basis = basis_elc_5d,
  evaluate = function(fld1, fld2)
    return fld1*fld2
  end
}
multiply5dCalcIon = Updater.FieldArithmeticUpdater5D {
  onGrid = grid_ion_5d,
  basis = basis_ion_5d,
  evaluate = function(fld1, fld2)
    return fld1*fld2
  end
}

-- initialize electron distribution function
initfElc = Updater.ProjectOnNodalBasis5D {
   onGrid = grid_elc_5d,
   basis = basis_elc_5d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,v,mu,t)
		 return n0*fElcProfile(x,y,z,v,mu)
	 end
}
runUpdater(initfElc, 0.0, 0.0, {}, {fElc})
fElc:sync()

-- to scale distribution function
scaleInitDistFElc = Updater.ETGInitializeDensity5D {
  -- 5D phase-space grid 
   onGrid = grid_elc_5d,
   -- 5D phase-space basis functions
   basis5d = basis_elc_5d,
   -- 3D spatial basis functions
   basis3d = basis_3d,
   -- Desired constant density
   constantDensity = n0,
   polyOrder = polyOrder,
}
-- to scale distribution function
scaleInitDistFIon = Updater.ETGInitializeDensity5D {
  -- 5D phase-space grid 
   onGrid = grid_ion_5d,
   -- 5D phase-space basis functions
   basis5d = basis_ion_5d,
   -- 3D spatial basis functions
   basis3d = basis_3d,
   -- Desired constant density
   constantDensity = n0,
   polyOrder = polyOrder,
}

-- initialize ion distribution function
initfIon = Updater.ProjectOnNodalBasis5D {
   onGrid = grid_ion_5d,
   basis = basis_ion_5d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,v,mu,t)
		 return n0*fIonProfile(x,y,z,v,mu)
	 end
}
runUpdater(initfIon, 0.0, 0.0, {}, {fIon})
fIon:sync()

function sourceTempProfile(x)
  return sourceTempInit*(math.exp(-(x-sourceX_UH)^2/sourceLambda_UH^2) + 
      math.exp(-(x-sourceX_EC)^2/sourceLambda_EC^2) + 0.1)
end

function sourceDensityProfile(x,y,z)
  return (0.95*sourceAmplitude_UH*math.exp(-(x-sourceX_UH)^2/sourceLambda_UH^2) + 
      0.95*sourceAmplitude_EC*math.exp(-(x-sourceX_EC)^2/sourceLambda_EC^2) + 
      0.05*sourceAmplitude_UH + 0.05*sourceAmplitude_EC)
end

ionSource = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
-- electron and ion sources
initSourceIon = Updater.EvalOnNodes5D {
   onGrid = grid_ion_5d,
   basis = basis_ion_5d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,v,mu,t)
		 return math.exp(-ionMass*v^2/(2*sourceTempProfile(x)*eV))*
      math.exp(-math.abs(mu)*bFieldProfile(x)/(sourceTempProfile(x)*eV))
	 end
}
runUpdater(initSourceIon, 0.0, 0.0, {}, {ionSource})

-- to scale source term
scaleInitSourceIon = Updater.SOLInitializeDensity {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/ionMass,
  evaluate = function(x,y,z,t)
    return sourceDensityProfile(x,y,z)
  end
}

elcSource = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
initSourceElc = Updater.EvalOnNodes5D {
   onGrid = grid_elc_5d,
   basis = basis_elc_5d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,v,mu,t)
		 return math.exp(-elcMass*v^2/(2*sourceTempProfile(x)*eV))*
      math.exp(-math.abs(mu)*bFieldProfile(x)/(sourceTempProfile(x)*eV))
	 end
}
runUpdater(initSourceElc, 0.0, 0.0, {}, {elcSource})

-- to scale source term
scaleInitSourceElc = Updater.SOLInitializeDensity {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/elcMass,
  evaluate = function(x,y,z,t)
    return sourceDensityProfile(x,y,z)
  end
}

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
  evaluate = function (x,y,z,t)
    return bFieldProfile(x)
  end
}
runUpdater(initb, 0.0, 0.0, {}, {bField3d})
bField3d:sync()
-- Jacobian Factor for electrons (5D)
jacobianFieldElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
initJacobianElc = Updater.EvalOnNodes5D {
  onGrid = grid_elc_5d,
  basis = basis_elc_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return elcMass*elcMass*bFieldProfile(x)
  end
}
-- Fill out jacobian
runUpdater(initJacobianElc, 0.0, 0.0, {}, {jacobianFieldElc})
jacobianFieldElc:sync()
-- Jacobian Factor for ions (5D)
jacobianFieldIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
initJacobianIon = Updater.EvalOnNodes5D {
  onGrid = grid_ion_5d,
  basis = basis_ion_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return ionMass*ionMass*bFieldProfile(x)
  end
}
-- Fill out jacobian
runUpdater(initJacobianIon, 0.0, 0.0, {}, {jacobianFieldIon})
jacobianFieldIon:sync()

-- B_y^* field for electrons (5D)
bStarYFieldElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarY for electrons
initbStarYElc = Updater.EvalOnNodes5D {
  onGrid = grid_elc_5d,
  basis = basis_elc_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return -elcMass*vPara/(elcCharge*x)
  end
}
runUpdater(initbStarYElc, 0.0, 0.0, {}, {bStarYFieldElc})
bStarYFieldElc:sync()
-- B_y^* field for ions (5D)
bStarYFieldIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarY for ions
initbStarYIon = Updater.EvalOnNodes5D {
  onGrid = grid_ion_5d,
  basis = basis_ion_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return -ionMass*vPara/(ionCharge*x)
  end
}
runUpdater(initbStarYIon, 0.0, 0.0, {}, {bStarYFieldIon})
bStarYFieldIon:sync()
-- B_z^* field for electrons (5D)
bStarZFieldElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarZ for electrons
initbStarZElc = Updater.EvalOnNodes5D {
  onGrid = grid_elc_5d,
  basis = basis_elc_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return bFieldProfile(x)
  end
}
runUpdater(initbStarZElc, 0.0, 0.0, {}, {bStarZFieldElc})
bStarZFieldElc:sync()
-- B_z^* field for ions (5D)
bStarZFieldIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarZ for ions
initbStarZIon = Updater.EvalOnNodes5D {
  onGrid = grid_ion_5d,
  basis = basis_ion_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return bFieldProfile(x)
  end
}
runUpdater(initbStarZIon, 0.0, 0.0, {}, {bStarZFieldIon})
bStarZFieldIon:sync()

-- Zero inflow boundary conditions for ions
-- Outflow BCs
function getRepTbl(pOrder, val)
  if pOrder == 1 then
    -- return 32 of these
    return {val, val, val, val, val, val, val, val,
            val, val, val, val, val, val, val, val,
            val, val, val, val, val, val, val, val,
            val, val, val, val, val, val, val, val}
  end
end

function getCountTbl(pOrder, val)
  if pOrder == 1 then
    return {0,1,2,3,4,5,6,7,
            8,9,10,11,12,13,14,15,
            16,17,18,19,20,21,22,23,
            24,25,26,27,28,29,30,31}
  end
end

bcConst = BoundaryCondition.Const { 
  components = getCountTbl(polyOrder),
  values = getRepTbl(polyOrder, 0.0),
}

bcLowerZIon = Updater.Bc5D {
  onGrid = grid_ion_5d,
  boundaryConditions = {bcConst},
  dir = 2,
  edge = "lower",
}

bcUpperZIon = Updater.Bc5D {
  onGrid = grid_ion_5d,
  boundaryConditions = {bcConst},
  dir = 2,
  edge = "upper",
}

bcLowerZElc = Updater.Bc5D {
  onGrid = grid_elc_5d,
  boundaryConditions = {bcConst},
  dir = 2,
  edge = "lower",
}

bcUpperZElc = Updater.Bc5D {
  onGrid = grid_elc_5d,
  boundaryConditions = {bcConst},
  dir = 2,
  edge = "upper",
}
-- to compute ion flux at domain boundaries in z
edgeMomentsIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
edgeMomentCalcIon = Updater.MomentsAtEdges5D {
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  onGrid = grid_ion_5d,
  polyOrder = polyOrder,
  integrateGhosts = true,
}
-- to compute electron flux at domain boundaries in z
edgeMomentsElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
edgeMomentsDelta = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
edgeMomentCalcElc = Updater.MomentsAtEdges5D {
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  onGrid = grid_elc_5d,
  polyOrder = polyOrder,
  integrateGhosts = true,
}
-- to store discontinuous sheath potential
phiSLower = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
phiSUpper = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- to average phiSLower and phiSUpper
phiSAverage = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- to compute sheath potential on lower and upper surfaces in z (reflect electrons)
biasedSheathElcCalc = Updater.BiasedSheath5D {
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  onGrid = grid_elc_5d,
  polyOrder = polyOrder,
  speciesMass = elcMass,
  speciesCharge = -eV,
}
-- to compute sheath potential on lower and upper surfaces in z (reflect ions)
biasedSheathIonCalc = Updater.BiasedSheath5D {
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  onGrid = grid_ion_5d,
  polyOrder = polyOrder,
  speciesMass = ionMass,
  speciesCharge = eV,
}

-- moments of distribution function
-- electron guiding center number density
numDensityElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron first parallel velocity moment
mom1dir3Elc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron second parallel velocity moment
mom2dir3Elc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron first mu moment
mom1dir4Elc = DataStruct.Field3D {
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
-- ion guiding center number density
numDensityIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- ion first parallel velocity moment
mom1dir3Ion = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- ion second parallel velocity moment
mom2dir3Ion = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- ion first mu moment
mom1dir4Ion = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- used to calculate number density without B weight
oneFieldElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
oneFieldElc:clear(1.0)
oneFieldIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
oneFieldIon:clear(1.0)
-- used to store number density without B weight
kineticEnergyIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
kineticEnergyElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- ion total tempearture
temperatureIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- store difference between elc and ion gc densities
numDensityDelta = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to compute number density for electrons
numDensityElcCalc = Updater.DistFuncMomentCalcWeighted3D {
   onGrid = grid_elc_5d,
   basis5d = basis_elc_5d,
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}

generalIntegralElcCalc = Updater.SOLGeneralIntegralAtNodeCalc { --WeightedProjectionCalc {
   onGrid = grid_elc_5d,
   basis5d = basis_elc_5d,
   basis3d = basis_3d,
   basis2d = basis_2d,
   scaleFactor = 2*math.pi/elcMass,
}

generalIntegralIonCalc = Updater.SOLGeneralIntegralAtNodeCalc { --WeightedProjectionCalc {
   onGrid = grid_ion_5d,
   basis5d = basis_ion_5d,
   basis3d = basis_3d,
   basis2d = basis_2d,
   scaleFactor = 2*math.pi/ionMass,
}

-- to compute first parallel velocity moment for electrons
mom1VParaElcCalc = Updater.DistFuncMomentCalcWeighted3D {
   onGrid = grid_elc_5d,
   basis5d = basis_elc_5d,
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 1,
   -- moment direction
   momentDirection = 3,
}
-- to compute second parallel velocity moment for electrons
mom2VParaElcCalc = Updater.DistFuncMomentCalcWeighted3D {
   onGrid = grid_elc_5d,
   basis5d = basis_elc_5d,
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 2,
   -- moment direction
   momentDirection = 3,
}
-- to compute first mu moment for electrons
mom1MuElcCalc = Updater.DistFuncMomentCalcWeighted3D {
   onGrid = grid_elc_5d,
   basis5d = basis_elc_5d,
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 1,
   -- moment direction
   momentDirection = 4,
}
-- to compute temperature for electrons
temperatureElcCalc = Updater.SOLTemperatureAtNodeCalc {
   onGrid = grid_3d,
   basis = basis_3d,
   speciesMass = elcMass,
}
-- to compute number density for ions
numDensityIonCalc = Updater.DistFuncMomentCalcWeighted3D {
   onGrid = grid_ion_5d,
   basis5d = basis_ion_5d,
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}
-- to compute first parallel velocity moment for electrons
mom1VParaIonCalc = Updater.DistFuncMomentCalcWeighted3D {
   onGrid = grid_ion_5d,
   basis5d = basis_ion_5d,
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 1,
   -- moment direction
   momentDirection = 3,
}
-- to compute second parallel velocity moment for electrons
mom2VParaIonCalc = Updater.DistFuncMomentCalcWeighted3D {
   onGrid = grid_ion_5d,
   basis5d = basis_ion_5d,
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 2,
   -- moment direction
   momentDirection = 3,
}
-- to compute first mu moment for electrons
mom1MuIonCalc = Updater.DistFuncMomentCalcWeighted3D {
   onGrid = grid_ion_5d,
   basis5d = basis_ion_5d,
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 1,
   -- moment direction
   momentDirection = 4,
}
-- to compute temperature for electrons
temperatureIonCalc = Updater.SOLTemperatureAtNodeCalc {
   onGrid = grid_3d,
   basis = basis_3d,
   speciesMass = ionMass,
}

-- to calculate the average of a 3d field
fieldInt3d = Updater.IntegrateGeneralField3D {
  onGrid = grid_3d,
  basis = basis_3d,
}
-- stores average of ion density
avgIonDensity = DataStruct.DynVector { numComponents = 1, }
avgElcDensity = DataStruct.DynVector { numComponents = 1, }
totalIonCounter = DataStruct.DynVector { numComponents = 1, }
totalElcCounter = DataStruct.DynVector { numComponents = 1, }

-- to calculate the average of the sheath potential
fieldInt2d = Updater.IntegrateGeneralField2D {
  onGrid = grid_2d,
  basis = basis_2d,
}
-- stores average of ion density
avgPhiSLower = DataStruct.DynVector { numComponents = 1, }

-- to help with positivity issues
positivityElcUpdater = Updater.SOLPositivityUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
}
positivityIonUpdater = Updater.SOLPositivityUpdater {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
}

localPositivityIonUpdater = Updater.SOLLocalPositivityUpdater {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
}

localPositivityElcUpdater = Updater.SOLLocalPositivityUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
}

positivityDragIonUpdater = Updater.SOLPositivityDragCellUpdater {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
}

positivityDragElcUpdater = Updater.SOLPositivityDragCellUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
}

positivityScaleIonUpdater = Updater.SOLPositivityScaleCellUpdater {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
}

positivityScaleElcUpdater = Updater.SOLPositivityScaleCellUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
}

-- diagnostic for flux across edge
totalEdgeFluxElcCalc = Updater.SOLFluxAcrossEdgeCalc {
  onGrid = grid_elc_5d,
  basis = basis_elc_5d,
  scaleFactor = 2*math.pi/elcMass^3,
  integrateGhosts = true,
}
totalEdgeFluxElc = DataStruct.DynVector { numComponents = 2, }
-- diagnostic for flux across edge
totalEdgeFluxIonCalc = Updater.SOLFluxAcrossEdgeCalc {
  onGrid = grid_ion_5d,
  basis = basis_ion_5d,
  scaleFactor = 2*math.pi/ionMass^3,
  integrateGhosts = true,
}
totalEdgeFluxIon = DataStruct.DynVector { numComponents = 2, }

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
  basis5d = basis_elc_5d,
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

-- dynvectors to store average heat flux on each surface
totalHeatFluxUpperElc = DataStruct.DynVector { numComponents = 1, }
totalHeatFluxLowerElc = DataStruct.DynVector { numComponents = 1, }
totalHeatFluxUpperIon = DataStruct.DynVector { numComponents = 1, }
totalHeatFluxLowerIon = DataStruct.DynVector { numComponents = 1, }

-- diagnostic to look at phi3d in a single cell at all times
recordFieldInCellPhiCalc = Updater.RecordFieldInCell3D {
  onGrid = grid_3d,
  cellIndex = {0,0,0},
}
fieldInCellPhi = DataStruct.DynVector { numComponents = basis_3d:numNodes(), }
fieldInCellNumElc = DataStruct.DynVector { numComponents = basis_3d:numNodes(), }
fieldInCellU = DataStruct.DynVector { numComponents = basis_3d:numNodes(), }

maxwellianSearchTarget = DataStruct.DynVector { numComponents = basis_3d:numNodes(), }
maxwellianSearchGuess = DataStruct.DynVector { numComponents = basis_3d:numNodes(), }
maxwellianSearchActual = DataStruct.DynVector { numComponents = basis_3d:numNodes(), }

-- total energy vs time dynvectors
elcEnergyHistory = DataStruct.DynVector { numComponents = 1, }
ionEnergyHistory = DataStruct.DynVector { numComponents = 1, }

-- to calculate ion energy integrated over a cell
ionEnergyAtCellCalc = Updater.SOLEnergyAtCellCalc {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  scaleFactor = 2*math.pi/ionMass^3
}

-- to calculate elc energy integrated over a cell
elcEnergyAtCellCalc = Updater.SOLEnergyAtCellCalc {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  scaleFactor = 2*math.pi/elcMass^3
}

ionEnergyAtCellAfterPositivity = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
ionEnergyAtCellBeforePositivity = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
ionEnergyAtCellAfterMaxDrag = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
elcEnergyAtCellAfterPositivity = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
elcEnergyAtCellBeforePositivity = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
elcEnergyAtCellAfterMaxDrag = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

ionSourceTempHistory = DataStruct.DynVector { numComponents = 1, }
elcSourceTempHistory = DataStruct.DynVector { numComponents = 1, }
ionTempHistory = DataStruct.DynVector { numComponents = 1, }
elcTempHistory = DataStruct.DynVector { numComponents = 1, }

-- collision operator code
function getElcAlpha(t)
  -- see http://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node35.html
  local n_e0 = avgElcDensity:lastInsertedData()/totalVol
  local logLambda = 6.6-0.5*math.log(n_e0/10^(20)) + 1.5*math.log(elcTemp)
  -- alpha needs to be multiplied by average n/(kTe)^(3/2) in C++ code
  return logLambda*elcCharge^4/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(elcMass))
end

maxwellianEval = Updater.SOLMaxwellianAtNodeCalc {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/elcMass,
  speciesMass = elcMass,
}

fElcCollisions = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}

fElcMaxwellian = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}

fElcMaxwellianForward = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}

mom1dir3ElcMaxwellianInput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
mom1dir3ElcMaxwellianForwardInput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
mom1dir3ElcMaxwellianIncrement = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
mom1dir3ElcMaxwellianNumerical = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
mom1dir3ElcIncrementOneNumerical = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
mom1dir3ElcIncrementTwoNumerical = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
mom1dir3ElcMaxwellianOutput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureElcMaxwellianInput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureElcMaxwellianForwardInput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureElcMaxwellianIncrement = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureElcMaxwellianNumerical = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureElcIncrementOneNumerical = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureElcIncrementTwoNumerical = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureElcMaxwellianOutput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

bgkCollisionOperator = Updater.SOLBGKCollisionUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  alpha = function(t)
    return getElcAlpha(t)
  end,
  onlyIncrement = true,
}

maxwellianNewtonSearch = Updater.SOLMaxwellianParameterNewtonCalc {
  onGrid = grid_3d,
  basis3d = basis_3d,
  epsilonU = epsilonU,
  epsilonT = epsilonT,
}

function applyPositivity(tCurr, myDt, distfElc, distfIon)  -- apply positivity to ions
  local ionStatus = runUpdater(positivityIonUpdater, tCurr, myDt, {bField3d}, {distfIon})
  -- apply positivity to elcs
  local elcStatus = runUpdater(positivityElcUpdater, tCurr, myDt, {bField3d}, {distfElc})

  if (elcStatus == false or ionStatus == false) then
    Lucee.logInfo(string.format("-- ending execution because an entire cell was negative --"))
    endExecution = true
    return false
  end
  
  return true
end

function calcHamiltonianElc(phiIn, hamilOut)
  -- copy phi3d into a 5d field
  runUpdater(copy3dTo5dElc, 0.0, 0.0, {phiIn}, {hamilOut})
  hamilOut:scale(elcCharge)
  -- add the kinetic energy component
  hamilOut:accumulate(1.0, hamilKeElc)
  hamilOut:sync()
end

function calcHamiltonianIon(phiIn, hamilOut)
  -- copy phi3d into a 5d field
  runUpdater(copy3dTo5dIon, 0.0, 0.0, {phiIn}, {hamilOut})
  hamilOut:scale(ionCharge)
  -- add the kinetic energy component
  hamilOut:accumulate(1.0, hamilKeIon)
  hamilOut:sync()
end

function applyBcToDistF(tCurr, myDt, fElc, fIon, phiIn)
  -- apply periodic bc in y to ions
  fIon:sync()
  -- apply periodic bc in y to electrons
  fElc:sync()
  -- apply zero inflow boundary conditions in z to ions (not sure if needed)
  runUpdater(bcLowerZIon, tCurr, myDt, {}, {fIon})
  runUpdater(bcUpperZIon, tCurr, myDt, {}, {fIon})
  -- apply zero inflow boundary conditions in z to electrons (not sure if needed)
  runUpdater(bcLowerZElc, tCurr, myDt, {}, {fElc})
  runUpdater(bcUpperZElc, tCurr, myDt, {}, {fElc})
  -- apply biased sheath boundary conditions to electron distribution function
  sheathStatus = runUpdater(biasedSheathElcCalc, tCurr, myDt, {phiIn, hamilDerivKeElc}, {fElc})
  -- apply biased sheath boundary conditions to ion distribution function
  sheathStatus = runUpdater(biasedSheathIonCalc, tCurr, myDt, {phiIn, hamilDerivKeIon}, {fIon})
  -- sync cells in distribution function (not sure if needed)
  fElc:sync()
  fIon:sync()
end

function calcCollisionsElc(tCurr, myDt, fInitial, fFinal)
  -- fFinal f(t_{n+1}) is what to add increment to
  -- fInitial is f(t_n)

  -- compute number density, temperature, and mom1dir3 of input
  runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fInitial, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fInitial, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fInitial, hamilKeElc, bField3d}, {kineticEnergyElc})
  local negativeTemperature = runUpdater(temperatureElcCalc, 0.0, 0.0, {numDensityElc,kineticEnergyElc,mom1dir3Elc},
    {temperatureElc})
  if (negativeTemperature == false) then
    endExecution = true
    Lucee.logInfo(string.format("-- Ending execution because negative temperature. --"))
    return
  end
  
  -- starting guess for Newton's method will be previous iteration final state
  local convergenceStatus = false
  local iterCount = 0
  local maxIter = 40
  local maxwellianValid = true

  -- loop to find maxwellian with desired density, <v>, and temperature
  while convergenceStatus == false and iterCount < maxIter do

    Lucee.logInfo (string.format("-- Iteration %d", iterCount))

    -- evaluate numerical maxwellian at current point
    maxwellianValid = runUpdater(maxwellianEval, 0.0, 0.0, {numDensityElc, mom1dir3ElcMaxwellianInput,
      temperatureElcMaxwellianInput, bField3d}, {fElcMaxwellian})

    if maxwellianValid == false then
      Lucee.logInfo (string.format("-- Maxwellian invalid"))
      endExecution = true
      break
    end
    fElcMaxwellian:sync()
    -- compute parameters of maxwellian (it should already have correct number density)
    runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElcMaxwellian, hamilDerivKeElc, bField3d},
      {mom1dir3ElcMaxwellianNumerical})
    runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElcMaxwellian, hamilKeElc, bField3d},
      {kineticEnergyElc})
    -- compute elc temperature
    runUpdater(temperatureElcCalc, 0.0, 0.0, {numDensityElc,kineticEnergyElc,mom1dir3ElcMaxwellianNumerical},
      {temperatureElcMaxwellianNumerical})

    -- increment mom1dir3ElcMaxwellianInput by epsilonU factor 
    mom1dir3ElcMaxwellianIncrement:copy(mom1dir3ElcMaxwellianInput)
    mom1dir3ElcMaxwellianIncrement:abs(1.0) -- if field is 0.0 at any point, replace by 1.0
    mom1dir3ElcMaxwellianIncrement:scale(epsilonU)
    -- increment temperatureElcMaxwellianInput by epsilonT factor
    temperatureElcMaxwellianIncrement:copy(temperatureElcMaxwellianInput)
    temperatureElcMaxwellianIncrement:abs()
    temperatureElcMaxwellianIncrement:scale(epsilonT)

    -- compute n*u and temperature for increment in mom1dir3
    mom1dir3ElcMaxwellianForwardInput:copy(mom1dir3ElcMaxwellianInput)
    mom1dir3ElcMaxwellianForwardInput:accumulate(1.0, mom1dir3ElcMaxwellianIncrement)
    temperatureElcMaxwellianForwardInput:copy(temperatureElcMaxwellianInput)
    -- evaluate numerical maxwellian at current point
    runUpdater(maxwellianEval, 0.0, 0.0, {numDensityElc, mom1dir3ElcMaxwellianForwardInput,
      temperatureElcMaxwellianForwardInput, bField3d}, {fElcMaxwellianForward})
    fElcMaxwellianForward:sync()
    -- compute parameters of maxwellian (it should already have correct number density)
    runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElcMaxwellianForward, hamilDerivKeElc, bField3d},
      {mom1dir3ElcIncrementOneNumerical})
    runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElcMaxwellianForward, hamilKeElc, bField3d},
      {kineticEnergyElc})
    -- compute elc temperature
    runUpdater(temperatureElcCalc, 0.0, 0.0, {numDensityElc,kineticEnergyElc,mom1dir3ElcIncrementOneNumerical},
      {temperatureElcIncrementOneNumerical})

    -- compute n*u and temperature for increment in temperature
    mom1dir3ElcMaxwellianForwardInput:copy(mom1dir3ElcMaxwellianInput)
    temperatureElcMaxwellianForwardInput:copy(temperatureElcMaxwellianInput)
    temperatureElcMaxwellianForwardInput:accumulate(1.0, temperatureElcMaxwellianIncrement)
    -- evaluate numerical maxwellian at current point
    runUpdater(maxwellianEval, 0.0, 0.0, {numDensityElc, mom1dir3ElcMaxwellianForwardInput,
      temperatureElcMaxwellianForwardInput, bField3d}, {fElcMaxwellianForward})
    fElcMaxwellianForward:sync()
    -- compute parameters of maxwellian (it should already have correct number density)
    runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElcMaxwellianForward, hamilDerivKeElc, bField3d},
      {mom1dir3ElcIncrementTwoNumerical})
    runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElcMaxwellianForward, hamilKeElc, bField3d},
      {kineticEnergyElc})
    -- compute elc temperature
    runUpdater(temperatureElcCalc, 0.0, 0.0, {numDensityElc,kineticEnergyElc,mom1dir3ElcIncrementTwoNumerical},
      {temperatureElcIncrementTwoNumerical})

    -- compute input parameters for next iteration using a newton's method
    convergenceStatus = runUpdater(maxwellianNewtonSearch, 0.0, 0.0, {mom1dir3Elc, temperatureElc,
      mom1dir3ElcMaxwellianNumerical,temperatureElcMaxwellianNumerical,
      mom1dir3ElcIncrementOneNumerical,temperatureElcIncrementOneNumerical,
      mom1dir3ElcIncrementTwoNumerical,temperatureElcIncrementTwoNumerical,
      mom1dir3ElcMaxwellianInput,temperatureElcMaxwellianInput},
      {mom1dir3ElcMaxwellianOutput, temperatureElcMaxwellianOutput})

    mom1dir3ElcMaxwellianInput:copy(mom1dir3ElcMaxwellianOutput)
    temperatureElcMaxwellianInput:copy(temperatureElcMaxwellianOutput)
    iterCount = iterCount + 1
  end

  if convergenceStatus == false then
    Lucee.logInfo (string.format("-- Unable to find Maxwellian after %d iterations.", iterCount))
    --endExecution = true
    --return
  end

  -- calculate total electron density
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityElc}, {avgElcDensity})
  -- unneeded?
  fElcCollisions:clear(0.0)
  runUpdater(bgkCollisionOperator, 0.0, 0.0, {fInitial, fElcMaxwellian, numDensityElc,
    temperatureElcMaxwellianNumerical}, {fElcCollisions})
  fFinal:accumulate(myDt, fElcCollisions)
end

endExecution = false
-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
  local netDtSuggested = myDt
  -- RK stage 1
  -- add on collisions
  f1Ion:copy(fIon)
  f1Elc:copy(fElc)
  calcCollisionsElc(tCurr, myDt, fElc, f1Elc)
  local positivityStatus = applyPositivity(tCurr, myDt, f1Elc, f1Ion)
  if (positivityStatus == false) then
    fElc:copy(f1Elc)
    fIon:copy(f1Ion)
    return false, netDtSuggested
  end
  
  if (endExecution == true) then
    endExecution = true
    calcContinuousDiagnostics(tCurr, myDt)
    return false, netDtSuggested
  end

  -- RK stage 2
  -- add on collisions
  fNewIon:copy(f1Ion)
  fNewElc:copy(f1Elc)
  calcCollisionsElc(tCurr, myDt, f1Elc, fNewElc)
  positivityStatus = applyPositivity(tCurr, myDt, fNewElc, fNewIon)
  if (positivityStatus == false) then
    fElc:copy(fNewElc)
    fIon:copy(fNewIon)
    return false, netDtSuggested
  end
  
  f1Elc:combine(3.0/4.0, fElc, 1.0/4.0, fNewElc)
  f1Ion:combine(3.0/4.0, fIon, 1.0/4.0, fNewIon)

  if (endExecution == true) then
    endExecution = true
    calcContinuousDiagnostics(tCurr, myDt)
    return false, math.min(myDtSuggestedElc, myDtSuggestedIon)
  end
  -- RK stage 3
  -- add on particle source and collisions
  fNewIon:copy(f1Ion)
  fNewElc:copy(f1Elc)
  calcCollisionsElc(tCurr, myDt, f1Elc, fNewElc)
  
  positivityStatus = applyPositivity(tCurr, myDt, fNewElc, fNewIon)
  if (positivityStatus == false) then
    fElc:copy(fNewElc)
    fIon:copy(fNewIon)
    return false, math.min(myDtSuggestedElc, myDtSuggestedIon)
  end
  f1Elc:combine(1.0/3.0, fElc, 2.0/3.0, fNewElc)
  f1Ion:combine(1.0/3.0, fIon, 2.0/3.0, fNewIon)

  -- solve for potential
  if (endExecution == true) then
    endExecution = true
    calcContinuousDiagnostics(tCurr, myDt)
    return false, math.min(myDtSuggestedElc, myDtSuggestedIon)
  end

  fElc:copy(f1Elc)
  fIon:copy(f1Ion)

  -- calculate continuous time diagnostics
  calcContinuousDiagnostics(tCurr, myDt)

  return true, netDtSuggested
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
  local step = 1
  local tCurr = tStart
  local myDt = initDt
  local status, dtSuggested

  while tCurr<=tEnd do
    -- Store fields that might need to be reverted
    fDupElc:copy(fElc)
    fDupIon:copy(fIon)
    hamilDupElc:copy(hamilElc)
    hamilDupIon:copy(hamilIon)

    -- if needed adjust dt to hit tEnd exactly
    if (tCurr+myDt > tEnd) then
      myDt = tEnd-tCurr
    end

    Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
    status, dtSuggested = rk3(tCurr, myDt)

    if (status == false) then
      if (endExecution == true) then
        --calcDiagnostics(tCurr, myDt)
        tCurr = tCurr + myDt
        Lucee.logInfo (string.format("** Something failed, ending execution "))
        break
      else
        -- time-step too large
        Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
        myDt = dtSuggested
        -- Revert fields to previous value
        fElc:copy(fDupElc)
        fIon:copy(fDupIon)
        hamilElc:copy(hamilDupElc)
        hamilIon:copy(hamilDupIon)
      end
    else
      calcDiagnostics(tCurr, myDt)
      tCurr = tCurr + myDt
      myDt = dtSuggested
      step = step + 1

      if (tCurr >= tEnd) then
        break
      end
    end
  end

  return dtSuggested, tCurr
end

-- function containing diagnostics that must be called at the end of every frame
function calcDiagnostics(tCurr, myDt)
  -- compute electron moments
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityElc}, {avgElcDensity})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilKeElc, bField3d}, {kineticEnergyElc})
  -- compute electron temperature
  negativeTemperature = runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc,
    kineticEnergyElc, mom1dir3Elc}, {temperatureElc})
  -- convert to ev
  temperatureElc:scale(1/eV)

  -- compute electron moments
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, oneFieldIon, bField3d}, {numDensityIon})
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityIon}, {avgIonDensity})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, hamilDerivKeIon, bField3d}, {mom1dir3Ion})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, hamilKeIon, bField3d}, {kineticEnergyIon})
  -- compute electron temperature
  negativeTemperature = runUpdater(temperatureIonCalc, tCurr, myDt, {numDensityIon,
    kineticEnergyIon, mom1dir3Ion}, {temperatureIon})
  -- convert to ev
  temperatureIon:scale(1/eV)

  -- compute average ion temperature (volume factor should be removed in post-processing)
  runUpdater(fieldInt3d, tCurr, myDt, {temperatureElc}, {elcTempHistory})
  runUpdater(fieldInt3d, tCurr, myDt, {temperatureIon}, {ionTempHistory})
end

-- function to compute source temperature
function calcSourceDiagnostics(tCurr, myDt, elcSource, ionSource)
  -- compute electron moments
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {elcSource, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {elcSource, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {elcSource, hamilKeElc, bField3d}, {kineticEnergyElc})
  -- compute electron temperature
  runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc,
    kineticEnergyElc, mom1dir3Elc}, {temperatureElc})
  -- convert to ev
  temperatureElc:scale(1/eV)

  -- compute electron moments
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {ionSource, oneFieldIon, bField3d}, {numDensityIon})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {ionSource, hamilDerivKeIon, bField3d}, {mom1dir3Ion})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {ionSource, hamilKeIon, bField3d}, {kineticEnergyIon})
  -- compute electron temperature
  runUpdater(temperatureIonCalc, tCurr, myDt, {numDensityIon,
    kineticEnergyIon, mom1dir3Ion}, {temperatureIon})
  -- convert to ev
  temperatureIon:scale(1/eV)

  numDensityElc:write( string.format("elcSource_%d.h5", 0), 0)
  numDensityIon:write( string.format("ionSource_%d.h5", 0), 0)

  -- compute average ion temperature
  local volume = (X_UPPER-X_LOWER)*(Y_UPPER-Y_LOWER)*(Z_UPPER-Z_LOWER)
  runUpdater(fieldInt3d, tCurr, myDt, {temperatureElc}, {elcSourceTempHistory})
  temperatureElc:write( string.format("elcSourceTemp_%d.h5", 0), 0)
  print(string.format( 'Average source elc temp = %g eV',elcSourceTempHistory:lastInsertedData()/volume ))
  runUpdater(fieldInt3d, tCurr, myDt, {temperatureIon}, {ionSourceTempHistory})
  temperatureIon:write( string.format("ionSourceTemp_%d.h5", 0), 0)
  print(string.format( 'Average source ion temp = %g eV',ionSourceTempHistory:lastInsertedData()/volume ))
end

-- function containing diagnostics that must be called at the end of every time step
function calcContinuousDiagnostics(tCurr, myDt)
  -- store phi in a cell
  --runUpdater(recordFieldInCellPhiCalc, tCurr, myDt, {phi3d}, {fieldInCellPhi})
  runUpdater(recordFieldInCellPhiCalc, tCurr, myDt, {numDensityElc}, {fieldInCellNumElc})
  runUpdater(recordFieldInCellPhiCalc, tCurr, myDt, {mom1dir3Elc}, {fieldInCellU})
  -- compute total guiding center density by integrating over whole domain
  -- note: includes a volume factor that should be removed in post-processing
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {totalElcCounter})
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityIon}, {totalIonCounter})
end

-- write data to H5 file
function writeFields(frameNum, tCurr)
  numDensityElc:write( string.format("nElc_%d.h5", frameNum), tCurr)
  --temperatureElc:write( string.format("tElc_%d.h5", frameNum), tCurr)
  numDensityIon:write( string.format("nIon_%d.h5", frameNum), tCurr)
  temperatureIon:write( string.format("tIon_%d.h5", frameNum), tCurr)
  --numDensityDelta:write( string.format("nDelta_%d.h5", frameNum), tCurr)
  
  --phi3d:write( string.format("phi_%d.h5", frameNum), tCurr)

  fElc:write( string.format("fElc_%d.h5", frameNum), tCurr)
  fIon:write( string.format("fIon_%d.h5", frameNum), tCurr)
  fElcMaxwellian:write( string.format("fElcMaxwellian_%d.h5", frameNum), tCurr)
  
  --edgeMomentsIon:write( string.format("edgeMomentsIon_%d.h5", frameNum), tCurr)
  --edgeMomentsElc:write( string.format("edgeMomentsElc_%d.h5", frameNum), tCurr)
  --edgeMomentsDelta:write( string.format("edgeMomentsDelta_%d.h5", frameNum), tCurr)

  totalIonCounter:write( string.format("totalIonCounter_%d.h5", frameNum), tCurr)
  totalElcCounter:write( string.format("totalElcCounter_%d.h5", frameNum), tCurr)
  --totalEdgeFluxIon:write( string.format("totalEdgeFluxIon_%d.h5", frameNum), tCurr)
  --totalEdgeFluxElc:write( string.format("totalEdgeFluxElc_%d.h5", frameNum), tCurr)

  --heatFluxUpperElc:write( string.format("heatFluxUpperElc_%d.h5", frameNum), tCurr)
  --heatFluxLowerElc:write( string.format("heatFluxLowerElc_%d.h5", frameNum), tCurr)
  --heatFluxUpperIon:write( string.format("heatFluxUpperIon_%d.h5", frameNum), tCurr)
  --heatFluxLowerIon:write( string.format("heatFluxLowerIon_%d.h5", frameNum), tCurr)

  elcTempHistory:write( string.format("elcTempHistory_%d.h5", frameNum), tCurr)
  ionTempHistory:write( string.format("ionTempHistory_%d.h5", frameNum), tCurr)

  mom1dir3Elc:write( string.format("mom1dir3Elc_%d.h5", frameNum), tCurr)
  mom2dir3Elc:write( string.format("mom2dir3Elc_%d.h5", frameNum), tCurr)
  mom1dir3Ion:write( string.format("mom1dir3Ion_%d.h5", frameNum), tCurr)

  --fieldInCellPhi:write( string.format("phiInCell_%d.h5", frameNum), tCurr)
  fieldInCellNumElc:write( string.format("nElcInCell_%d.h5", frameNum), tCurr)
  fieldInCellU:write( string.format("uInCell_%d.h5", frameNum), tCurr)
end

-- apply positivity preservation to the sources
applyPositivity(0.0, 0.0, elcSource, ionSource)
-- scale source terms to have correct density
runUpdater(scaleInitSourceIon, 0.0, 0.0, {bField3d}, {ionSource})
runUpdater(scaleInitSourceElc, 0.0, 0.0, {bField3d}, {elcSource})
-- ensure that source terms are in exact balance
runUpdater(numDensityElcCalc, 0.0, 0.0, {elcSource, bField3d}, {numDensityElc})
numDensityElc:scale(2*math.pi/elcMass)
-- store electron density temporarily in avgIonDensity structure
runUpdater(fieldInt3d, 0.0, 0.0, {numDensityElc}, {avgIonDensity})
totalSourceElc = avgIonDensity:lastInsertedData()
runUpdater(numDensityIonCalc, 0.0, 0.0, {ionSource, bField3d}, {numDensityIon})
numDensityIon:scale(2*math.pi/ionMass)
-- store ion density temporarily in avgIonDensity structure
runUpdater(fieldInt3d, 0.0, 0.0, {numDensityIon}, {avgIonDensity})
totalSourceIon = avgIonDensity:lastInsertedData()
-- scale ions to match electrons
ionSource:scale(totalSourceElc/totalSourceIon)
--print( string.format("ionSource scaled by %g (%g/%g)",totalSourceElc/totalSourceIon,totalSourceElc,totalSourceIon) )
ionSource:sync()
elcSource:sync()

-- write out some source fields
calcSourceDiagnostics(0.0, 0.0, elcSource, ionSource)

if Lucee.IsRestarting then
  fElc:read("fElc_" .. Lucee.RestartFrame .. ".h5")
  fElc:sync()
  fIon:read("fIon_" .. Lucee.RestartFrame .. ".h5")
  fIon:sync()

  startFrame = Lucee.RestartFrame + 1
  tCurr = tStart + tFrame*Lucee.RestartFrame
  runUpdater(numDensityElcCalc, 0.0, 0.0, {fElc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  runUpdater(numDensityIonCalc, 0.0, 0.0, {fIon, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- calculate new hamiltonians
  calcHamiltonianElc(phi3d, hamilElc)
  calcHamiltonianIon(phi3d, hamilIon)
  -- apply boundary conditions to distribution functions
  applyBcToDistF(0.0, 0.0, fElc, fIon, phi3d)
else
  -- Make sure distribution functions are positive
  applyPositivity(0.0, 0.0, fElc, fIon)
  -- first make sure background distribution functions have uniform density
  runUpdater(numDensityElcCalc, 0.0, 0.0, {fElc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  runUpdater(scaleInitDistFElc, 0.0, 0.0, {numDensityElc}, {fElc})
  runUpdater(numDensityIonCalc, 0.0, 0.0, {fIon, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  runUpdater(scaleInitDistFIon, 0.0, 0.0, {numDensityIon}, {fIon})
  -- multiply some small random density perturbations
  runUpdater(multiply5dCalcIon, 0.0, 0.0, {fIon, fIonInitialPerturb}, {fDupIon})
  fIon:copy(fDupIon)
  fIon:sync()
  runUpdater(multiply5dCalcElc, 0.0, 0.0, {fElc, fElcInitialPerturb}, {fDupElc})
  fElc:copy(fDupElc)
  fElc:sync()
  -- make sure total ions and electrons are equal
  runUpdater(numDensityElcCalc, 0.0, 0.0, {fElc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  -- store electron density temporarily in avgIonDensity structure
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityElc}, {avgIonDensity})
  totalInitElc = avgIonDensity:lastInsertedData()
  runUpdater(numDensityIonCalc, 0.0, 0.0, {fIon, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- store ion density temporarily in avgIonDensity structure
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityIon}, {avgIonDensity})
  totalInitIon = avgIonDensity:lastInsertedData()
  -- scale ions to match electrons
  fIon:scale(totalInitElc/totalInitIon)
  fIon:sync()
  fElc:sync()
  
  tCurr = tStart
  startFrame = 1
  runUpdater(numDensityElcCalc, 0.0, 0.0, {fElc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  runUpdater(numDensityIonCalc, 0.0, 0.0, {fIon, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
end

calcDiagnostics(tCurr, 0.0)
calcContinuousDiagnostics(tCurr, 0.0)
writeFields(startFrame-1,tCurr)

-- Newton's method initial guess
runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElc, oneFieldElc, bField3d}, {numDensityElc})
runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElc, hamilKeElc, bField3d}, {kineticEnergyElc})
runUpdater(temperatureElcCalc, 0.0, 0.0, {numDensityElc,kineticEnergyElc,mom1dir3Elc},
  {temperatureElc})
mom1dir3ElcMaxwellianInput:copy(mom1dir3Elc)
mom1dir3ElcMaxwellianInput:scale(0.5) -- half since this tends to overestimate
temperatureElcMaxwellianInput:copy(temperatureElc)

for frame = startFrame, nFrames do
  if endExecution == false then
    dtSuggested = 1e-7
    Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
    dtSuggested, tCurr = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
    writeFields(frame, tCurr)
    Lucee.logInfo ("")
  else
    --Lucee.logInfo (string.format("-- Not advancing solution from %g", tCurr))
  end
end
