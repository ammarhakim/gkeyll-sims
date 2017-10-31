-- Input file for a 3X2V SOL Simulation
-- 01-28-2016 Creation
-- 03-15-2016 First attempt at input file with kinetic elections and ions
-- Command line parameters
-- mpiexec -n 1 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i maxwellian.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist
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
-- 5-16-16: testing diffusion
-- 5-27-16: testing maxwellian search

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
tEnd = 30000e-6
dtSuggested = 10e-6 --0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 5
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

-- physical parameters
eV        = Lucee.ElementaryCharge
elcCharge = -eV
ionCharge = eV
eps0 = Lucee.Epsilon0 -- permittivity of free space
elcMass   = Lucee.ElectronMass
ionMass   = 2.014*Lucee.ProtonMass 
elcTemp   = 3.5 -- [eV]
ionTemp   = 3.5 -- [eV]
B0 = 0.076  -- [T]
n0 = 10^16 -- [1/m^3]
N = 4 -- field line turns
kPerpTimesRhoS = 0.2
sourceTemp     = 3.5 -- [eV] -- temperature to base grid off of
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
fIonMaxwellian = DataStruct.Field5D {
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
hamilDerivKeExactIon = DataStruct.Field5D {
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
  return B0--*R/x
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

-- compute v-parallel
initHamilKeDerivExactIon = Updater.EvalOnNodes5D {
   onGrid = grid_ion_5d,
   basis = basis_ion_5d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,vPara,mu,t)
      return vPara
   end
}
runUpdater(initHamilKeDerivExactIon, 0.0, 0.0, {}, {hamilDerivKeExactIon})


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
  if math.abs(v) < 2*vtIon and mu < ionMass*4*(vtIon*vtIon)/(2*B0) then
    return 100
  else
    return 10
  end
  --return math.exp(-ionMass*v^2/(2*ionTempProfile(x)*eV))*
  --  math.exp(-math.abs(mu)*bFieldProfile(x)/(ionTempProfile(x)*eV))
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
  return sourceTempInit
end

function sourceDensityProfile(x,y,z)
  return sourceAmplitude_UH
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
numDensityIonAlt = DataStruct.Field3D {
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
-- to compute first parallel moment and hamiltonian integral for ions
generalIntegralIonCalc = Updater.SOLWeightedProjectionCalc {
   onGrid = grid_ion_5d,
   basis5d = basis_ion_5d,
   basis3d = basis_3d,
   scaleFactor = 2*math.pi/ionMass,
}
-- to compute first parallel moment and hamiltonian integral for electrons
generalIntegralElcCalc = Updater.SOLWeightedProjectionCalc {
   onGrid = grid_elc_5d,
   basis5d = basis_elc_5d,
   basis3d = basis_3d,
   scaleFactor = 2*math.pi/elcMass,
}
-- fields for kinetic energy of ion and electrons at nodes
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

lenardBernsteinScaleIonUpdater = Updater.SOLPositivityScaleCellUpdater {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  positivityChecks = false,
}

positivityScaleElcUpdater = Updater.SOLPositivityScaleCellUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
}

-- diagnostic to look at phi3d in a single cell at all times
recordFieldInCellPhiCalc = Updater.RecordFieldInCell3D {
  onGrid = grid_3d,
  cellIndex = {1,1,1},
}
fieldInCellPhi = DataStruct.DynVector { numComponents = basis_3d:numNodes(), }

-- total energy vs time dynvectors
elcEnergyHistory = DataStruct.DynVector { numComponents = 1, }
ionEnergyHistory = DataStruct.DynVector { numComponents = 1, }
ionEnergyHistoryWrite = DataStruct.DynVector { numComponents = 1, }

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
ionEnergyAtCellAfterPositivity = DataStruct.Field3D {
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
elcEnergyAtCellAfterPositivity = DataStruct.Field3D {
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

-- lenard berstein code
function getIonAlpha(t)
  -- see http://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node35.html
  local n_i0 = avgIonDensity:lastInsertedData()/totalVol
  local logLambda = 6.6-0.5*math.log(n_i0/10^(20)) + 1.5*math.log(elcTemp)
  -- alpha needs to be multiplied by average n/(kTi)^(3/2) in C++ code
  return logLambda*ionCharge^4/(12*math.pi^(3/2)*eps0^2*math.sqrt(ionMass))
end
diffSlvrIon = Updater.LenardBernsteinDiff5DUpdater {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  cfl = cfl,
  onlyIncrement = true,
  speciesMass = ionMass,
  alpha = function(t)
    return getIonAlpha(t)
  end
}
dragSlvrIon = Updater.LenardBernsteinDrag5DUpdater {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  cfl = cfl,
  onlyIncrement = true,
  speciesMass = ionMass,
  alpha = function(t)
    return getIonAlpha(t)
  end
}

fIonCollisions = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
gVParIonCollisions = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
gMuIonCollisions = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
ionEnergyAtCellDrag = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
ionEnergyAtCellDiff = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

function applyPositivity(tCurr, myDt, distfElc, distfIon)
  -- apply positivity to ions
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

function applyPositivityAdvanced(tCurr, myDt, distfElc, distfIon)
  --Lucee.logInfo(string.format("-- ion positivity start --"))
  -- calculate initial ion energy
  runUpdater(ionEnergyAtCellCalc, tCurr, myDt, {distfIon, jacobianFieldIon, hamilParaIon},
    {ionEnergyAtCellBeforePositivity, ionEnergyHistory})
  -- apply positivity to ions
  local ionStatus = runUpdater(positivityIonUpdater, tCurr, myDt, {bField3d}, {distfIon})
  
  --Lucee.logInfo(string.format("-- elc positivity start --"))
  -- calculate initial elc energy
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilParaElc},
    {elcEnergyAtCellBeforePositivity, elcEnergyHistory})
  -- apply positivity to elcs
  local elcStatus = runUpdater(positivityElcUpdater, tCurr, myDt, {bField3d}, {distfElc})

  if (elcStatus == false or ionStatus == false) then
    Lucee.logInfo(string.format("-- ending execution because an entire cell was negative --"))
    endExecution = true
    return false
  end
  
  -- compute ion parallel energy at nodes after positivity is applied
  runUpdater(ionEnergyAtCellCalc, tCurr, myDt, {distfIon, jacobianFieldIon, hamilParaIon},
    {ionEnergyAtCellAfterPositivity})
  -- need to sync ghost cells for drag update
  distfIon:sync()
  runUpdater(positivityDragIonUpdater, tCurr, myDt, {distfIon, bField3d, ionEnergyAtCellAfterPositivity,
    ionEnergyAtCellBeforePositivity}, {fIonPositivityDrag})
  -- compute ion parallel energy at nodes after maximum drag step
  runUpdater(ionEnergyAtCellCalc, tCurr, myDt, {fIonPositivityDrag, jacobianFieldIon, hamilParaIon},
    {ionEnergyAtCellAfterMaxDrag})
  -- Compute the difference between distf post-positivity and distf post-drag
  fIonPositivityDelta:combine(1.0, fIonPositivityDrag, -1.0, distfIon)
  -- Use {ionEnergyAtCellBeforePositivity, ionEnergyAtCellAfterPositivity, ionEnergyAtCellAfterMaxDrag}
  -- to get correct energy at each node
  runUpdater(positivityScaleIonUpdater, tCurr, myDt, {fIonPositivityDelta, ionEnergyAtCellBeforePositivity,
    ionEnergyAtCellAfterPositivity, ionEnergyAtCellAfterMaxDrag}, {distfIon})
  distfIon:sync()

  -- compute elc parallel energy at nodes after positivity is applied
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilParaElc},
    {elcEnergyAtCellAfterPositivity})
  -- need to sync ghost cells for drag update
  distfElc:sync()
  -- take a maximum drag step
  runUpdater(positivityDragElcUpdater, tCurr, myDt, {distfElc, bField3d, elcEnergyAtCellAfterPositivity,
    elcEnergyAtCellBeforePositivity}, {fElcPositivityDrag})
  -- compute energy in cell after maximum drag step
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcPositivityDrag, jacobianFieldElc, hamilParaElc},
    {elcEnergyAtCellAfterMaxDrag, elcEnergyHistory})
  -- Compute the difference between distf post-positivity and distf post-drag
  fElcPositivityDelta:combine(1.0, fElcPositivityDrag, -1.0, distfElc)
  -- Use {elcEnergyAtCellBeforePositivity, elcEnergyAtCellAfterPositivity, elcEnergyAtCellAfterMaxDrag}
  -- to get correct energy at each node
  runUpdater(positivityScaleElcUpdater, tCurr, myDt, {fElcPositivityDelta, elcEnergyAtCellBeforePositivity,
    elcEnergyAtCellAfterPositivity, elcEnergyAtCellAfterMaxDrag}, {distfElc})
  distfElc:sync()

  return true
end

function applyPositivityDiagnostics(tCurr, myDt, distfElc, distfIon)
  -- calculate initial ion energy
  runUpdater(ionEnergyAtCellCalc, tCurr, myDt, {distfIon, jacobianFieldIon, hamilParaIon},
    {ionEnergyAtCellBeforePositivity, ionEnergyHistory})
  runUpdater(numDensityIonCalc, tCurr, myDt, {distfIon, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- store ion density temporarily in avgIonDensity structure
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityIon}, {avgIonDensity})
  local ionEnergyStart = ionEnergyHistory:lastInsertedData()
  local ionNumberStart = avgIonDensity:lastInsertedData()
  -- calculate initial elc energy
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilParaElc},
    {elcEnergyAtCellBeforePositivity, elcEnergyHistory})
  runUpdater(numDensityElcCalc, tCurr, myDt, {distfElc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  -- store elc density temporarily in avgElcDensity structure
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {avgElcDensity})
  local elcEnergyStart = elcEnergyHistory:lastInsertedData()
  local elcNumberStart = avgElcDensity:lastInsertedData()

  -- apply positivity to ions
  --Lucee.logInfo(string.format("-- ion positivity start --"))
  local ionStatus = runUpdater(positivityIonUpdater, tCurr, myDt, {bField3d}, {distfIon})
  
  -- apply positivity to elcs
  --Lucee.logInfo(string.format("-- elc positivity start --"))
  local elcStatus = runUpdater(positivityElcUpdater, tCurr, myDt, {bField3d}, {distfElc})

  if (elcStatus == false or ionStatus == false) then
    Lucee.logInfo(string.format("-- ending execution because an entire cell was negative --"))
    endExecution = true
    return false
  end
  
  -- compute ion parallel energy at nodes after positivity is applied
  runUpdater(ionEnergyAtCellCalc, tCurr, myDt, {distfIon, jacobianFieldIon, hamilParaIon},
    {ionEnergyAtCellAfterPositivity})

  -- temporary diagnostic
  runUpdater(numDensityIonCalc, tCurr, myDt, {distfIon, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- store ion density temporarily in avgIonDensity structure
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityIon}, {avgIonDensity})
  local ionEnergyPos = ionEnergyHistory:lastInsertedData()
  local ionNumberPos = avgIonDensity:lastInsertedData()

  -- need to sync ghost cells for drag update
  distfIon:sync()
  runUpdater(positivityDragIonUpdater, tCurr, myDt, {distfIon, bField3d, ionEnergyAtCellAfterPositivity,
    ionEnergyAtCellBeforePositivity}, {fIonPositivityDrag})
  -- compute ion parallel energy at nodes after maximum drag step
  runUpdater(ionEnergyAtCellCalc, tCurr, myDt, {fIonPositivityDrag, jacobianFieldIon, hamilParaIon},
    {ionEnergyAtCellAfterMaxDrag})

  -- temporary diagnostic
  runUpdater(numDensityIonCalc, tCurr, myDt, {fIonPositivityDrag, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- store ion density temporarily in avgIonDensity structure
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityIon}, {avgIonDensity})
  local ionEnergyDrag = ionEnergyHistory:lastInsertedData()
  local ionNumberDrag = avgIonDensity:lastInsertedData()

  -- Compute the difference between distf post-positivity and distf post-drag
  fIonPositivityDelta:combine(1.0, fIonPositivityDrag, -1.0, distfIon)
  -- Use {ionEnergyAtCellBeforePositivity, ionEnergyAtCellAfterPositivity, ionEnergyAtCellAfterMaxDrag}
  -- to get correct energy at each node
  runUpdater(positivityScaleIonUpdater, tCurr, myDt, {fIonPositivityDelta, ionEnergyAtCellBeforePositivity,
    ionEnergyAtCellAfterPositivity, ionEnergyAtCellAfterMaxDrag}, {distfIon})
  distfIon:sync()

  -- compute elc parallel energy at nodes after positivity is applied
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilParaElc},
    {elcEnergyAtCellAfterPositivity})

  -- temporary diagnostic
  runUpdater(numDensityElcCalc, tCurr, myDt, {distfElc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  -- store elc density temporarily in avgElcDensity structure
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {avgElcDensity})
  local elcEnergyPos = elcEnergyHistory:lastInsertedData()
  local elcNumberPos = avgElcDensity:lastInsertedData()

  -- need to sync ghost cells for drag update
  distfElc:sync()
  -- take a maximum drag step
  runUpdater(positivityDragElcUpdater, tCurr, myDt, {distfElc, bField3d, elcEnergyAtCellAfterPositivity,
    elcEnergyAtCellBeforePositivity}, {fElcPositivityDrag})
  -- compute energy in cell after maximum drag step
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcPositivityDrag, jacobianFieldElc, hamilParaElc},
    {elcEnergyAtCellAfterMaxDrag, elcEnergyHistory})

  -- temporary diagnostic
  runUpdater(numDensityElcCalc, tCurr, myDt, {fElcPositivityDrag, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  -- store elc density temporarily in avgElcDensity structure
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {avgElcDensity})
  local elcEnergyDrag = elcEnergyHistory:lastInsertedData()
  local elcNumberDrag = avgElcDensity:lastInsertedData()

  -- Compute the difference between distf post-positivity and distf post-drag
  fElcPositivityDelta:combine(1.0, fElcPositivityDrag, -1.0, distfElc)
  -- Use {elcEnergyAtCellBeforePositivity, elcEnergyAtCellAfterPositivity, elcEnergyAtCellAfterMaxDrag}
  -- to get correct energy at each node
  runUpdater(positivityScaleElcUpdater, tCurr, myDt, {fElcPositivityDelta, elcEnergyAtCellBeforePositivity,
    elcEnergyAtCellAfterPositivity, elcEnergyAtCellAfterMaxDrag}, {distfElc})
  distfElc:sync()

  -- calculate final ion energy
  runUpdater(ionEnergyAtCellCalc, tCurr, myDt, {distfIon, jacobianFieldIon, hamilParaIon},
    {ionEnergyAtCellAfterPositivity, ionEnergyHistory})
  runUpdater(numDensityIonCalc, tCurr, myDt, {distfIon, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- store ion density temporarily in avgIonDensity structure
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityIon}, {avgIonDensity})
  local ionEnergyEnd = ionEnergyHistory:lastInsertedData()
  local ionNumberEnd = avgIonDensity:lastInsertedData()
  -- calculate final elc energy
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilParaElc},
    {elcEnergyAtCellAfterPositivity, elcEnergyHistory})
  runUpdater(numDensityElcCalc, tCurr, myDt, {distfElc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  -- store elc density temporarily in avgElcDensity structure
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {avgElcDensity})
  local elcEnergyEnd = elcEnergyHistory:lastInsertedData()
  local elcNumberEnd = avgElcDensity:lastInsertedData()

  Lucee.logInfo(string.format("-- ion diagnostics --"))
  Lucee.logInfo(string.format("-- start N = %.10g --", ionNumberStart))
  Lucee.logInfo(string.format("-- pos N = %.10g --", ionNumberPos))
  Lucee.logInfo(string.format("-- drag N = %.10g --", ionNumberDrag))
  Lucee.logInfo(string.format("-- end N = %.10g --", ionNumberEnd))
  Lucee.logInfo(string.format("-- start E = %.10g --", ionEnergyStart))
  Lucee.logInfo(string.format("-- pos E = %.10g --", ionEnergyPos))
  Lucee.logInfo(string.format("-- drag E = %.10g --", ionEnergyDrag))
  Lucee.logInfo(string.format("-- end E = %.10g --", ionEnergyEnd))

  Lucee.logInfo(string.format("-- elc diagnostics --"))
  Lucee.logInfo(string.format("-- start N = %.10g --", elcNumberStart))
  Lucee.logInfo(string.format("-- pos N = %.10g --", elcNumberPos))
  Lucee.logInfo(string.format("-- drag N = %.10g --", elcNumberDrag))
  Lucee.logInfo(string.format("-- end N = %.10g --", elcNumberEnd))
  Lucee.logInfo(string.format("-- start E = %.10g --", elcEnergyStart))
  Lucee.logInfo(string.format("-- pos E = %.10g --", elcEnergyPos))
  Lucee.logInfo(string.format("-- drag E = %.10g --", elcEnergyDrag))
  Lucee.logInfo(string.format("-- end E = %.10g --", elcEnergyEnd))

  -- check that we haven't modified the total number of ions
  if math.abs((ionNumberEnd-ionNumberStart)/ionNumberStart) > 1e-10 or
    math.abs((ionEnergyEnd-ionEnergyStart)/ionEnergyStart) > 1e-10 then
    Lucee.logInfo(string.format("-- positivity modified ion number or energy --"))
  end

  -- check that we haven't modified the total number of electrons
  if math.abs((elcNumberEnd-elcNumberStart)/elcNumberStart) > 1e-10 or
    math.abs((elcEnergyEnd-elcEnergyStart)/elcEnergyStart) > 1e-10 then
    Lucee.logInfo(string.format("-- positivity modified elc number or energy --"))
  end

  return true
end

endExecution = false

function calcCollisionsIon(tCurr, myDt, fInitial, fFinal)
  -- fFinal f(t_{n+1}) is what to add increment to
  -- fInitial is f(t_n)

  -- compute number density, temperature, and mom1dir3 of input
  runUpdater(generalIntegralIonCalc, 0.0, 0.0, {fInitial, oneFieldIon, bField3d}, {numDensityIon})
  runUpdater(generalIntegralIonCalc, 0.0, 0.0, {fInitial, hamilDerivKeIon, bField3d}, {mom1dir3Ion})
  runUpdater(generalIntegralIonCalc, 0.0, 0.0, {fInitial, hamilKeIon, bField3d}, {kineticEnergyIon})
  local negativeTemperature = runUpdater(temperatureIonCalc, 0.0, 0.0, {numDensityIon,kineticEnergyIon,mom1dir3Ion,bField3d},
    {temperatureIon})

  if (negativeTemperature == false) then
    endExecution = true
    numDensityIon:write( string.format("nIonBreak_%d.h5", 0), 0.0)
    Lucee.logInfo(string.format("-- Ending execution because negative temperature. --"))
    return
  end
  
  mom1dir3IonMaxwellianInput:copy(mom1dir3Ion)
  temperatureIonMaxwellianInput:copy(temperatureIon)

  local convergenceStatus = false
  local iterCount = 0
  local maxIter = 100

  local mom1dir3Nan = false
  local temperatureNan = false
  -- loop to find maxwellian with desired density, <v>, and temperature
  while convergenceStatus == false and iterCount < maxIter do

    --Lucee.logInfo (string.format("-- Iteration %d", iterCount))

    runUpdater(maxwellianEval, 0.0, 0.0, {numDensityIon, mom1dir3IonMaxwellianInput,
      temperatureIonMaxwellianInput, bField3d}, {fIonMaxwellian})
    fIonMaxwellian:sync()
    -- compute parameters of maxwellian (it should already have correct number density)
    runUpdater(generalIntegralIonCalc, 0.0, 0.0, {fIonMaxwellian, hamilDerivKeIon, bField3d}, {mom1dir3IonMaxwellianNumerical})
    runUpdater(generalIntegralIonCalc, 0.0, 0.0, {fIonMaxwellian, hamilKeIon, bField3d}, {kineticEnergyIon})
    -- compute ion temperature
    runUpdater(temperatureIonCalc, 0.0, 0.0, {numDensityIon,kineticEnergyIon,mom1dir3IonNumerical,bField3d},
      {temperatureIonMaxwellianNumerical})

    -- compute input parameters for next iteration
    convergenceStatus = runUpdater(maxwellianParameterFinder, 0.0, 0.0, {mom1dir3Ion, temperatureIon,
      mom1dir3IonMaxwellianInput, temperatureIonMaxwellianInput,
      mom1dir3IonMaxwellianNumerical, temperatureIonMaxwellianNumerical}, {mom1dir3IonMaxwellianOutput,
      temperatureIonMaxwellianOutput})

    mom1dir3IonMaxwellianInput:copy(mom1dir3IonMaxwellianOutput)
    temperatureIonMaxwellianInput:copy(temperatureIonMaxwellianOutput)
    iterCount = iterCount + 1
  end

  if convergenceStatus == false then
    Lucee.logInfo (string.format("-- Unable to find Maxwellian after %d iterations.", iterCount))
  end

  -- calculate total electron density
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityIon}, {avgIonDensity})
  runUpdater(bgkCollisionOperator, 0.0, 0.0, {fInitial, fIonMaxwellian, numDensityIon,
    temperatureIonMaxwellianNumerical}, {fIonCollisions})
  fFinal:accumulate(myDt, fIonCollisions)
end

endExecution = false
-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
  -- RK stage 1
  f1Ion:copy(fIon)
  local myStatusCollisionIon = calcCollisionsIon(tCurr, myDt, fIon, f1Ion)

  if (myStatusCollisionIon == false or endExecution == true) then
    return false, myDt
  end

  fIon:copy(f1Ion)

  calcContinuousDiagnostics(tCurr, myDt)
  
  return true, myDt
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
        calcDiagnostics(tCurr, myDt)
        tCurr = tCurr + myDt
        Lucee.logInfo (string.format("** Something failed, ending execution "))
        break
      else
        -- time-step too large
        myDt = dtSuggested
        Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
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
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilKeElc, bField3d}, {kineticEnergyElc})
  -- compute electron temperature
  runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc,kineticEnergyElc,mom1dir3Elc,bField3d},
    {temperatureElc})
  -- convert to ev
  temperatureElc:scale(1/eV)

  -- compute ion moments
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, oneFieldIon, bField3d}, {numDensityIonAlt})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, hamilDerivKeIon, bField3d}, {mom1dir3Ion})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, hamilKeIon, bField3d}, {kineticEnergyIon})
  -- compute ion temperature
  runUpdater(temperatureIonCalc, tCurr, myDt, {numDensityIonAlt,kineticEnergyIon,mom1dir3Ion,bField3d},
    {temperatureIon})
  -- convert to ev
  temperatureIon:scale(1/eV)

  -- compute average ion temperature (volume factor should be removed in post-processing)
  runUpdater(fieldInt3d, tCurr, myDt, {temperatureElc}, {elcTempHistory})
  runUpdater(fieldInt3d, tCurr, myDt, {temperatureIon}, {ionTempHistory})
end

-- function containing diagnostics that must be called at the end of every time step
function calcContinuousDiagnostics(tCurr, myDt)
  -- store phi in a cell
  runUpdater(recordFieldInCellPhiCalc, tCurr, myDt, {mom1dir3Ion}, {fieldInCellPhi})
  -- compute total guiding center density by integrating over whole domain
  -- note: includes a volume factor that should be removed in post-processing
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {totalElcCounter})
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityIon}, {totalIonCounter})
  -- total energy
  runUpdater(ionEnergyAtCellCalc, tCurr, myDt, {fIon, jacobianFieldIon, hamilKeIon},
    {ionEnergyAtCellBeforePositivity, ionEnergyHistoryWrite})

  return true
end

-- write data to H5 file
function writeFields(frameNum, tCurr)
  numDensityElc:write( string.format("nElc_%d.h5", frameNum), tCurr)
  temperatureElc:write( string.format("tElc_%d.h5", frameNum), tCurr)
  numDensityIon:write( string.format("nIon_%d.h5", frameNum), tCurr)
  
  temperatureIon:write( string.format("tIon_%d.h5", frameNum), tCurr)
  --numDensityDelta:write( string.format("nDelta_%d.h5", frameNum), tCurr)

  fElc:write( string.format("fElc_%d.h5", frameNum), tCurr)
  fIon:write( string.format("fIon_%d.h5", frameNum), tCurr)
  
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

  fieldInCellPhi:write( string.format("phiInCell_%d.h5", frameNum), tCurr)

  ionEnergyHistoryWrite:write( string.format("ionEnergyHistory_%d.h5", frameNum), tCurr)
end


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
--runUpdater(multiply5dCalcIon, 0.0, 0.0, {fIon, fIonInitialPerturb}, {fDupIon})
--fIon:copy(fDupIon)
--fIon:sync()
--runUpdater(multiply5dCalcElc, 0.0, 0.0, {fElc, fElcInitialPerturb}, {fDupElc})
--fElc:copy(fDupElc)
--fElc:sync()
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

calcDiagnostics(tCurr, 0.0)
calcContinuousDiagnostics(tCurr, 0.0)
writeFields(startFrame-1,tCurr)

maxwellianEval = Updater.SOLMaxwellianAtNodeCalc {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/ionMass,
  speciesMass = ionMass,
}
mom1dir3IonMaxwellianInput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureIonMaxwellianInput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
mom1dir3IonMaxwellianNumerical = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureIonMaxwellianNumerical = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
mom1dir3IonMaxwellianOutput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
temperatureIonMaxwellianOutput = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

maxwellianParameterFinder = Updater.SOLMaxwellianParameterCalc {
  onGrid = grid_3d,
  basis3d = basis_3d,
}

-- distribution function for ions
fIonCollision = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
-- apply BGK operator
bgkCollisionOperator = Updater.SOLBGKCollisionUpdater {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  alpha = function(t)
    return getIonAlpha(t)
  end,
  onlyIncrement = true,
}

for frame = startFrame, nFrames do
  if endExecution == false then
    dtSuggested = 10e-6
    Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
    dtSuggested, tCurr = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
    writeFields(frame, tCurr)
    Lucee.logInfo ("")
  else
    Lucee.logInfo (string.format("-- Not advancing solution from %g", tCurr))
  end
end
