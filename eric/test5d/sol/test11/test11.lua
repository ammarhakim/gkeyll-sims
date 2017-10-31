-- Input file for a 3X2V SOL Simulation
-- 01-28-2016 Creation
-- 03-15-2016 First attempt at input file with kinetic elections and ions
-- Command line parameters
-- mpiexec -n 4 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i test11.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist
-- added positivity back in
-- number of parallel velocity space cells changed to 6, max v = 6*vt
-- changed boundary condition on lhs to be the same as that on rhs

-- phase-space decomposition (5d)
phaseDecomp = DecompRegionCalc5D.CartProd { cuts = {2, 2, 1, 1, 1} }
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
tEnd = 3000e-6
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 3000
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

-- physical parameters
eV        = Lucee.ElementaryCharge
elcCharge = -eV
ionCharge = eV
elcMass   = Lucee.ElectronMass
ionMass   = 2.014*Lucee.ProtonMass 
elcTemp   = 3 -- [eV]
ionTemp   = 1 -- [eV]
B0 = 0.076  -- [T]
n0 = 10^16 -- [1/m^3]
N = 4 -- field line turns
kPerpTimesRhoS = 0.2
sourceTemp    = 3.5 -- [eV]
-- derived parameters
vtElc   = math.sqrt(sourceTemp*eV/elcMass)
omega_e = math.abs(elcCharge*B0/elcMass)
rho_e   = vtElc/omega_e

vtIon   = math.sqrt(sourceTemp*eV/ionMass)
omega_i = math.abs(ionCharge*B0/ionMass)
rho_i   = vtIon/omega_i

rho_s = math.sqrt(sourceTemp*eV/ionMass)/omega_i
deltaX  = 64*rho_s
deltaY  = 64*rho_s/N
R          = 200*rho_s -- [m]
L_parallel = 2*math.pi*R*N -- [m]
-- grid parameters: number of cells
N_X = 8
N_Y = 8
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
sourceLambda    = 20*rho_s
sourceX         = X_LOWER + 32*rho_s
sourceAmplitude = 1.5e20 -- [1/(m^3*sec)] 
-- total volume
totalVol = (X_UPPER-X_LOWER)*(Y_UPPER-Y_LOWER)*(Z_UPPER-Z_LOWER)

VPARA_UPPER_ELC = 6*vtElc
VPARA_LOWER_ELC = -VPARA_UPPER_ELC
MU_LOWER_ELC = 0
MU_UPPER_ELC = elcMass*(VPARA_UPPER_ELC*VPARA_UPPER_ELC)/(2*B0)

VPARA_UPPER_ION = 6*vtIon
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
hamilDupIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
-- Magnetic field
function bFieldProfile(x)
  return B0*R/x
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

function elcTempProfile(x)
  return elcTemp--+ elcTemp*(1 - (x-X_LOWER)/(X_UPPER-X_LOWER))
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
		 return sourceAmplitude*math.exp(-(x-sourceX)^2/sourceLambda^2)*
		  (2*math.pi*sourceTemp*eV/ionMass)^(-3/2)*
      math.exp(-ionMass*v^2/(2*sourceTemp*eV))*
      math.exp(-math.abs(mu)*bFieldProfile(x)/(sourceTemp*eV))
	 end
}
runUpdater(initSourceIon, 0.0, 0.0, {}, {ionSource})

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
		 return sourceAmplitude*math.exp(-(x-sourceX)^2/sourceLambda^2)*
		  (2*math.pi*sourceTemp*eV/elcMass)^(-3/2)*
      math.exp(-elcMass*v^2/(2*sourceTemp*eV))*
      math.exp(-math.abs(mu)*bFieldProfile(x)/(sourceTemp*eV))
	 end
}
runUpdater(initSourceElc, 0.0, 0.0, {}, {elcSource})

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
  integrateGhosts = false, -- can be false, but use true to make sure fIon = 0 in ghosts
}
-- to compute electron flux at domain boundaries in z
edgeMomentsElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
edgeMomentCalcElc = Updater.MomentsAtEdges5D {
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  onGrid = grid_elc_5d,
  polyOrder = polyOrder,
  integrateGhosts = true, -- only want outward flux
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
-- to compute sheath potential on lower and upper surfaces in z
logicalSheathCalc = Updater.LogicalSheath5D {
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  onGrid = grid_elc_5d,
  polyOrder = polyOrder,
  ionMass = ionMass,
  elcMass = elcMass,
  eV = eV,
}
-- Use to verify that sheath bc has been properly applied to electrons
-- Should be used only for debugging
edgeMomentCalcElcPostBC = Updater.MomentsAtEdges5D {
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  onGrid = grid_elc_5d,
  polyOrder = polyOrder,
  integrateGhosts = true,
}
-- Sets potential on upper x surface to a constant value
-- Also averages solution on boundaries to make sure potential smoothing operation
-- has smooth boundary conditions
upperXPotentialBc = Updater.SOLUpperXPotentialBcUpdater {
  onGrid = grid_2d,
  basis = basis_2d,
  applyLowerEdge = true,
}
-- Stores value of phi_s(x = x_uppper)
phiSAtEdge = DataStruct.DynVector { numComponents = 2, }
-- Stores smoothed sheath potential on lower surface
phiSLowerSmoothed = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- Stores smoothed sheath potential on upper surface
phiSUpperSmoothed = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- Used to compute a continuous potential field
fieldSmoother = Updater.FemPoisson2D {
   onGrid = grid_2d,
   basis = basis_2d,
   --periodicDirs = {1},
   sourceNodesShared = false, -- default true
   solutionNodesShared = false, -- default true
   writeStiffnessMatrix = false,
   modifierConstant = 1.0,
   laplacianWeight = 0.0, -- indicates smoothing
   -- use variable dirichlet bc in all directions
   bcLeft = { T ="D_VAR", V = 0.0},
   bcRight = { T ="D_VAR", V = 0.0},
   bcBottom = { T ="D_VAR", V = 0.0},
   bcTop = { T ="D_VAR", V = 0.0},
}
-- Used to store boundary conditions for potential solve
dirichletBoundaryField = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- Used to fill out boundary condition field for potential solve
dirichletFieldAssembler = Updater.SOLSetPotentialAtBoundary {
  onGrid = grid_3d,
  basis = basis_3d,
  applyLowerEdge = true,
}
-- stores result of potential solve
phi3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- solves gyrokinetic poission equation
fieldSolver = Updater.FemPoisson3D {
   onGrid = grid_3d,
   basis = basis_3d,
   periodicDirs = {1},
   sourceNodesShared = false, -- default true
   solutionNodesShared = false, -- default true
   writeStiffnessMatrix = false,
   modifierConstant = 0,--kPerpTimesRhoS^2,
   laplacianWeight = rho_s^2,
   isGyroKineticPoisson = true,
   -- boundary conditions to apply
   bcLeft = { T ="D_VAR", V = 0.0},
   bcRight = { T ="D_VAR", V = 0.0},
   --bcBottom = { T ="D_VAR", V = 0.0},
   --bcTop = { T ="D_VAR", V = 0.0},
   bcBack = { T ="D_VAR", V = 0.0},
   bcFront = { T ="D_VAR", V = 0.0},
}
-- gyrokinetic equation for ions
gyroEqnIon = PoissonBracketEquation.GyroEquation5D {
  speciesMass = ionMass,
  speciesCharge = ionCharge,
  bStarY = bStarYFieldIon,
  bStarZ = bStarZFieldIon,
}
-- gyrokinetic equation for electrons
gyroEqnElc = PoissonBracketEquation.GyroEquation5D {
  speciesMass = elcMass,
  speciesCharge = elcCharge,
  bStarY = bStarYFieldElc,
  bStarZ = bStarZFieldElc,
}
pbSlvrIon = Updater.PoissonBracketImp5D {
   onGrid = grid_ion_5d,
   basis = basis_ion_5d,
   cfl = cfl,
   equation = gyroEqnIon,
   -- let solver know about additional jacobian factor
   jacobianField = jacobianFieldIon,
   updateDirections = {0,1,2,3},
   zeroFluxDirectionsLower = {0,3,4},
   zeroFluxDirectionsUpper = {0,3,4},
   fluxType = "upwind",
}
pbSlvrElc = Updater.PoissonBracketImp5D {
   onGrid = grid_elc_5d,
   basis = basis_elc_5d,
   cfl = cfl,
   equation = gyroEqnElc,
   -- let solver know about additional jacobian factor
   jacobianField = jacobianFieldElc,
   updateDirections = {0,1,2,3},
   zeroFluxDirectionsLower = {0,3,4},
   zeroFluxDirectionsUpper = {0,3,4},
   fluxType = "upwind",
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
temperatureElcCalc = Updater.SOLTemperatureCalc {
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
temperatureIonCalc = Updater.SOLTemperatureCalc {
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

function applyPositivity(distfElc, distfIon)
  runUpdater(positivityElcUpdater, 0.0, 0.0, {bField3d}, {distfElc})
  runUpdater(positivityIonUpdater, 0.0, 0.0, {bField3d}, {distfIon})
end

function computePotential(phiSLower, phiSUpper, nIon, nElc, phiOut)
  -- average lower and upper sheath potential fields
  phiSAverage:combine(0.5, phiSLower, 0.5, phiSUpper)
  phiSAverage:sync()
  -- smooth sheath potential fields at domain boundaries and replace
  -- potential with an averaged potential on the upper-x boundary
  runUpdater(upperXPotentialBc, 0.0, 0.0, {phiSAverage}, {phiSLower, phiSUpper, phiSAtEdge})
  phiSLower:sync()
  phiSUpper:sync()
  -- make discontinuous potentials continuous via non-local solve
  local statusPhiLower = runUpdater(fieldSmoother, 0.0, 0.0, {phiSLower, phiSLower}, {phiSLowerSmoothed})
  local statusPhiUpper = runUpdater(fieldSmoother, 0.0, 0.0, {phiSUpper, phiSUpper}, {phiSUpperSmoothed})
  if statusPhiLower == false or statusPhiUpper == false then
    Lucee.logInfo(string.format("-- phi lower or upper solve failed --"))
  end
  phiSLowerSmoothed:sync()
  phiSUpperSmoothed:sync()
  -- combine sheath potential data into a single field (no sync required)
  runUpdater(dirichletFieldAssembler, 0.0, 0.0, {phiSLowerSmoothed, phiSUpperSmoothed, phiSAtEdge},
    {dirichletBoundaryField})
  dirichletBoundaryField:sync()
  -- calculate average ion density
  runUpdater(fieldInt3d, 0.0, 0.0, {nIon}, {avgIonDensity})
  local n_i0 = avgIonDensity:lastInsertedData()/totalVol
  -- solve for potential
  numDensityDelta:combine(1.0, nElc, -1.0, nIon)
  numDensityDelta:scale(elcTemp/(n_i0))

  local statusPhi = runUpdater(fieldSolver, 0.0, 0.0, {numDensityDelta, dirichletBoundaryField}, {phiOut})
  if statusPhi == false then
    Lucee.logInfo(string.format("-- phi3d solve failed --"))
  end
end

function computeHamiltonianElc(phiIn, hamilOut)
  -- copy phi3d into a 5d field
  runUpdater(copy3dTo5dElc, 0.0, 0.0, {phiIn}, {hamilOut})
  hamilOut:scale(elcCharge)
  -- add the kinetic energy component
  hamilOut:accumulate(1.0, hamilKeElc)
  hamilOut:sync()
end

function computeHamiltonianIon(phiIn, hamilOut)
  -- copy phi3d into a 5d field
  runUpdater(copy3dTo5dIon, 0.0, 0.0, {phiIn}, {hamilOut})
  hamilOut:scale(ionCharge)
  -- add the kinetic energy component
  hamilOut:accumulate(1.0, hamilKeIon)
  hamilOut:sync()
end

function applyBcToDistF(tCurr, myDt, fElc, fIon)
  -- apply periodic bc in y to ions
  fIon:sync()
  -- apply periodic bc in y to electrons
  fElc:sync()
  -- apply zero inflow boundary conditions in z to ions
  runUpdater(bcLowerZIon, tCurr, myDt, {}, {fIon})
  runUpdater(bcUpperZIon, tCurr, myDt, {}, {fIon})
  -- calculate first vPar moment of the ion distribution function
  runUpdater(edgeMomentCalcIon, tCurr, myDt, {fIon}, {edgeMomentsIon})
  -- apply logical sheath boundary conditions to electron distribution function
  sheathStatus = runUpdater(logicalSheathCalc, tCurr, myDt, {edgeMomentsIon}, {fElc, phiSLower, phiSUpper})
  if sheathStatus == false then
    Lucee.logInfo(string.format("-- cutoff velocity solve failed --"))
  end
  -- sync cells in distribution function
  fElc:sync()
  -- sync ghost cells in phiS
  phiSLower:sync()
  phiSUpper:sync()

  return sheathStatus
end

endExecution = false
-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
  -- RK stage 1
  local myStatusElc, myDtSuggestedElc = runUpdater(pbSlvrElc, tCurr, myDt, {fElc, hamilElc}, {f1Elc})
  local myStatusIon, myDtSuggestedIon = runUpdater(pbSlvrIon, tCurr, myDt, {fIon, hamilIon}, {f1Ion})
  if (myStatusElc == false or myStatusIon == false) then
    return false, math.min(myDtSuggestedElc, myDtSuggestedIon)
  end
  -- add on particle source
  f1Elc:accumulate(myDt, elcSource)
  f1Ion:accumulate(myDt, ionSource)
  applyPositivity(f1Elc, f1Ion)
  -- apply boundary conditions to distribution functions
  local bcStatus = applyBcToDistF(tCurr, myDt, f1Elc, f1Ion)
  if (bcStatus == false) then
    endExecution = true
    fElc:copy(f1Elc)
    fIon:copy(f1Ion)
    return false, math.min(myDtSuggestedElc, myDtSuggestedIon)
  end
  -- calculate density of each species for potential solve
  runUpdater(numDensityElcCalc, tCurr, myDt, {f1Elc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  runUpdater(numDensityIonCalc, tCurr, myDt, {f1Ion, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- solve for potential
  computePotential(phiSLower, phiSUpper, numDensityIon, numDensityElc, phi3d)
  -- calculate new hamiltonians
  computeHamiltonianElc(phi3d, hamilElc)
  computeHamiltonianIon(phi3d, hamilIon)

  -- RK stage 2
  myStatusElc, myDtSuggestedElc = runUpdater(pbSlvrElc, tCurr, myDt, {f1Elc, hamilElc}, {fNewElc})
  myStatusIon, myDtSuggestedIon = runUpdater(pbSlvrIon, tCurr, myDt, {f1Ion, hamilIon}, {fNewIon})
  if (myStatusElc == false or myStatusIon == false) then
    return false, math.min(myDtSuggestedElc, myDtSuggestedIon)
  end
  -- add on particle source
  fNewElc:accumulate(myDt, elcSource)
  fNewIon:accumulate(myDt, ionSource)
  applyPositivity(fNewElc, fNewIon)
  f1Elc:combine(3.0/4.0, fElc, 1.0/4.0, fNewElc)
  f1Ion:combine(3.0/4.0, fIon, 1.0/4.0, fNewIon)
  -- apply boundary conditions to distribution functions
  bcStatus = applyBcToDistF(tCurr, myDt, f1Elc, f1Ion)
  if (bcStatus == false) then
    endExecution = true
    fElc:copy(f1Elc)
    fIon:copy(f1Ion)
    return false, math.min(myDtSuggestedElc, myDtSuggestedIon)
  end
  -- calculate density of each species for potential solve
  runUpdater(numDensityElcCalc, tCurr, myDt, {f1Elc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  runUpdater(numDensityIonCalc, tCurr, myDt, {f1Ion, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- solve for potential
  computePotential(phiSLower, phiSUpper, numDensityIon, numDensityElc, phi3d)
  -- calculate new hamiltonians
  computeHamiltonianElc(phi3d, hamilElc)
  computeHamiltonianIon(phi3d, hamilIon)

  -- RK stage 3
  myStatusElc, myDtSuggestedElc = runUpdater(pbSlvrElc, tCurr, myDt, {f1Elc, hamilElc}, {fNewElc})
  myStatusIon, myDtSuggestedIon = runUpdater(pbSlvrIon, tCurr, myDt, {f1Ion, hamilIon}, {fNewIon})
  if (myStatusElc == false or myStatusIon == false) then
    return false, math.min(myDtSuggestedElc, myDtSuggestedIon)
  end
  -- add on particle source
  fNewElc:accumulate(myDt, elcSource)
  fNewIon:accumulate(myDt, ionSource)
  applyPositivity(fNewElc, fNewIon)
  f1Elc:combine(1.0/3.0, fElc, 2.0/3.0, fNewElc)
  f1Ion:combine(1.0/3.0, fIon, 2.0/3.0, fNewIon)
  -- apply boundary conditions to distribution functions
  bcStatus = applyBcToDistF(tCurr, myDt, f1Elc, f1Ion)
  if (bcStatus == false) then
    endExecution = true
    fElc:copy(f1Elc)
    fIon:copy(f1Ion)
    return false, math.min(myDtSuggestedElc, myDtSuggestedIon)
  end
  -- calculate density of each species for potential solve
  runUpdater(numDensityElcCalc, tCurr, myDt, {f1Elc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  runUpdater(numDensityIonCalc, tCurr, myDt, {f1Ion, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- solve for potential
  computePotential(phiSLower, phiSUpper, numDensityIon, numDensityElc, phi3d)
  -- calculate new hamiltonians
  computeHamiltonianElc(phi3d, hamilElc)
  computeHamiltonianIon(phi3d, hamilIon)

  fElc:copy(f1Elc)
  fIon:copy(f1Ion)

  return true, math.min(myDtSuggestedElc, myDtSuggestedIon)
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
        tCurr = tCurr + myDt
        Lucee.logInfo (string.format("** Boundary conditions failed, ending execution "))
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

function calcDiagnostics(tCurr, myDt)
  runUpdater(fieldInt2d, tCurr, myDt, {phiSLower}, {avgPhiSLower})
  -- compute total guiding center density by integrating over whole domain
  -- note: includes a volume factor that should be removed in post-processing
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {totalElcCounter})
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityIon}, {totalIonCounter})

  -- calculate first vPar moment of the electron distribution function
  runUpdater(edgeMomentCalcElc, tCurr, myDt, {fElc}, {edgeMomentsElc})

  -- compute electron moments
  runUpdater(mom1VParaElcCalc, tCurr, myDt, {fElc, bField3d}, {mom1dir3Elc})
  mom1dir3Elc:scale(2*math.pi/elcMass)
  runUpdater(mom2VParaElcCalc, tCurr, myDt, {fElc, bField3d}, {mom2dir3Elc})
  mom2dir3Elc:scale(2*math.pi/elcMass)
  runUpdater(mom1MuElcCalc, tCurr, myDt, {fElc, bField3d}, {mom1dir4Elc})
  mom1dir4Elc:scale(2*math.pi/elcMass)
  -- compute elctron temperature
  runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc,mom1dir3Elc,mom2dir3Elc,mom1dir4Elc,bField3d},
    {temperatureElc})
  -- convert to ev
  temperatureElc:scale(1/eV)

  -- compute ion moments
  runUpdater(mom1VParaIonCalc, tCurr, myDt, {fIon, bField3d}, {mom1dir3Ion})
  mom1dir3Ion:scale(2*math.pi/ionMass)
  runUpdater(mom2VParaIonCalc, tCurr, myDt, {fIon, bField3d}, {mom2dir3Ion})
  mom2dir3Ion:scale(2*math.pi/ionMass)
  runUpdater(mom1MuIonCalc, tCurr, myDt, {fIon, bField3d}, {mom1dir4Ion})
  mom1dir4Ion:scale(2*math.pi/ionMass)
  -- compute ion temperature
  runUpdater(temperatureIonCalc, tCurr, myDt, {numDensityIon,mom1dir3Ion,mom2dir3Ion,mom1dir4Ion,bField3d},
    {temperatureIon})
  -- convert to ev
  temperatureIon:scale(1/eV)
end

-- write data to H5 file
function writeFields(frameNum, tCurr)
  Lucee.logInfo( string.format("phiSE = %g\n",phiSAtEdge:lastInsertedData()) )
  numDensityElc:write( string.format("nElc_%d.h5", frameNum), tCurr)
  temperatureElc:write( string.format("tElc_%d.h5", frameNum), tCurr)
  numDensityIon:write( string.format("nIon_%d.h5", frameNum), tCurr)
  temperatureIon:write( string.format("tIon_%d.h5", frameNum), tCurr)
  --nInitialPerturbDG:write( string.format("nPerturb_%d.h5", frameNum), tCurr)
   --fieldEnergy:write( string.format("fieldEnergy_%d.h5", frameNum), tCurr)
  numDensityDelta:write( string.format("nDelta_%d.h5", frameNum), tCurr)
  
  phi3d:write( string.format("phi_%d.h5", frameNum), tCurr)
  phiSAverage:write( string.format("phiSA_%d.h5", frameNum), tCurr)
  phiSLower:write( string.format("phiSL_%d.h5", frameNum), tCurr)
  phiSLowerSmoothed:write( string.format("phiSLS_%d.h5", frameNum), tCurr)
  avgPhiSLower:write( string.format("avgPhiSL_%d.h5", frameNum), tCurr)
  phiSUpper:write( string.format("phiSU_%d.h5", frameNum), tCurr)

  dirichletBoundaryField:write( string.format("phiD_%d.h5", frameNum), tCurr)
  fElc:write( string.format("fElc_%d.h5", frameNum), tCurr)
  fIon:write( string.format("fIon_%d.h5", frameNum), tCurr)
  edgeMomentsIon:write( string.format("edgeMomentsIon_%d.h5", frameNum), tCurr)
  edgeMomentsElc:write( string.format("edgeMomentsElc_%d.h5", frameNum), tCurr)

  totalIonCounter:write( string.format("totalIonCounter_%d.h5", frameNum), tCurr)
  totalElcCounter:write( string.format("totalElcCounter_%d.h5", frameNum), tCurr)
end

-- ensure that source terms are in balance
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
print( string.format("ionSource scaled by %g (%g/%g)",totalSourceElc/totalSourceIon,totalSourceElc,totalSourceIon) )
ionSource:sync()
elcSource:sync()

if Lucee.IsRestarting then
  fElc:read("fElc_" .. Lucee.RestartFrame .. ".h5")
  fElc:sync()
  fIon:read("fIon_" .. Lucee.RestartFrame .. ".h5")
  fIon:sync()

  startFrame = Lucee.RestartFrame + 1
  tCurr = tStart + tFrame*Lucee.RestartFrame
  -- apply boundary conditions to distribution functions
  applyBcToDistF(0.0, 0.0, fElc, fIon)
  runUpdater(numDensityElcCalc, 0.0, 0.0, {fElc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  runUpdater(numDensityIonCalc, 0.0, 0.0, {fIon, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- solve for potential
  computePotential(phiSLower, phiSUpper, numDensityIon, numDensityElc, phi3d)
else
  -- Make sure distribution functions are still positive
  applyPositivity(fIon, fElc)
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
  print( string.format("fIon scaled by %g (%g/%g)",totalInitElc/totalInitIon,totalInitElc,totalInitIon) )
  fIon:sync()
  fElc:sync()
  tCurr = tStart
  startFrame = 1
  -- apply boundary conditions to distribution functions
  applyBcToDistF(0.0, 0.0, fElc, fIon)
  runUpdater(numDensityElcCalc, 0.0, 0.0, {fElc, bField3d}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  runUpdater(numDensityIonCalc, 0.0, 0.0, {fIon, bField3d}, {numDensityIon})
  numDensityIon:scale(2*math.pi/ionMass)
  -- solve for potential
  computePotential(phiSLower, phiSUpper, numDensityIon, numDensityElc, phi3d)
  -- zero out potential on first time step
  --phi3d:clear(0.0)
end

-- calculate new hamiltonians
computeHamiltonianElc(phi3d, hamilElc)
computeHamiltonianIon(phi3d, hamilIon)

-- Compute diagnostics for t = 0
calcDiagnostics(tCurr, 0.0)
writeFields(startFrame-1,tCurr)

for frame = startFrame, nFrames do
  if endExecution == false then
    Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
    dtSuggested, tCurr = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
    --tCurr = tCurr+tFrame
    writeFields(frame, tCurr)
    Lucee.logInfo ("")
  else
    Lucee.logInfo (string.format("-- Not advancing solution from %g", tCurr))
  end
end
