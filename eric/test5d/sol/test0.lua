-- Input file for a 3X2V SOL Simulation
-- 01-28-2016 Creation
-- Command line parameters
-- -pc_type lu -pc_factor_mat_solver_package superlu_dist

-- phase-space decomposition
phaseDecomp = DecompRegionCalc5D.CartProd { cuts = {1, 1, 1, 1, 1} }
-- configuration space decomposition
confDecomp = DecompRegionCalc3D.SubCartProd5D {
   decomposition = phaseDecomp,
   collectDirections = {0, 1, 2},
}
-- logical sheath decomposition
sheathDecomp = DecompRegionCalc2D.SubCartProd3D {
   decomposition = confDecomp,
   collectDirections = {0, 1},
}
--DecompRegionCalc2D.CartProd { cuts = {1, 1} }

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.2
-- parameters to control time-stepping
tStart = 0.0
tEnd = 1e-10
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 1
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

-- physical parameters
eV        = Lucee.ElementaryCharge
elcCharge = -eV
ionCharge = eV
elcMass   = Lucee.ElectronMass
ionMass   = 2.014*Lucee.ProtonMass 
elcTemp   = 100 -- [eV]
ionTemp   = 100 -- [eV]
B0 = 1  -- [T]
n0 = 10^19 -- [1/m^3]
-- derived parameters
R          = 1 -- [m]
L_parallel = 2 -- [m]
vtElc  = math.sqrt(elcTemp*eV/elcMass)
vtIon  = math.sqrt(ionTemp*eV/ionMass)
c_e        = math.sqrt(elcTemp*eV/elcMass)
omega_e    = math.abs(elcCharge*B0/elcMass)
rho_e      = c_e/omega_e

c_i        = math.sqrt(ionTemp*eV/ionMass)
omega_i    = math.abs(ionCharge*B0/ionMass)
rho_i      = c_i/omega_i

deltaR     = 32*rho_i
-- grid parameters: number of cells
N_X = 4
N_Y = 4
N_Z = 2
N_VPARA = 4
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = 0
X_UPPER = 2
Y_LOWER = -1
Y_UPPER = 1
Z_LOWER = 0
Z_UPPER = 2

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
}
-- distribution function for ions
fIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}

function bFieldProfile(x)
  return B0 --*R/x
end

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

-- write data to H5 file
function writeFields(frameNum, tCurr)
   --numDensityKinetic:write( string.format("n_%d.h5", frameNum), tCurr)
   --fieldEnergy:write( string.format("fieldEnergy_%d.h5", frameNum), tCurr)
   --numDensityDelta:write( string.format("nDelta_%d.h5", frameNum), tCurr)
   --phi3dSmoothed:write( string.format("phi_%d.h5", frameNum), tCurr)
   --vParaSquaredKinetic:write( string.format("vParaSq_%d.h5", frameNum), tCurr)
   --hamil:write( string.format("hamil_%d.h5", frameNum), tCurr)
   --f:write( string.format("f_%d.h5", frameNum), tCurr)
   --phi3d:write( string.format("phi_%d.h5", frameNum), tCurr)
   --f:write( string.format("f_%d.h5", frameNum), tCurr)
end

startFrame = 1

-- moment at edge testing
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
runUpdater(bcLowerZIon, 0.0, 0.0, {}, {fIon})
runUpdater(bcUpperZIon, 0.0, 0.0, {}, {fIon})
print('Ions with Ghosts Post-BC')
runUpdater(edgeMomentCalcIon, 0.0, 0.0, {fIon}, {edgeMomentsIon})
-- sheath potential
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
  integrateGhosts = false,
}
print('Pre-BC Electrons')
runUpdater(edgeMomentCalcElc, 0.0, 0.0, {fElc}, {edgeMomentsElc})

runUpdater(logicalSheathCalc, 0.0, 0.0, {edgeMomentsIon}, {fElc, phiSLower, phiSUpper})
--phiSUpper:write( string.format("phi_%d.h5", 0), 0)
-- calculate average phiS
phiSAverage = phiSLower:duplicate()
phiSAverage:copy(phiSLower)
phiSAverage:accumulate(1.0,phiSUpper)
phiSAverage:scale(0.5)
phiSAverage:write( string.format("phiA_%d.h5", 0), 0)

-- Verify that electrons have correct edge moments
edgeMomentCalcElcPostBC = Updater.MomentsAtEdges5D {
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  onGrid = grid_elc_5d,
  polyOrder = polyOrder,
  integrateGhosts = true,
}
print('Post-BC Electrons')
runUpdater(edgeMomentCalcElcPostBC, 0.0, 0.0, {fElc}, {edgeMomentsElc})

-- set random seed based on processor number to avoid imprinting
-- domain decomp on ICs.
r = Lucee.getRank()
math.randomseed(100000*r+os.time())
-- Test smoothing
initPhi = Updater.ProjectOnNodalBasis2D {
   onGrid = grid_2d,
   -- basis functions to use
   basis = basis_2d,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,t)
      return math.random()--.cos( (y-Y_LOWER)/(Y_UPPER-Y_LOWER)*math.pi )
	 end
}
runUpdater(initPhi, 0.0, 0.0, {}, {phiSLower})
runUpdater(initPhi, 0.0, 0.0, {}, {phiSUpper})

upperXPotentialBc = Updater.SOLUpperXPotentialBcUpdater {
  onGrid = grid_2d,
  basis = basis_2d,
}
phiSAverage:copy(phiSLower)
phiSAverage:accumulate(1.0, phiSUpper)
phiSAverage:scale(0.5)

phiSAtEdge = DataStruct.DynVector { numComponents = 1, }

--phiSLower:write( string.format("phi_%d.h5", 0), 0)
-- Modify potential on upper and lower z surfaces so that the value on the upper 0
-- boundary is a constant
phiSLower:sync()
phiSUpper:sync()
runUpdater(upperXPotentialBc, 0.0, 0.0, {phiSAverage}, {phiSLower, phiSUpper, phiSAtEdge})
phiSUpper:write( string.format("phi_%d.h5", 0), 0)
--print( string.format("phiAvg = %g", phiSAtEdge:lastInsertedData()) )

-- sheath potential (smoothed)
phiSLowerSmoothed = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
phiSUpperSmoothed = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
fieldSmoother = Updater.FemPoisson2D {
   onGrid = grid_2d,
   basis = basis_2d,
   periodicDirs = {1},
   sourceNodesShared = false, -- default true
   solutionNodesShared = false, -- default true
   writeStiffnessMatrix = false,
   modifierConstant = 1.0,
   laplacianWeight = 0.0, -- indicates smoothing
   -- boundary conditions to apply
   bcLeft = { T ="D_VAR", V = 0.0},
   bcRight = { T ="D_VAR", V = 0.0},
   --bcBottom = { T ="D_VAR", V = 0.0},
   --bcTop = { T ="D_VAR", V = 0.0},
}
print('hello')
runUpdater(fieldSmoother, 0.0, 0.0, {phiSLower, phiSLower}, {phiSLowerSmoothed})
runUpdater(fieldSmoother, 0.0, 0.0, {phiSUpper, phiSUpper}, {phiSUpperSmoothed})

--phiSLowerSmoothed:write( string.format("phiS_%d.h5", 0), 0)
--phiSUpperSmoothed:write( string.format("phiS_%d.h5", 0), 0)

-- Combine into a single 3d field
dirichletBoundaryField = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
dirichletFieldAssembler = Updater.SOLSetPotentialAtBoundary {
  onGrid = grid_3d,
  basis = basis_3d,
}
runUpdater(dirichletFieldAssembler, 0.0, 0.0, {phiSLowerSmoothed, phiSUpperSmoothed, phiSAtEdge},
{dirichletBoundaryField})

dirichletBoundaryField:write( string.format("phi3d_%d.h5", 0), 0)

fieldSolver = Updater.FemPoisson3D {
   onGrid = grid_3d,
   basis = basis_3d,
   --periodicDirs = {1},
   sourceNodesShared = false, -- default true
   solutionNodesShared = false, -- default true
   writeStiffnessMatrix = false,
   modifierConstant = 0.0,
   isGyroKineticPoisson = true,
   -- boundary conditions to apply
   bcLeft = { T ="N", V = 0.0},
   bcRight = { T ="D_VAR", V = 0.0},
   bcBottom = { T ="D_VAR", V = 0.0},
   bcTop = { T ="D_VAR", V = 0.0},
   bcBack = { T ="D_VAR", V = 0.0},
   bcFront = { T ="D_VAR", V = 0.0},
}
--runUpdater(fieldSolver, 0.0, 0.0, {phiSLower, dirichletBoundaryField}, {phiSLowerSmoothed})

-- df/fx = 0 boundary condition
distfZeroNormalBc = BoundaryCondition.SOLZeroNormal5D {
  basis = basis_ion_5d,
}
distfZeroNormalBcUpdater = Updater.Bc5D {
  onGrid = grid_ion_5d,
  basis = basis_ion_5d,
  boundaryConditions = {distfZeroNormalBc},
  dir = 0, -- direction to apply
  edge = "lower"
}
runUpdater(distfZeroNormalBcUpdater, 0.0, 0.0, {}, {fIon})

-- test field average calculation
fieldInt3d = Updater.IntegrateNodalField3D {
  onGrid = grid_3d,
  basis = basis_3d,
}
-- Stores value of phi_s(x = x_uppper)
avgDensity = DataStruct.DynVector { numComponents = 1, }
bField3d:clear(1.0)
runUpdater(fieldInt3d, 0.0, 0.0, {bField3d}, {avgDensity})
print(avgDensity:lastInsertedData())
-- Compute diagnostics for t = 0
--calcDiagnostics(tCurr, 0.0)
--writeFields(startFrame-1,tCurr)

--for frame = startFrame, nFrames do
--  Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
--  dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
--  tCurr = tCurr+tFrame
--  writeFields(frame, tCurr)
--  Lucee.logInfo ("")
--end
