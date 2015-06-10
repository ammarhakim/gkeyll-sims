-- Input file for ETG test problem
-- Species are referred to as the 'kinetic' or 'adiabatic' specie
-- 4-4-2015: Split original full-F simulation into a delta F simulation
-- 4-10-2015: Added free energy calculation. Haven't debugged adjoint looping yet.
-- 4-22-2015: More testing
-- 4-25-2015: Starting with eigenmode
-- 5-1-2015: possibly another guess at eigenmode
-- 6-3-2015: new simulation.
-- 6-5-2015: constructs inner product matrix

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.05
-- parameters to control time-stepping
tStart = 0.0
tEnd = 1e-7
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
iterTotal = 40
nFrames = 25
tFrame = (tEnd-tStart)/nFrames -- time between frames

-- physical parameters
eV            = Lucee.ElementaryCharge
kineticCharge = -Lucee.ElementaryCharge
adiabaticCharge = Lucee.ElementaryCharge
kineticMass   = Lucee.ElectronMass
adiabaticMass = 2.014*Lucee.ProtonMass -- (deuterium ions)
kineticTemp   = 2072 -- [eV]
adiabaticTemp = 2072 -- [ev]
B0 = 1.91   -- [T]
R0 = 1.313  -- [m]
a  = 0.4701 -- [m]
n0 = 4.992*10^(19) -- [1/m^3]
-- derived parameters
R         = R0 + 0.5*a
vtKinetic = math.sqrt(kineticTemp*eV/kineticMass)
c_s       = math.sqrt(kineticTemp*eV/kineticMass)
omega_s   = math.abs(kineticCharge*B0/kineticMass)
rho_s     = c_s/omega_s
deltaR    = 32*rho_s
L_T       = R/2
ky_min    = 2*math.pi/deltaR
-- grid parameters: number of cells
N_X = 1
N_Y = 8
N_VPARA = 4
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = R
X_UPPER = R + deltaR
Y_LOWER = -deltaR/2
Y_UPPER = deltaR/2
VPARA_UPPER = math.min(4, 2.5*math.sqrt(N_VPARA/4))*vtKinetic
VPARA_LOWER = -VPARA_UPPER
MU_LOWER = 0
MU_UPPER = math.min(8, 4*math.sqrt(N_MU/2))*kineticMass*vtKinetic*vtKinetic/(2*B0)

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

-- GRID NOTES
-- Even though non-fluctuating components are stored on the fluctuating
-- grids (which have periodicity in certain directions), they are never
-- called with sync() which makes this okay.
-- Additionally, grid_back_2d can probably be eliminated since it is 0
-- but I am keeping this in case it is not 0 somehow.

-- full 4d phase space grid
grid_4d = Grid.RectCart4D {
   lower = {X_LOWER, Y_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_VPARA, N_MU},
   periodicDirs = {0, 1},
}
-- 2d spatial grid
grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
   periodicDirs = {0, 1},
}
grid_back_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
   periodicDirs = {1},
}
-- create 4d basis functions
basis_4d = NodalFiniteElement4D.SerendipityElement {
   onGrid = grid_4d,
   polyOrder = polyOrder,
}
-- create 2d basis functions
basis_2d = NodalFiniteElement2D.SerendipityElement {
   onGrid = grid_2d,
   polyOrder = polyOrder,
}

-- distribution function for electrons
f = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
g = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
h = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- to store background distfElc
fBackground = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}

-- to copy 2d field to a 4d field
copy2dTo4d = Updater.NodalCopy2DTo4DFieldUpdater {
  -- 4D phase-space grid 
  onGrid = grid_4d,
  -- 4D phase-space basis functions
  basis4d = basis_4d,
  -- 2D spatial basis functions
  basis2d = basis_2d,
  -- Basis function order
  polyOrder = polyOrder,
}

function bFieldProfile(x)
  return B0*R/x
end

function kineticTempProfile(x)
  return kineticTemp*(1 - (x-R)/L_T)
end

function perturbDensityProfile(x,y,v,mu)
  return 1e-3*(vtKinetic/omega_s)/L_T*( math.cos(ky_min*y)
  - math.sqrt(2)*((kineticMass*v^2 + 2*mu*bFieldProfile(x))/(2*kineticTempProfile(x)*eV) - 3/2)*math.sin(ky_min*y) ) 
end

function fProfile(x,y,v,mu)
  return (2*math.pi*kineticTempProfile(x)*eV/kineticMass)^(-3/2)*
    math.exp(-kineticMass*v^2/(2*kineticTempProfile(x)*eV))*
    math.exp(-math.abs(mu)*bFieldProfile(x)/(kineticTempProfile(x)*eV))
end

-- initialize electron distribution function
initKineticF = Updater.EvalOnNodes4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,v,mu,t)
		 return n0*fProfile(x,y,v,mu)
	 end
}
runUpdater(initKineticF, 0.0, 0.0, {}, {f})

-- Magnetic Field (2D)
bField2d = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- Magnetic Field (4D)
bField4d = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bfield2d
initb = Updater.EvalOnNodes2D {
  onGrid = grid_2d,
  basis = basis_2d,
  shareCommonNodes = false,
  evaluate = function (x,y,t)
    return bFieldProfile(x)
  end
}
runUpdater(initb, 0.0, 0.0, {}, {bField2d})
bField2d:sync()
-- Copy magnetic field to 4d
runUpdater(copy2dTo4d, 0.0, 0.0, {bField2d}, {bField4d})

-- Jacobian Factor (4D)
jacobianField = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
initJacobian = Updater.EvalOnNodes4D {
  onGrid = grid_4d,
  basis = basis_4d,
  shareCommonNodes = false,
  evaluate = function (x,y,t)
    return kineticMass*kineticMass*bFieldProfile(x)
  end
}
-- Fill out jacobian
runUpdater(initJacobian, 0.0, 0.0, {}, {jacobianField})

-- B_y^* field (4D)
bStarYField = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarY
initbStarY = Updater.EvalOnNodes4D {
  onGrid = grid_4d,
  basis = basis_4d,
  shareCommonNodes = false,
  evaluate = function (x,y,vPara,mu,t)
    return -kineticMass*vPara/(kineticCharge*x)
  end
}
runUpdater(initbStarY, 0.0, 0.0, {}, {bStarYField})

-- define equation to solve
gyroEqn = PoissonBracketEquation.GyroEquation4D {
  speciesMass = kineticMass,
  speciesCharge = kineticCharge,
  bStarY = bStarYField,
}
-- This one takes dF and H0
pbSlvrOne = Updater.PoissonBracketOpt4D {
   onGrid = grid_4d,
   -- basis functions to use
   basis = basis_4d,
   -- CFL number
   cfl = cfl,
   -- equation to solve
   equation = gyroEqn,
   -- let solver know about additional jacobian factor
   jacobianField = jacobianField,
   updateDirections = {0,1,2},
   zeroFluxDirections = {2,3},
   onlyIncrement = true,
   fluxType = "central",
}

-- This one takes F0 and dH
pbSlvrTwo = Updater.PoissonBracketOpt4D {
   onGrid = grid_4d,
   -- basis functions to use
   basis = basis_4d,
   -- CFL number
   cfl = cfl,
   -- equation to solve
   equation = gyroEqn,
   -- let solver know about additional jacobian factor
   jacobianField = jacobianField,
   updateDirections = {0,1,2},
   zeroFluxDirections = {2,3},
   onlyIncrement = true,
   fluxType = "central"
}

-- Perturbed Hamiltonian
hamilPerturbed = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Duplicate of hamil for reverting
hamilDup = hamilPerturbed:duplicate()
-- Static hamiltonian (KE only)
hamilBackground = hamilPerturbed:duplicate()
-- Updater to initialize KE part of hamil
initHamilBackground = Updater.EvalOnNodes4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   evaluate = function (x,y,vPara,mu,t)
      return 0.5*kineticMass*vPara*vPara + math.abs(mu)*bFieldProfile(x)
   end
}
runUpdater(initHamilBackground, 0.0, 0.0, {}, {hamilBackground})

-- to store number density
numDensityKinetic = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
gNumDensity = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
hNumDensity = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
numDensityKineticBackground = numDensityKinetic:duplicate()
numDensityKineticPerturbed = numDensityKinetic:duplicate()
numDensityAdiabatic = numDensityKinetic:duplicate()
-- to compute number density
numDensityCalc = Updater.DistFuncMomentCalcWeighted2D {
   -- 4D phase-space grid 
   onGrid = grid_4d,
   -- 4D phase-space basis functions
   basis4d = basis_4d,
   -- 2D spatial basis functions
   basis2d = basis_2d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}

-- to compute vPara moment
vParaCalc = Updater.DistFuncMomentCalcWeighted2D {
   -- 4D phase-space grid 
   onGrid = grid_4d,
   -- 4D phase-space basis functions
   basis4d = basis_4d,
   -- 2D spatial basis functions
   basis2d = basis_2d,
   -- desired moment (0, 1 or 2)
   moment = 1,
   -- direction to calculate moment
   momentDirection = 2,
}
vPara = numDensityKinetic:duplicate()

-- to compute mu moment
muCalc = Updater.DistFuncMomentCalcWeighted2D {
   -- 4D phase-space grid 
   onGrid = grid_4d,
   -- 4D phase-space basis functions
   basis4d = basis_4d,
   -- 2D spatial basis functions
   basis2d = basis_2d,
   -- desired moment (0, 1 or 2)
   moment = 1,
   -- direction to calculate moment
   momentDirection = 3,
}
muMoment = numDensityKinetic:duplicate()
muMomentTimesB = numDensityKinetic:duplicate()

-- to compute vParaSq moment
vParaSqCalc = Updater.DistFuncMomentCalcWeighted2D {
   -- 4D phase-space grid 
   onGrid = grid_4d,
   -- 4D phase-space basis functions
   basis4d = basis_4d,
   -- 2D spatial basis functions
   basis2d = basis_2d,
   -- desired moment (0, 1 or 2)
   moment = 2,
   -- direction to calculate moment
   momentDirection = 2,
}
vParaSq = numDensityKinetic:duplicate()

kineticTempField = numDensityKinetic:duplicate()
backgroundKineticTemp = numDensityKinetic:duplicate()

adjointPotential = numDensityKinetic:duplicate()
adjointPotentialSmoothed = numDensityKinetic:duplicate()
adjointPotential4d = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Duplicate of adjoint potential for reverting
adjointPotential4dDup = adjointPotential4d:duplicate()

-- to store the electrostatic potential on spatial grid
phi2d = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
phi2dSmoothedForG = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
phi2dSmoothedForH = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}

phi2dBackground = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
phi2dFluctuating = phi2d:duplicate()
-- to store electrostatic potential for addition to hamiltonian
phi4dSmoothed = DataStruct.Field4D {
  onGrid = grid_4d,
  numComponents = basis_4d:numNodes(),
  ghost = {1, 1},
}
phi4dSmoothedDup = phi4dSmoothed:duplicate()

-- Updater to compute potential
phiCalc = Updater.ETGAdiabaticPotentialUpdater {
  onGrid = grid_2d,
  basis = basis_2d,
  n0 = n0,
  adiabaticTemp = adiabaticTemp,
  adiabaticCharge = adiabaticCharge,
}

-- Updater to smooth out 2d field (phi)
smoothCalc = Updater.SimpleSmoothToC02D {
   onGrid = grid_back_2d,
   basis = basis_2d,
   polyOrder = polyOrder,
}

multiply4dCalc = Updater.FieldArithmeticUpdater4D {
  onGrid = grid_4d,
  basis = basis_4d,
  evaluate = function(fld1, fld2)
    return fld1*fld2
  end
}

multiply2dCalc = Updater.FieldArithmeticUpdater2D {
  onGrid = grid_2d,
  basis = basis_2d,
  evaluate = function(fld1, fld2)
    return fld1*fld2
  end
}

function applyBcToBackgroundPotential(fld)
  fld:applyCopyBc(0, "lower")
  fld:applyCopyBc(0, "upper")
  fld:sync()
end

function applyBcToPerturbedDistF(fld)
  -- Apply periodic boundary conditions to fluctuating component
  fld:sync()
end

function calcPerturbedHamiltonian(phi2dIn, hamilOut)
  -- copy 2d potential to 4d output field
  runUpdater(copy2dTo4d, 0.0, 0.0, {phi2dIn}, {hamilOut})
  -- scale by charge
  hamilOut:scale(kineticCharge)
end

-- dynvector for free energy
freeEnergy = DataStruct.DynVector { numComponents = 1, }
-- used with lastInsertedData() for adjoint iteration
tempFreeEnergy = DataStruct.DynVector { numComponents = 1, }
-- to compute field energy
freeEnergyCalc = Updater.ETGFreeEnergy {
   onGrid = grid_4d,
   basis4d = basis_4d,
   basis2d = basis_2d,
   adiabaticTemp = adiabaticTemp,
   kineticMass = kineticMass,
}

-- to scale distribution function
scaleInitDistF = Updater.ETGInitializeDensity {
  -- 4D phase-space grid 
   onGrid = grid_4d,
   -- 4D phase-space basis functions
   basis4d = basis_4d,
   -- 2D spatial basis functions
   basis2d = basis_2d,
   -- Desired constant density
   constantDensity = n0,
   polyOrder = polyOrder,
}

function calcNumDensity(fIn, nOut)
  runUpdater(numDensityCalc, 0.0, 0.0, {fIn, bField2d}, {nOut})
  nOut:scale(2*math.pi/kineticMass)
end

function calcTotalNumDensity(nBackground, nPerturbed, nTotal)
  nTotal:copy(nBackground)
  nTotal:accumulate(1.0, nPerturbed)
end

function calcDiagnostics(curr, dt)
  --runUpdater(vParaCalc, curr, dt, {f, bField2d}, {vPara})
  --vPara:scale(2*math.pi/kineticMass)
  --runUpdater(vParaSqCalc, curr, dt, {f, bField2d}, {vParaSq})
  --vParaSq:scale(2*math.pi/kineticMass)
  --runUpdater(fieldEnergyCalc, curr, dt, {phi2dSmoothed}, {fieldEnergy})
  calcFreeEnergy(curr, dt, f, freeEnergy)
end

-- Return temperautre (in eV?)
function calcTemperature(fIn, curr, dt, outputField)
  fFull:copy(fBackground)
  fFull:accumulate(fIn)

  runUpdater(muCalc, curr, dt, {fFull, bField2d}, {muMoment})
  muMoment:scale(2*math.pi/kineticMass)
  runUpdater(multiply2dCalc, curr, dt, {bField2d, muMoment}, {muMomentTimesB})

  runUpdater(vParaSqCalc, curr, dt, {fFull, bField2d}, {vParaSq})
  vParaSq:scale(2*math.pi/kineticMass)

  outputField:clear(0.0)
  outputField:accumulate(kineticMass/(3*n0), vParaSq)
  outputField:accumulate(2/(3*n0), muMomentTimesB)
  outputField:scale(1/eV)
end

-- Return temperature of a background distribution function (in eV?)
function calcBackgroundTemperature(fIn, curr, dt, outputField)
  runUpdater(muCalc, curr, dt, {fIn, bField2d}, {muMoment})
  muMoment:scale(2*math.pi/kineticMass)
  runUpdater(multiply2dCalc, curr, dt, {bField2d, muMoment}, {muMomentTimesB})

  runUpdater(vParaSqCalc, curr, dt, {fIn, bField2d}, {vParaSq})
  vParaSq:scale(2*math.pi/kineticMass)

  outputField:clear(0.0)
  outputField:accumulate(kineticMass/(3*n0), vParaSq)
  outputField:accumulate(2/(3*n0), muMomentTimesB)
  outputField:scale(1/eV)
end

function calcPotential(phiOut, numDensityPerturbed)
  phiOut:copy(numDensityPerturbed)
  phiOut:scale(-adiabaticTemp/n0)
  phiOut:sync()
end

-- (in eV?)
function calcFreeEnergy(curr, dt, fIn, freeEnergyOut)
  calcNumDensity(fIn, numDensityKineticPerturbed)
  --calcTemperature(f, curr, dt, kineticTempField)
  runUpdater(freeEnergyCalc, curr, dt, {bField2d, backgroundKineticTemp,
    numDensityKineticBackground, numDensityKineticPerturbed, fBackground, fIn}, {freeEnergyOut})
end

-- Compute initial kinetic density
calcNumDensity(f, numDensityKinetic)
-- Scale distribution function and apply bcs
runUpdater(scaleInitDistF, 0.0, 0.0, {numDensityKinetic}, {f})
-- Recalculate number density
calcNumDensity(f, numDensityKinetic)
-- Store static numDensityAdiabatic field
numDensityAdiabatic:copy(numDensityKinetic)
numDensityKineticBackground:copy(numDensityKinetic)
-- Store background f
fBackground:copy(f)
fBackground:sync()
-- Compute background f's temperature
calcBackgroundTemperature(fBackground, 0.0, 0.0, backgroundKineticTemp)
backgroundKineticTemp:sync()

-- keeps track of what node is being passed to initSingleNode
initIndex = 0
-- initialize perturbation to electron distribution function
initSingleNode = Updater.SetSingleNodeToOne4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   -- function to communicate what node to set to one
   evaluate = function (t)
		 return initIndex
	 end
}

-- figure out how many nodes are in the system
totalNodes = (N_X+2)*(N_Y+2)*(N_VPARA+2)*(N_MU+2)*basis_4d:numNodes()
print(string.format("-- Total nodes = %g",totalNodes))
nodesPerPosition = (N_VPARA+2)*(N_MU+2)*basis_4d:numNodes()

gNodeIndex = 0
hNodeIndex = 0
matrixConstructor = Updater.ETGInnerProduct {
  onGrid = grid_4d,
  basis4d = basis_4d,
  basis2d = basis_2d,
  shareCommonNodes = false,
  adiabaticTemp = adiabaticTemp,
  n0 = n0,
  eV = eV,
  kineticMass = kineticMass,
  totalNodes = totalNodes,
  nodesPerPosition = nodesPerPosition,
  -- function to communicate what column to write to
  evaluate = function (t)
   return gNodeIndex, hNodeIndex
  end
}

-- MAIN LOOP
-- Loop over lower triangular portion of matrix since we know it will be symmetric
for gNode = 0, totalNodes-1 do
  Lucee.logInfo (string.format("-- Step for row %g out of %g", gNode, totalNodes-1))
  -- initialize g distribution function
  gNodeIndex = gNode
  initIndex = gNode
  runUpdater(initSingleNode, 0.0, 0.0, {}, {g})
  -- Apply boundary conditions to distribution function
  g:sync()
  -- Compute initial potential with perturbation added
  calcNumDensity(g, gNumDensity)
  
  for hNode = math.max(0,gNode-nodesPerPosition), gNode do
    hNodeIndex = hNode
    -- initialize h distribution function
    initIndex = hNode
    runUpdater(initSingleNode, 0.0, 0.0, {}, {h})
    -- Apply boundary conditions to distribution function
    h:sync()
    -- Compute initial potential with perturbation added
    calcNumDensity(h, hNumDensity)
    -- Compute inner product
    runUpdater(matrixConstructor, 0.0, 0.0, {bField2d, backgroundKineticTemp, fBackground,
      gNumDensity, hNumDensity, g, h}, {})
  end
end
