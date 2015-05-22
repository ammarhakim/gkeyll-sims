-- Input file for ETG test problem
-- 4-4-2015: Split original full-F simulation into a delta F simulation
-- 4-10-2015: Added free energy calculation. Haven't debugged adjoint looping yet.
-- 4-22-2015: More testing
-- 4-23-2015: parallel version
-- 5-1-2015: debugged known issues so far
-- 5-11-2015: same as adjoint3.lua, but with different decomposition
-- 5-11-2015: cos perturbation with 2x k
-- 5-14-2015: lua code to determine iteration count based on accuracy
-- 5-22-2015: optimization times 0.2,0.4,0.6,0.8,1.0

-- phase-space decomposition
phaseDecomp = DecompRegionCalc4D.CartProd { cuts = {4, 4, 8, 1} }
-- configuration space decomposition
confDecomp = DecompRegionCalc2D.SubCartProd4D {
   decomposition = phaseDecomp,
   collectDirections = {0, 1},
}
-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.05
-- parameters to control time-stepping
tStart = 0.0
tEnd = 1e-6
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
iterTol = 0.0001 -- desired tolerance for convergence

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
ky_min    = 2*2*math.pi/deltaR
-- grid parameters: number of cells
N_X = 8
N_Y = 8
N_VPARA = 16
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = R
X_UPPER = R + deltaR
Y_LOWER = -deltaR/2
Y_UPPER = deltaR/2
VPARA_UPPER = math.min(4, 2.5*math.sqrt(N_VPARA/4))*vtKinetic
VPARA_LOWER = -VPARA_UPPER
MU_LOWER = 0
MU_UPPER = math.min(8, 4*math.sqrt(N_MU/2))*kineticMass*vtKinetic*vtKinetic/B0

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
   decomposition = phaseDecomp,
}
grid_back_4d = Grid.RectCart4D {
   lower = {X_LOWER, Y_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_VPARA, N_MU},
   decomposition = phaseDecomp,
}
-- 2d spatial grid
grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
   periodicDirs = {0, 1},
   decomposition = confDecomp,
}
grid_back_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
   decomposition = confDecomp,
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

-- perturbed distribution function for electrons
f = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- for RK time-stepping
f1 = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
f1Intermediate = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- updated solution
fNew = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- for use in time-stepping
fNewDup = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- to store background distfElc
fBackground = DataStruct.Field4D {
   onGrid = grid_back_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
fInitialPerturb = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
fFull = DataStruct.Field4D {
   onGrid = grid_back_4d,
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
  --return 1e-3*(vtKinetic/omega_s)/L_T*math.cos(ky_min*y)
  local x0 = (X_LOWER+X_UPPER)/2
  local y0 = (Y_LOWER+Y_UPPER)/2
  local sigma = deltaR/4
  return 1e-3*(vtKinetic/omega_s)/L_T*math.exp(-(y-y0)^2/(2*sigma^2))*math.exp(-(x-x0)^2/(2*sigma^2))
end

function fProfile(x,y,v,mu)
  return (2*math.pi*kineticTempProfile(x)*eV/kineticMass)^(-3/2)*
    math.exp(-kineticMass*v^2/(2*kineticTempProfile(x)*eV))*
    math.exp(-math.abs(mu)*bFieldProfile(x)/(kineticTempProfile(x)*eV))
end

-- initialize electron distribution function
initKineticF = Updater.EvalOnNodes4D {
   onGrid = grid_back_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,v,mu,t)
		 return n0*fProfile(x,y,v,mu)
	 end
}
runUpdater(initKineticF, 0.0, 0.0, {}, {f})

-- initialize perturbation to electron distribution function
initKineticFPerturb = Updater.EvalOnNodes4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,v,mu,t)
		 return perturbDensityProfile(x,y,v,mu)
	 end
}
runUpdater(initKineticFPerturb, 0.0, 0.0, {}, {fInitialPerturb})
fInitialPerturb:sync()

-- Magnetic Field (2D)
bField2d = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- Magnetic Field (4D)
bField4d = DataStruct.Field4D {
   onGrid = grid_back_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bfield2d
initb = Updater.EvalOnNodes2D {
  onGrid = grid_back_2d,
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
bField4d:sync()

-- Jacobian Factor (4D)
jacobianField = DataStruct.Field4D {
   onGrid = grid_back_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
initJacobian = Updater.EvalOnNodes4D {
  onGrid = grid_back_4d,
  basis = basis_4d,
  shareCommonNodes = false,
  evaluate = function (x,y,t)
    return kineticMass*kineticMass*bFieldProfile(x)
  end
}
-- Fill out jacobian
runUpdater(initJacobian, 0.0, 0.0, {}, {jacobianField})
jacobianField:sync()

-- B_y^* field (4D)
bStarYField = DataStruct.Field4D {
   onGrid = grid_back_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarY
initbStarY = Updater.EvalOnNodes4D {
  onGrid = grid_back_4d,
  basis = basis_4d,
  shareCommonNodes = false,
  evaluate = function (x,y,vPara,mu,t)
    return -kineticMass*vPara/(kineticCharge*x)
  end
}
runUpdater(initbStarY, 0.0, 0.0, {}, {bStarYField})
bStarYField:sync()

-- define equation to solve
gyroEqn = PoissonBracketEquation.GyroEquation4D {
  speciesMass = kineticMass,
  speciesCharge = kineticCharge,
  bStarY = bStarYField,
}
-- This one takes dF and H0
pbSlvr = Updater.PoissonBracketOpt4D {
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

-- Perturbed Hamiltonian
hamilPerturbed = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Duplicate of hamil for reverting
hamilDup = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Static hamiltonian (KE only)
hamilBackground = DataStruct.Field4D {
   onGrid = grid_back_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize KE part of hamil
initHamilBackground = Updater.EvalOnNodes4D {
   onGrid = grid_back_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   evaluate = function (x,y,vPara,mu,t)
      return 0.5*kineticMass*vPara*vPara + math.abs(mu)*bFieldProfile(x)
   end
}
runUpdater(initHamilBackground, 0.0, 0.0, {}, {hamilBackground})
hamilBackground:sync()

-- to store number density
numDensityKinetic = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
numDensityKineticBackground = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
numDensityKineticPerturbed = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
numDensityAdiabatic = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
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
vPara = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}

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
muMoment = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
muMomentTimesB = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}

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
vParaSq = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}

kineticTempField = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
backgroundKineticTemp = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
backgroundKineticTemp4d = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
adjointPotential = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
adjointPotentialSmoothed = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
adjointPotential4d = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Duplicate of adjoint potential for reverting
adjointPotential4dDup = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}

-- to store the electrostatic potential on spatial grid
phi2d = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
phi2dSmoothed = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- to store electrostatic potential for addition to hamiltonian
phi4dSmoothed = DataStruct.Field4D {
  onGrid = grid_4d,
  numComponents = basis_4d:numNodes(),
  ghost = {1, 1},
}
phi4dSmoothedDup = DataStruct.Field4D {
  onGrid = grid_4d,
  numComponents = basis_4d:numNodes(),
  ghost = {1, 1},
}

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

function calcPerturbedHamiltonian(phi2dIn, hamilOut)
  -- copy 2d potential to 4d output field
  runUpdater(copy2dTo4d, 0.0, 0.0, {phi2dIn}, {hamilOut})
  -- scale by charge
  hamilOut:scale(kineticCharge)
  hamilOut:sync()
end

-- dynvector for field energy
fieldEnergy = DataStruct.DynVector { numComponents = 1, }
-- to compute field energy
fieldEnergyCalc = Updater.IntegrateGeneralField2D {
   onGrid = grid_2d,
   basis = basis_2d,
   moment = 2,
}

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

-- adjoint equation source
adjointSource = Updater.ETGAdjointSource4D {
   onGrid = grid_4d,
   -- basis functions to use
   basis = basis_4d,
   -- CFL number
   cfl = cfl,
   -- let solver know about additional jacobian factor
   jacobianField = bField4d,
   --updateDirections = {0,1,2},
   --zeroFluxDirections = {2,3},
   onlyIncrement = true,
   kineticMass = kineticMass,
   eV = eV,
   tauOverZi = adiabaticTemp/kineticTemp,
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

-- Return temperature (in eV?)
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

function calcPotential(phiOut, fIn)
  calcNumDensity(fIn, numDensityKineticPerturbed)
  phiOut:copy(numDensityKineticPerturbed)
  phiOut:scale(-adiabaticTemp/n0)
  phiOut:sync()
end

-- calculate potential for adjoint update
function calcAdjointPotential(curr, dt, fIn, adjointPotentialOut)
  -- basically, calculate the v^2 moment of the perturbed electron distribution
  -- function
  
  -- <2B*mu/m_e>
  runUpdater(muCalc, curr, dt, {fIn, bField2d}, {muMoment})
  muMoment:scale(2*math.pi/kineticMass*2/kineticMass)
  runUpdater(multiply2dCalc, curr, dt, {bField2d, muMoment}, {muMomentTimesB})

  -- <v_para^2>
  runUpdater(vParaSqCalc, curr, dt, {fIn, bField2d}, {vParaSq})
  vParaSq:scale(2*math.pi/kineticMass)

  adjointPotentialOut:clear(0.0)
  adjointPotentialOut:accumulate(1.0, muMomentTimesB)
  adjointPotentialOut:accumulate(1.0, vParaSq)
  adjointPotentialOut:scale(-adiabaticTemp/n0)
  adjointPotentialOut:sync()

  runUpdater(smoothCalc, curr, dt, {adjointPotentialOut}, {adjointPotentialSmoothed})
  adjointPotentialOut:copy(adjointPotentialSmoothed)
  adjointPotentialOut:sync()
end

-- (in eV?)
function calcFreeEnergy(curr, dt, fIn, freeEnergyOut)
  calcNumDensity(fIn, numDensityKineticPerturbed)
  --calcTemperature(fIn, curr, dt, kineticTempField)
  runUpdater(freeEnergyCalc, curr, dt, {bField2d, backgroundKineticTemp,
    numDensityKineticBackground, numDensityKineticPerturbed, fBackground, fIn}, {freeEnergyOut})
end

-- function to take a time-step using SSP-RK3 time-stepping scheme (ADJOINT STEP)
function rk3adjoint(tCurr, myDt)
  -- RK stage 1
  f1:copy(f)
  local myStatusOne, myDtSuggestedOne = runUpdater(pbSlvr, tCurr, myDt, {f, hamilBackground}, {f1Intermediate})
  f1:accumulate(-myDt, f1Intermediate)
  local myStatusTwo, myDtSuggestedTwo = runUpdater(pbSlvr, tCurr, myDt, {fBackground, hamilPerturbed}, {f1Intermediate})
  f1:accumulate(-myDt, f1Intermediate)
  local myStatusThree, myDtSuggestedThree = runUpdater(adjointSource, tCurr, myDt,
    {f1, phi4dSmoothed, adjointPotential4d, backgroundKineticTemp4d, fBackground, bField4d}, {f1Intermediate})
  f1:accumulate(-myDt, f1Intermediate)

  if (myStatusOne == false or myStatusTwo == false) then
    return false, math.min(myDtSuggestedOne, myDtSuggestedTwo)
  end
  f1:sync()
  -- compute potential
  calcPotential(phi2d, f1)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi2d}, {phi2dSmoothed})
  phi2dSmoothed:sync()
  -- calculate adjoint potential
  calcAdjointPotential(tCurr, myDt, f1, adjointPotential)
  runUpdater(copy2dTo4d, tCurr, myDt, {adjointPotential}, {adjointPotential4d})
  runUpdater(copy2dTo4d, tCurr, myDt, {phi2dSmoothed}, {phi4dSmoothed})
  adjointPotential4d:sync()
  phi4dSmoothed:sync()
  -- Compute hamiltonian
  calcPerturbedHamiltonian(phi2dSmoothed, hamilPerturbed)

  -- RK stage 2
  fNew:copy(f1)
  local myStatusOne, myDtSuggestedOne = runUpdater(pbSlvr, tCurr, myDt, {f1, hamilBackground}, {f1Intermediate})
  fNew:accumulate(-myDt, f1Intermediate)
  local myStatusTwo, myDtSuggestedTwo = runUpdater(pbSlvr, tCurr, myDt, {fBackground, hamilPerturbed}, {f1Intermediate})
  fNew:accumulate(-myDt, f1Intermediate)
  local myStatusThree, myDtSuggestedThree = runUpdater(adjointSource, tCurr, myDt,
    {fNew, phi4dSmoothed, adjointPotential4d, backgroundKineticTemp4d, fBackground, bField4d}, {f1Intermediate})
  fNew:accumulate(-myDt, f1Intermediate)

  if (myStatusOne == false or myStatusTwo == false) then
    return false, math.min(myDtSuggestedOne, myDtSuggestedTwo)
  end

  f1:combine(3.0/4.0, f, 1.0/4.0, fNew)
  f1:sync()
  -- compute potential
  calcPotential(phi2d, f1)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi2d}, {phi2dSmoothed})
  phi2dSmoothed:sync()
  -- calculate adjoint potential
  calcAdjointPotential(tCurr, myDt, f1, adjointPotential)
  runUpdater(copy2dTo4d, myDt, myDt, {adjointPotential}, {adjointPotential4d})
  runUpdater(copy2dTo4d, tCurr, myDt, {phi2dSmoothed}, {phi4dSmoothed})
  adjointPotential4d:sync()
  phi4dSmoothed:sync()
  -- Compute hamiltonian
  calcPerturbedHamiltonian(phi2dSmoothed, hamilPerturbed)

  -- RK stage 3
  fNew:copy(f1)
  local myStatusOne, myDtSuggestedOne = runUpdater(pbSlvr, tCurr, myDt, {f1, hamilBackground}, {f1Intermediate})
  fNew:accumulate(-myDt, f1Intermediate)
  local myStatusTwo, myDtSuggestedTwo = runUpdater(pbSlvr, tCurr, myDt, {fBackground, hamilPerturbed}, {f1Intermediate})
  fNew:accumulate(-myDt, f1Intermediate)
  local myStatusThree, myDtSuggestedThree = runUpdater(adjointSource, tCurr, myDt,
    {fNew, phi4dSmoothed, adjointPotential4d, backgroundKineticTemp4d, fBackground, bField4d}, {f1Intermediate})
  fNew:accumulate(-myDt, f1Intermediate)

  if (myStatusOne == false or myStatusTwo == false) then
    return false, math.min(myDtSuggestedOne, myDtSuggestedTwo)
  end

  f1:combine(1.0/3.0, f, 2.0/3.0, fNew)
  f1:sync()
  -- compute potential
  calcPotential(phi2d, f1)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi2d}, {phi2dSmoothed})
  phi2dSmoothed:sync()
  -- calculate adjoint potential
  calcAdjointPotential(tCurr, myDt, f1, adjointPotential)
  runUpdater(copy2dTo4d, myDt, myDt, {adjointPotential}, {adjointPotential4d})
  runUpdater(copy2dTo4d, tCurr, myDt, {phi2dSmoothed}, {phi4dSmoothed})
  adjointPotential4d:sync()
  phi4dSmoothed:sync()
  -- Compute hamiltonian
  calcPerturbedHamiltonian(phi2dSmoothed, hamilPerturbed)

  f:copy(f1)

  return true, math.min(myDtSuggestedOne, myDtSuggestedTwo)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
  -- RK stage 1
  f1:copy(f)
  local myStatusOne, myDtSuggestedOne = runUpdater(pbSlvr, tCurr, myDt, {f, hamilBackground}, {f1Intermediate})
  f1:accumulate(myDt, f1Intermediate)
  local myStatusTwo, myDtSuggestedTwo = runUpdater(pbSlvr, tCurr, myDt, {fBackground, hamilPerturbed}, {f1Intermediate})
  f1:accumulate(myDt, f1Intermediate)

  if (myStatusOne == false or myStatusTwo == false) then
    return false, math.min(myDtSuggestedOne, myDtSuggestedTwo)
  end
  f1:sync()
  -- compute potential
  calcPotential(phi2d, f1)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi2d}, {phi2dSmoothed})
  phi2dSmoothed:sync()
  -- Compute hamiltonian
  calcPerturbedHamiltonian(phi2dSmoothed, hamilPerturbed)

  -- RK stage 2
  fNew:copy(f1)
  local myStatusOne, myDtSuggestedOne = runUpdater(pbSlvr, tCurr, myDt, {f1, hamilBackground}, {f1Intermediate})
  fNew:accumulate(myDt, f1Intermediate)
  local myStatusTwo, myDtSuggestedTwo = runUpdater(pbSlvr, tCurr, myDt, {fBackground, hamilPerturbed}, {f1Intermediate})
  fNew:accumulate(myDt, f1Intermediate)

  if (myStatusOne == false or myStatusTwo == false) then
    return false, math.min(myDtSuggestedOne, myDtSuggestedTwo)
  end

  f1:combine(3.0/4.0, f, 1.0/4.0, fNew)
  f1:sync()
  -- compute potential
  calcPotential(phi2d, f1)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi2d}, {phi2dSmoothed})
  phi2dSmoothed:sync()
  -- Compute hamiltonian
  calcPerturbedHamiltonian(phi2dSmoothed, hamilPerturbed)

  -- RK stage 3
  fNew:copy(f1)
  local myStatusOne, myDtSuggestedOne = runUpdater(pbSlvr, tCurr, myDt, {f1, hamilBackground}, {f1Intermediate})
  fNew:accumulate(myDt, f1Intermediate)
  local myStatusTwo, myDtSuggestedTwo = runUpdater(pbSlvr, tCurr, myDt, {fBackground, hamilPerturbed}, {f1Intermediate})
  fNew:accumulate(myDt, f1Intermediate)

  if (myStatusOne == false or myStatusTwo == false) then
    return false, math.min(myDtSuggestedOne, myDtSuggestedTwo)
  end

  f1:combine(1.0/3.0, f, 2.0/3.0, fNew)
  f1:sync()
  -- compute potential
  calcPotential(phi2d, f1)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi2d}, {phi2dSmoothed})
  phi2dSmoothed:sync()
  -- Compute hamiltonian
  calcPerturbedHamiltonian(phi2dSmoothed, hamilPerturbed)

  f:copy(f1)

  return true, math.min(myDtSuggestedOne, myDtSuggestedTwo)
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
  local step = 1
  local tCurr = tStart
  local myDt = initDt
  local status, dtSuggested

  while tCurr<=tEnd do
    -- Store fields that might need to be reverted
    fNewDup:copy(f)
    hamilDup:copy(hamilPerturbed)

    -- if needed adjust dt to hit tEnd exactly
    if (tCurr+myDt > tEnd) then
      myDt = tEnd-tCurr
    end

    Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
    status, dtSuggested = rk3(tCurr, myDt)

    if (status == false) then
      -- time-step too large
      Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
      -- Revert fields to previous value
      f:copy(fNewDup)
      hamilPerturbed:copy(hamilDup)

      myDt = dtSuggested
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

  return dtSuggested
end

-- function to advance solution from tStart to tEnd
function advanceFrameAdjoint(tStart, tEnd, initDt)
  local step = 1
  local tCurr = tStart
  local myDt = initDt
  local status, dtSuggested

  while tCurr<=tEnd do
    -- Store fields that might need to be reverted
    fNewDup:copy(f)
    hamilDup:copy(hamilPerturbed)
    adjointPotential4dDup:copy(adjointPotential4d)
    phi4dSmoothedDup:copy(phi4dSmoothed)

    -- if needed adjust dt to hit tEnd exactly
    if (tCurr+myDt > tEnd) then
      myDt = tEnd-tCurr
    end

    Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
    status, dtSuggested = rk3adjoint(tCurr, myDt)

    if (status == false) then
      -- time-step too large
      Lucee.logInfo (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
      -- Revert fields to previous value
      f:copy(fNewDup)
      hamilPerturbed:copy(hamilDup)
      adjointPotential4d:copy(adjointPotential4dDup)
      phi4dSmoothed:copy(phi4dSmoothedDup)

      myDt = dtSuggested
    else
      tCurr = tCurr + myDt
      myDt = dtSuggested
      step = step + 1

      if (tCurr >= tEnd) then
        break
      end
    end
  end

  return dtSuggested
end

-- Compute initial kinetic density
calcNumDensity(f, numDensityKinetic)
-- Scale distribution function and apply bcs
runUpdater(scaleInitDistF, 0.0, 0.0, {numDensityKinetic}, {f})
--f:sync()
-- Recalculate number density
calcNumDensity(f, numDensityKinetic)
-- Store static numDensityAdiabatic field
numDensityAdiabatic:copy(numDensityKinetic)
numDensityKineticBackground:copy(numDensityKinetic)
-- Store background f
fBackground:copy(f)
-- Compute background f's temperature
calcBackgroundTemperature(fBackground, 0.0, 0.0, backgroundKineticTemp)
runUpdater(copy2dTo4d, 0.0, 0.0, {backgroundKineticTemp}, {backgroundKineticTemp4d})

-- Compute total perturbed distribution function using scaled fBackground
runUpdater(multiply4dCalc, 0.0, 0.0, {f, fInitialPerturb}, {fNew})
-- Set f to perturbed f
f:copy(fNew)
-- Apply boundary conditions
f:sync()
-- Copy to fInitialPetrub again
fInitialPerturb:copy(f)

-- function to take one complete step of adjoint algorithm
function adjointSingleIteration(currIter, tEndIter)
  -- Compute initial potential with perturbation added
  calcPotential(phi2d, f)
  runUpdater(smoothCalc, 0.0, 0.0, {phi2d}, {phi2dSmoothed})
  phi2dSmoothed:sync()
  calcPerturbedHamiltonian(phi2dSmoothed, hamilPerturbed)
  -- calculate free energy at time 0
  if (currIter == 0) then
    calcFreeEnergy(2*currIter, 0.0, f, tempFreeEnergy)
    W_0 = tempFreeEnergy:lastInsertedData()
  end
  calcDiagnostics(0.0, 0.0)
  -- perform standard iteration to time tEndIter
  Lucee.logInfo (string.format("--Iteration %g: Advancing solution from %g to %g", currIter, tStart, tEndIter))
  dtSuggested = advanceFrame(tStart, tEndIter, dtSuggested)
  Lucee.logInfo ("")
  W_T_PREV = W_T
  W_T = freeEnergy:lastInsertedData()
  
  -- recompute initial potential and hamiltonian with scaled f
  calcPotential(phi2d, f)
  runUpdater(smoothCalc, 0.0, 0.0, {phi2d}, {phi2dSmoothed})
  phi2dSmoothed:sync()
  calcPerturbedHamiltonian(phi2dSmoothed, hamilPerturbed)
  -- calculate initial adjoint potential
  calcAdjointPotential(0.0, 0.0, f, adjointPotential)
  runUpdater(copy2dTo4d, 0.0, 0.0, {adjointPotential}, {adjointPotential4d})
  runUpdater(copy2dTo4d, 0.0, 0.0, {phi2dSmoothed}, {phi4dSmoothed})
  adjointPotential4d:sync()
  phi4dSmoothed:sync()

  Lucee.logInfo (string.format("--Iteration %g: (Adjoint) Advancing solution from %g to %g", currIter, tEndIter, tStart))
  dtSuggested = advanceFrameAdjoint(tStart, tEndIter, dtSuggested)
  Lucee.logInfo ("")
  
  calcFreeEnergy(2*currIter+1, 0.0, f, tempFreeEnergy)
  W_0_NEW = tempFreeEnergy:lastInsertedData()
  -- calculate f at time 0
  f:scale(math.sqrt(W_0/W_0_NEW))
  
  -- check for convergence
  if currIter == 0 or math.abs((W_T-W_T_PREV)/W_T_PREV) > iterTol then
    if currIter ~= 0 then
      Lucee.logInfo (string.format("-- Convergence factor = %g", math.abs((W_T-W_T_PREV)/W_T_PREV)))
    end
    return false
  else
    return true
  end
end

timePoints = 5 -- number of points to optimize free energy growth
-- build array of times to optimize free energy growth
tEndArray = {}
for i = 1,timePoints do
  tEndArray[i] = i*tEnd/timePoints
end

for i,tEndVal in ipairs(tEndArray) do

  iterCounter = 0
  repeat
    convergenceStatus = adjointSingleIteration(iterCounter, tEndVal)
    -- need to write to temporary file to clear freeEnergy dynvector
    freeEnergy:write( string.format("freeEnergy.h5"), tEndVal)
    --freeEnergy:write( string.format("freeEnergy_%d.h5", iterCounter), tEndVal)
    iterCounter = iterCounter + 1
  until convergenceStatus == true
  -- final iteration to write out a free energy vs time curve for optimized curve
  calcDiagnostics(0.0, 0.0)
  calcPotential(phi2d, f)
  runUpdater(smoothCalc, 0.0, 0.0, {phi2d}, {phi2dSmoothed})
  phi2dSmoothed:sync()
  calcPerturbedHamiltonian(phi2dSmoothed, hamilPerturbed)
  -- perform standard iteration to time tEnd
  Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tStart, tEnd))
  dtSuggested = advanceFrame(tStart, tEnd, dtSuggested)
  Lucee.logInfo ("")
  freeEnergy:write( string.format("freeEnergy_%d.h5", i-1), tEnd)

  -- Revert f to initial condition before any iterations
  f:copy(fInitialPerturb)
end
