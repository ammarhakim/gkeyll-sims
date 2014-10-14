-- Input file for ETG test problem
-- Species are referred to as the 'kinetic' or 'adiabatic' species

--decomp4d = DecompRegionCalc4D.CartProd {cuts = {1,1,1,1}, periodicDirs={0,1}}
--decomp2d = DecompRegionCalc2D.CartProd {cuts = {1,1}}

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.1
-- parameters to control time-stepping
tStart = 0.0
tEnd = 2e-9
--tEnd = 6.60129e-09
dtSuggested = 0.1*1e-6 -- initial time-step to use (will be adjusted)
nFrames = 4
tFrame = (tEnd-tStart)/nFrames -- time between frames

-- physical parameters
eV            = Lucee.ElementaryCharge
kineticCharge = -Lucee.ElementaryCharge
ionCharge     = Lucee.ElementaryCharge
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
L_T       = R/10
ky_min    = 2*math.pi/deltaR
-- grid parameters: number of cells
N_X = 4
N_Y = 8
N_VPARA = 4
N_MU = 2
-- grid parameters: domain extent
X_LOWER = R
X_UPPER = R + deltaR
Y_LOWER = -deltaR/2
Y_UPPER = deltaR/2
VPARA_LOWER = -2*vtKinetic
VPARA_UPPER = 2*vtKinetic
MU_LOWER = 0
MU_UPPER = 4*kineticMass*VPARA_UPPER*VPARA_UPPER/(2*B0)

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

-- full 4d phase space grid
grid_4d = Grid.RectCart4D {
   lower = {X_LOWER, Y_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_VPARA, N_MU},
}
grid_fluct_4d = Grid.RectCart4D {
   lower = {X_LOWER, Y_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_VPARA, N_MU},
   periodicDirs = {0, 1},
}
grid_back_4d = Grid.RectCart4D {
   lower = {X_LOWER, Y_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_VPARA, N_MU},
   periodicDirs = {1},
}
-- 2d spatial grid
grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
}
grid_fluct_2d = Grid.RectCart2D {
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
   onGrid = grid_back_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- for RK time-stepping
f1 = f:duplicate()
-- updated solution
fNew = f:duplicate()
-- for use in time-stepping
fNewDup = fNew:duplicate()
-- to store background distfElc
fBackground = f:duplicate()
fFluctuating = DataStruct.Field4D {
   onGrid = grid_fluct_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
fInitialPerturb = f:duplicate()

function bFieldProfile(x)
  return B0*R/x
end

function kineticTempProfile(x)
  return kineticTemp*(1 - (x-R)/L_T)
end

function perturbDensityProfile(x,y)
  return n0*1e-3*(vtKinetic/omega_s)/L_T*math.cos(ky_min*y)
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

-- initialize electron distribution function
initKineticFPerturb = Updater.EvalOnNodes4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,v,mu,t)
		 return perturbDensityProfile(x,y)*fProfile(x,y,v,mu)
	 end
}
runUpdater(initKineticFPerturb, 0.0, 0.0, {}, {fInitialPerturb})

-- Magnetic Field (2D)
bField2d = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
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

-- to copy phi to a 4d field
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
    return kineticMass*vPara/(kineticCharge*x)
  end
}
runUpdater(initbStarY, 0.0, 0.0, {}, {bStarYField})

-- define equation to solve
gyroEqn = PoissonBracketEquation.GyroEquation4D {
  speciesMass = kineticMass,
  speciesCharge = kineticCharge,
  bStarY = bStarYField,
}

pbSlvr = Updater.PoissonBracket4D {
   onGrid = grid_4d,
   -- basis functions to use
   basis = basis_4d,
   -- CFL number
   cfl = cfl,
   -- equation to solve
   equation = gyroEqn,
   -- let solver know about additional jacobian factor
   hasJacobian = true,
   jacobianField = jacobianField,
}

-- Hamiltonian
hamil = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Duplicate of hamil for reverting
hamilDup = hamil:duplicate()
-- Static hamiltonian (KE only)
hamilKE = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize KE part of hamil
initHamilKE = Updater.EvalOnNodes4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   evaluate = function (x,y,vPara,mu,t)
      return 0.5*kineticMass*vPara*vPara + math.abs(mu)*bFieldProfile(x)
   end
}
runUpdater(initHamilKE, 0.0, 0.0, {}, {hamilKE})

-- to store number density
numDensityKinetic = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
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

-- dynvector for total particle count
totalPtcl = DataStruct.DynVector { numComponents = 1, }
-- to compute total number of particles in domain
totalPtclCalc = Updater.IntegrateGeneralField2D {
   onGrid = grid_2d,
   basis = basis_2d,
}

-- to store the electrostatic potential on spatial grid
phi2d = DataStruct.Field2D {
   onGrid = grid_back_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
phi2dSmoothed = phi2d:duplicate()
phi2dBackground = phi2d:duplicate()
phi2dFluctuating = DataStruct.Field2D {
   onGrid = grid_fluct_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- to store electrostatic potential for addition to hamiltonian
phi4d = DataStruct.Field4D {
  onGrid = grid_4d,
  numComponents = basis_4d:numNodes(),
  ghost = {1, 1},
}

-- Updater to smooth out 2d field
smoothCalc = Updater.SimpleSmoothToC02D {
   onGrid = grid_2d,
   basis = basis_2d,
}

function applyBcToTotalPotential(fld)
  -- Subtract out fluctuating component from phi
  phi2dFluctuating:copy(fld)
  phi2dFluctuating:accumulate(-1.0, phi2dBackground)
  fld:accumulate(-1.0, phi2dFluctuating)
  -- Apply boundary conditions to fluctuating component
  applyBcToFluctuatingPotential(phi2dFluctuating)
  -- Add back to total field
  fld:accumulate(1.0, phi2dFluctuating)
end

function applyBcToFluctuatingPotential(fld)
  fld:sync()
end

function applyBcToBackgroundPotential(fld)
  fld:applyCopyBc(0, "lower")
  fld:applyCopyBc(0, "upper")
  fld:sync()
end

function applyBcToTotalDistF(fld)
  -- Subtract out fluctuating component from f
  fFluctuating:copy(fld)
  fFluctuating:accumulate(-1.0, fBackground)
  fld:accumulate(-1.0, fFluctuating)
  -- Apply boundary conditions to fluctuating component
  applyBcToFluctuatingDistF(fFluctuating)
  -- Add back to total field
  fld:accumulate(1.0, fFluctuating)
end

function applyBcToFluctuatingDistF(fld)
  fld:applyCopyBc(2, "lower")
  fld:applyCopyBc(2, "upper")
  fld:applyCopyBc(3, "lower")
  fld:applyCopyBc(3, "upper")
  fld:sync()
end

function applyBcToBackgroundDistF(fld)
  fld:applyCopyBc(0, "lower")
  fld:applyCopyBc(0, "upper")
  fld:applyCopyBc(2, "lower")
  fld:applyCopyBc(2, "upper")
  fld:applyCopyBc(3, "lower")
  fld:applyCopyBc(3, "upper")
  fld:sync()
end

function calcPotential(kineticN, adiabaticN, nAtCenter, phiOut)
  phiOut:copy(kineticN)
  phiOut:accumulate(-1.0, adiabaticN)
  phiOut:scale(adiabaticTemp/nAtCenter)
end

function calcHamiltonian(hamilKeIn, phi2dIn, hamilOut)
  -- accumulate static kinetic energy component
  hamilOut:copy(hamilKeIn)
  -- copy 2d potential to 4d field
  runUpdater(copy2dTo4d, 0.0, 0.0, {phi2dIn}, {phi4d})
  hamilOut:accumulate(kineticCharge, phi4d)
end

-- dynvector for field energy
fieldEnergy = DataStruct.DynVector { numComponents = 1, }
-- to compute field energy
fieldEnergyCalc = Updater.IntegrateGeneralField2D {
   onGrid = grid_2d,
   basis = basis_2d,
   moment = 2,
}

function calcNumDensity(fIn, nOut)
  runUpdater(numDensityCalc, 0.0, 0.0, {fIn, bField2d}, {nOut})
  nOut:scale(2*math.pi/kineticMass)
end

function calcDiagnostics(curr, dt)
  runUpdater(fieldEnergyCalc, curr, dt, {phi2dSmoothed}, {fieldEnergy})
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
  -- RK stage 1
  local myStatus, myDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {f, hamil}, {f1})

  if (myStatus == false) then
    return myStatus, myDtSuggested
  end

  applyBcToTotalDistF(f1)
  -- compute number density
  calcNumDensity(f1, numDensityKinetic)
  -- compute potential
  calcPotential(numDensityKinetic, numDensityAdiabatic, n0, phi2d)
  -- Apply boundary conditions
  applyBcToTotalPotential(phi2d)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi2d}, {phi2dSmoothed})
  -- Compute hamiltonian
  calcHamiltonian(hamilKE, phi2dSmoothed, hamil)

  -- RK stage 2
  local myStatus, myDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {f1, hamil}, {fNew})

  if (myStatus == false) then
    return myStatus, myDtSuggested
  end

  f1:combine(3.0/4.0, f, 1.0/4.0, fNew)
  applyBcToTotalDistF(f1)
  -- compute number density
  calcNumDensity(f1, numDensityKinetic)
  -- compute potential
  calcPotential(numDensityKinetic, numDensityAdiabatic, n0, phi2d)
  -- Apply boundary conditions
  applyBcToTotalPotential(phi2d)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi2d}, {phi2dSmoothed})
  -- Compute hamiltonian
  calcHamiltonian(hamilKE, phi2dSmoothed, hamil)

  -- RK stage 3
  local myStatus, myDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {f1, hamil}, {fNew})

  if (myStatus == false) then
    return myStatus, myDtSuggested
  end

  f1:combine(1.0/3.0, f, 2.0/3.0, fNew)
  applyBcToTotalDistF(f1)
  -- compute number density
  calcNumDensity(f1, numDensityKinetic)
  -- compute potential
  calcPotential(numDensityKinetic, numDensityAdiabatic, n0, phi2d)
  -- Apply boundary conditions
  applyBcToTotalPotential(phi2d)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi2d}, {phi2dSmoothed})
  -- Compute hamiltonian
  calcHamiltonian(hamilKE, phi2dSmoothed, hamil)

  f:copy(f1)

  return myStatus, myDtSuggested
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
  local step = 1
  local tCurr = tStart
  local myDt = initDt
  local status, dtSuggested

  while tCurr<=tEnd do
    -- Store fields that might need to be reverted
    fNewDup:copy(fNew)
    hamilDup:copy(hamil)

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
      fNew:copy(fNewDup)
      hamil:copy(hamilDup)
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

-- write data to H5 file
function writeFields(frameNum, tCurr)
   numDensityKinetic:write( string.format("n_%d.h5", frameNum), tCurr)
   fieldEnergy:write( string.format("fieldEnergy_%d.h5", frameNum), tCurr)
   phi2dSmoothed:write( string.format("phi_%d.h5", frameNum), tCurr)
end

applyBcToBackgroundDistF(f)

-- Compute initial kinetic density
calcNumDensity(f, numDensityKinetic)
-- Compute total number of particles
runUpdater(totalPtclCalc, 0.0, 0.0, {numDensityKinetic}, {totalPtcl})
-- Compute fraction to multiply distribution function by to get desired density
scaleFactor = deltaR*deltaR*n0/totalPtcl:lastInsertedData()
-- Scale fields and recalculate density
f:scale(scaleFactor)
Lucee.logInfo (string.format("-- Scaling distribution function by  %f", scaleFactor))
calcNumDensity(f, numDensityKinetic)
-- Store static numDensityAdiabatic field
numDensityAdiabatic:copy(numDensityKinetic)
-- Compute background phi
calcPotential(numDensityKinetic, numDensityAdiabatic, n0, phi2d)
-- Apply bouncary conditions
applyBcToBackgroundPotential(phi2d)
-- smooth potential
runUpdater(smoothCalc, 0.0, 0.0, {phi2d}, {phi2dSmoothed})

-- Store background phi
phi2dBackground:copy(phi2dSmoothed)
-- Store background f
fBackground:copy(f)

-- Scale f perturbation and add to f
f:accumulate(scaleFactor, fInitialPerturb)
-- Apply boundary conditions
applyBcToTotalDistF(f)
-- Compute potential
calcNumDensity(f, numDensityKinetic)
calcPotential(numDensityKinetic, numDensityAdiabatic, n0, phi2d)
applyBcToTotalPotential(phi2d)
runUpdater(smoothCalc, 0.0, 0.0, {phi2d}, {phi2dSmoothed})
calcHamiltonian(hamilKE, phi2dSmoothed, hamil)

-- Compute diagnostics for t = 0
calcDiagnostics(0.0, 0.0)
writeFields(0,0)

tCurr = tStart
fNew:copy(f)

for frame = 1, nFrames do
  Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
  dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
  tCurr = tCurr+tFrame
  writeFields(frame, tCurr)
  Lucee.logInfo ("")
end
