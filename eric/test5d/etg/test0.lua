-- Input file for ETG test problem
-- Species are referred to as the 'kinetic' or 'adiabatic' species

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.05
-- parameters to control time-stepping
tStart = 0.0
tEnd = 40e-6
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 100
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
q = 8
-- derived parameters
R          = R0 + 0.5*a
L_parallel = 2*math.pi*R*q
vtKinetic  = math.sqrt(kineticTemp*eV/kineticMass)
c_s        = math.sqrt(kineticTemp*eV/kineticMass)
omega_s    = math.abs(kineticCharge*B0/kineticMass)
rho_s      = c_s/omega_s
deltaR     = 32*rho_s
L_T        = R/6
ky_min     = 2*math.pi/deltaR
kz_min     = 2*math.pi/L_parallel
-- grid parameters: number of cells
N_X = 8
N_Y = 8
N_Z = 4
N_VPARA = 4
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = R
X_UPPER = R + deltaR
Y_LOWER = -deltaR/2
Y_UPPER = deltaR/2
Z_LOWER = 0
Z_UPPER = L_parallel
VPARA_UPPER = math.min(4, 2.5*math.sqrt(N_VPARA/4))*vtKinetic
VPARA_LOWER = -VPARA_UPPER
MU_LOWER = 0
MU_UPPER = math.min(16, 8*math.sqrt(N_MU/2))*kineticMass*vtKinetic*vtKinetic/(2*B0)

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
-- Additionally, grid_back_3d can probably be eliminated since it is 0
-- but I am keeping this in case it is not 0 somehow.

-- full 5d phase space grid
grid_5d = Grid.RectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
   periodicDirs = {0, 1, 2},
}
-- 3d spatial grid
grid_3d = Grid.RectCart3D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER},
   cells = {N_X, N_Y, N_Z},
   periodicDirs = {0, 1, 2},
}
grid_back_3d = Grid.RectCart3D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER},
   cells = {N_X, N_Y, N_Z},
   periodicDirs = {1,2},
}
-- create 5d basis functions
basis_5d = NodalFiniteElement5D.SerendipityElement {
   onGrid = grid_5d,
   polyOrder = polyOrder,
}
-- create 3d basis functions
basis_3d = NodalFiniteElement3D.SerendipityElement {
   onGrid = grid_3d,
   polyOrder = polyOrder,
}

-- distribution function for electrons
f = DataStruct.Field5D {
   onGrid = grid_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- for RK time-stepping
f1 = f:duplicate()
-- updated solution
fNew = f:duplicate()
-- for use in time-stepping
fNewDup = f:duplicate()
-- to store background distfElc
fBackground = f:duplicate()
fFluctuating = f:duplicate()
fInitialPerturb = fFluctuating:duplicate()

function bFieldProfile(x)
  return B0*R/x
end

function kineticTempProfile(x)
  return kineticTemp*(1 - (x-R)/L_T)
end

function perturbDensityProfile(x,y,z)
  local x0 = (X_LOWER+X_UPPER)/2
  local sigma = deltaR/4
  return 1e-2*(vtKinetic/omega_s)/L_T*math.cos(ky_min*y)*
    math.exp(-(x-x0)^2/(2*sigma^2))*(1 + math.cos(kz_min*z))
end

function fProfile(x,y,z,v,mu)
  return (2*math.pi*kineticTempProfile(x)*eV/kineticMass)^(-3/2)*
    math.exp(-kineticMass*v^2/(2*kineticTempProfile(x)*eV))*
    math.exp(-math.abs(mu)*bFieldProfile(x)/(kineticTempProfile(x)*eV))
end

-- initialize electron distribution function
initKineticF = Updater.EvalOnNodes5D {
   onGrid = grid_5d,
   basis = basis_5d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,v,mu,t)
		 return n0*fProfile(x,y,z,v,mu)
	 end
}
runUpdater(initKineticF, 0.0, 0.0, {}, {f})

-- initialize perturbation to electron distribution function
-- does not contain f_0, instead this will be multiplied with f_0
initKineticFPerturb = Updater.EvalOnNodes5D {
   onGrid = grid_5d,
   basis = basis_5d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,v,mu,t)
		 return 1 + perturbDensityProfile(x,y,z)
	 end
}
runUpdater(initKineticFPerturb, 0.0, 0.0, {}, {fInitialPerturb})

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

-- to copy phi to a 5d field
copy3dTo5d = Updater.NodalCopy3DTo5DFieldUpdater {
  -- 5D phase-space grid 
  onGrid = grid_5d,
  -- 5D phase-space basis functions
  basis5d = basis_5d,
  -- 3D spatial basis functions
  basis3d = basis_3d,
  -- Basis function order
  polyOrder = polyOrder,
}

-- Jacobian Factor (5D)
jacobianField = DataStruct.Field5D {
   onGrid = grid_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
initJacobian = Updater.EvalOnNodes5D {
  onGrid = grid_5d,
  basis = basis_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return kineticMass*kineticMass*bFieldProfile(x)
  end
}
-- Fill out jacobian
runUpdater(initJacobian, 0.0, 0.0, {}, {jacobianField})

-- B_y^* field (5D)
bStarYField = DataStruct.Field5D {
   onGrid = grid_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarY
initbStarY = Updater.EvalOnNodes5D {
  onGrid = grid_5d,
  basis = basis_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return -kineticMass*vPara/(kineticCharge*x)
  end
}
runUpdater(initbStarY, 0.0, 0.0, {}, {bStarYField})

-- B_z^* field (5D)
bStarZField = DataStruct.Field5D {
   onGrid = grid_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarZ
initbStarZ = Updater.EvalOnNodes5D {
  onGrid = grid_5d,
  basis = basis_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return bFieldProfile(x)
  end
}
runUpdater(initbStarZ, 0.0, 0.0, {}, {bStarZField})

-- define equation to solve
gyroEqn = PoissonBracketEquation.GyroEquation5D {
  speciesMass = kineticMass,
  speciesCharge = kineticCharge,
  bStarY = bStarYField,
  bStarZ = bStarZField,
}
pbSlvr = Updater.PoissonBracketOpt5D {
   onGrid = grid_5d,
   -- basis functions to use
   basis = basis_5d,
   -- CFL number
   cfl = cfl,
   -- equation to solve
   equation = gyroEqn,
   -- let solver know about additional jacobian factor
   jacobianField = jacobianField,
   updateDirections = {0,1,2,3},
   zeroFluxDirections = {3,4},
}

-- Hamiltonian
hamil = DataStruct.Field5D {
   onGrid = grid_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- Duplicate of hamil for reverting
hamilDup = hamil:duplicate()
-- Static hamiltonian (KE only)
hamilKE = hamil:duplicate()
-- Updater to initialize KE part of hamil
initHamilKE = Updater.EvalOnNodes5D {
   onGrid = grid_5d,
   basis = basis_5d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,vPara,mu,t)
      return 0.5*kineticMass*vPara*vPara + math.abs(mu)*bFieldProfile(x)
   end
}
runUpdater(initHamilKE, 0.0, 0.0, {}, {hamilKE})

-- to store number density
numDensityKinetic = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
numDensityAdiabatic = numDensityKinetic:duplicate()
-- to compute number density
numDensityCalc = Updater.DistFuncMomentCalcWeighted3D {
   -- 5D phase-space grid 
   onGrid = grid_5d,
   -- 5D phase-space basis functions
   basis5d = basis_5d,
   -- 5D spatial basis functions
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}

-- For temp calc
vParaSquaredKinetic = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
vParaSquaredCalc = Updater.DistFuncMomentCalcWeighted3D {
   -- 5D phase-space grid 
   onGrid = grid_5d,
   -- 5D phase-space basis functions
   basis5d = basis_5d,
   -- 3D spatial basis functions
   basis3d = basis_3d,
   -- desired moment (0, 1 or 2)
   moment = 2,
}

-- to store the electrostatic potential on spatial grid
phi3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
phi3dSmoothed = phi3d:duplicate()
phi3dBackground = DataStruct.Field3D {
   onGrid = grid_back_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
phi3dFluctuating = phi3d:duplicate()
-- to store electrostatic potential for addition to hamiltonian
phi5d = DataStruct.Field5D {
  onGrid = grid_5d,
  numComponents = basis_5d:numNodes(),
  ghost = {1, 1},
}

-- Updater to compute potential
phiCalc = Updater.ETGAdiabaticPotentialUpdater3D {
  onGrid = grid_3d,
  basis = basis_3d,
  adiabaticTemp = adiabaticTemp,
  adiabaticCharge = adiabaticCharge,
}

-- Updater to smooth out 3d field (phi)
smoothCalc = Updater.SimpleSmoothToC03D {
   onGrid = grid_back_3d,
   basis = basis_3d,
}

multiply5dCalc = Updater.FieldArithmeticUpdater5D {
  onGrid = grid_5d,
  basis = basis_5d,
  evaluate = function(fld1, fld2)
    return fld1*fld2
  end
}

function applyBcToTotalPotential(fld)
  -- Subtract out fluctuating component from phi
  phi3dFluctuating:copy(fld)
  phi3dFluctuating:accumulate(-1.0, phi3dBackground)
  fld:accumulate(-1.0, phi3dFluctuating)
  -- Apply boundary conditions to fluctuating component
  phi3dFluctuating:sync()
  -- Add back to total field
  fld:accumulate(1.0, phi3dFluctuating)
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
  -- Apply periodic boundary conditions to fluctuating component
  fFluctuating:sync()
  -- Add back to total field
  fld:accumulate(1.0, fFluctuating)
end

function calcHamiltonian(hamilKeIn, phi3dIn, hamilOut)
  -- accumulate static kinetic energy component
  hamilOut:copy(hamilKeIn)
  -- copy 3d potential to 5d field
  runUpdater(copy3dTo5d, 0.0, 0.0, {phi3dIn}, {phi5d})
  hamilOut:accumulate(kineticCharge, phi5d)
end

-- dynvector for field energy
fieldEnergy = DataStruct.DynVector { numComponents = 1, }
-- to compute field energy
fieldEnergyCalc = Updater.IntegrateGeneralField3D {
   onGrid = grid_3d,
   basis = basis_3d,
   moment = 2,
}

-- to scale distribution function
scaleInitDistF = Updater.ETGInitializeDensity5D {
  -- 5D phase-space grid 
   onGrid = grid_5d,
   -- 5D phase-space basis functions
   basis5d = basis_5d,
   -- 3D spatial basis functions
   basis3d = basis_3d,
   -- Desired constant density
   constantDensity = n0,
   polyOrder = polyOrder,
}

function calcNumDensity(fIn, nOut)
  runUpdater(numDensityCalc, 0.0, 0.0, {fIn, bField3d}, {nOut})
  nOut:scale(2*math.pi/kineticMass)
end

function calcDiagnostics(curr, dt)
  runUpdater(fieldEnergyCalc, curr, dt, {phi3dSmoothed}, {fieldEnergy})

  -- For use in determining temp
  runUpdater(vParaSquaredCalc, curr, dt, {f, bField3d}, {vParaSquaredKinetic})
  vParaSquaredKinetic:scale(2*math.pi/kineticMass)
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
  runUpdater(phiCalc, tCurr, myDt, {numDensityKinetic, numDensityAdiabatic}, {phi3d})
  -- Apply boundary conditions
  applyBcToTotalPotential(phi3d)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi3d}, {phi3dSmoothed})
  phi3dSmoothed:sync()
  -- Compute hamiltonian
  calcHamiltonian(hamilKE, phi3dSmoothed, hamil)

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
  runUpdater(phiCalc, tCurr, myDt, {numDensityKinetic, numDensityAdiabatic}, {phi3d})
  -- Apply boundary conditions
  applyBcToTotalPotential(phi3d)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi3d}, {phi3dSmoothed})
  phi3dSmoothed:sync()
  -- Compute hamiltonian
  calcHamiltonian(hamilKE, phi3dSmoothed, hamil)

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
  runUpdater(phiCalc, tCurr, myDt, {numDensityKinetic, numDensityAdiabatic}, {phi3d})
  -- Apply boundary conditions
  applyBcToTotalPotential(phi3d)
  -- smooth potential
  runUpdater(smoothCalc, tCurr, myDt, {phi3d}, {phi3dSmoothed})
  phi3dSmoothed:sync()
  -- Compute hamiltonian
  calcHamiltonian(hamilKE, phi3dSmoothed, hamil)

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
   phi3dSmoothed:write( string.format("phi_%d.h5", frameNum), tCurr)
   vParaSquaredKinetic:write( string.format("vParaSq_%d.h5", frameNum), tCurr)
   hamil:write( string.format("hamil_%d.h5", frameNum), tCurr)
   f:write( string.format("f_%d.h5", frameNum), tCurr)
   --phi3d:write( string.format("phiUnsmoothed_%d.h5", frameNum), tCurr)
end

-- Compute initial kinetic density
calcNumDensity(f, numDensityKinetic)
-- Scale distribution function and apply bcs
runUpdater(scaleInitDistF, 0.0, 0.0, {numDensityKinetic}, {f})
-- Recalculate number density
calcNumDensity(f, numDensityKinetic)
-- Store static numDensityAdiabatic field
numDensityAdiabatic:copy(numDensityKinetic)
-- Compute background phi
runUpdater(phiCalc, 0.0, 0.0, {numDensityKinetic, numDensityAdiabatic}, {phi3d})
-- Apply bouncary conditions
applyBcToBackgroundPotential(phi3d)
-- smooth potential
runUpdater(smoothCalc, 0.0, 0.0, {phi3d}, {phi3dSmoothed})

-- Store background phi
phi3dBackground:copy(phi3dSmoothed)
applyBcToBackgroundPotential(phi3dBackground)

-- Store background f
fBackground:copy(f)

-- Add perturbation to f
runUpdater(multiply5dCalc, 0.0, 0.0, {f, fInitialPerturb}, {fNew})
f:copy(fNew)
-- Apply boundary conditions
applyBcToTotalDistF(f)
-- Compute potential with perturbation added
calcNumDensity(f, numDensityKinetic)
runUpdater(phiCalc, 0.0, 0.0, {numDensityKinetic, numDensityAdiabatic}, {phi3d})
applyBcToTotalPotential(phi3d)
runUpdater(smoothCalc, 0.0, 0.0, {phi3d}, {phi3dSmoothed})
phi3dSmoothed:sync()
calcHamiltonian(hamilKE, phi3dSmoothed, hamil)

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
