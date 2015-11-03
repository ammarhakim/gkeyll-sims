-- Input file for ETG test problem in 3x2v slab
-- Uses gyrokinetic poisson equation
-- Species are referred to as the 'kinetic' or 'adiabatic' species
-- 10-31-2015: 5d simulation with R/L_T = 4, ITG poisson equation solve
-- Box size = dR x 4*dR x 2*pi*R*q
-- Random initial conditions

-- phase-space decomposition
phaseDecomp = DecompRegionCalc5D.CartProd { cuts = {4, 8, 4, 1, 1} }
-- configuration space decomposition
confDecomp = DecompRegionCalc3D.SubCartProd5D {
   decomposition = phaseDecomp,
   collectDirections = {0, 1, 2},
}
-- 2d (x,z) configuration space decomposition
confDecomp2d = DecompRegionCalc2D.SubCartProd5D {
   decomposition = phaseDecomp,
   collectDirections = {0, 2},
}

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.05
-- parameters to control time-stepping
tStart = 0.0
tEnd = 10e-6
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 1000
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

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
L_T        = R/4
ky_min     = 2*math.pi/deltaR
kz_min     = 2*math.pi/L_parallel
-- grid parameters: number of cells
N_X = 8
N_Y = 32
N_Z = 8
N_VPARA = 4
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = R
X_UPPER = R + deltaR
Y_LOWER = -4*deltaR/2
Y_UPPER = 4*deltaR/2
Z_LOWER = 0
Z_UPPER = L_parallel
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

-- full 5d phase space grid
grid_5d = Grid.RectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
   periodicDirs = {0, 1, 2},
   decomposition = phaseDecomp,
}

grid_back_5d = Grid.RectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
   decomposition = phaseDecomp,
}
-- 3d spatial grid
grid_3d = Grid.RectCart3D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER},
   cells = {N_X, N_Y, N_Z},
   periodicDirs = {0, 1, 2},
   decomposition = confDecomp,
}
grid_back_3d = Grid.RectCart3D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER},
   cells = {N_X, N_Y, N_Z},
   decomposition = confDecomp,
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
   onGrid = grid_back_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- for RK time-stepping
f1 = DataStruct.Field5D {
   onGrid = grid_back_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- updated solution
fNew = DataStruct.Field5D {
   onGrid = grid_back_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- for use in time-stepping
fDup = DataStruct.Field5D {
   onGrid = grid_back_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- to store background distfElc
fBackground = DataStruct.Field5D {
   onGrid = grid_back_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
fFluctuating = DataStruct.Field5D {
   onGrid = grid_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}

function bFieldProfile(x)
  return B0*R/x
end

function kineticTempProfile(x)
  return kineticTemp*(1 - (x-R)/L_T)
end

-- set random seed based on processor number to avoid imprinting
-- domain decomp on ICs.
r = Lucee.getRank()
math.randomseed(100000*r+os.time())
function perturbDensityProfile(x,y,z)
  --local x0 = (X_LOWER+X_UPPER)/2
  --local sigma = deltaR/4
  --return 1e-3*(vtKinetic/omega_s)/L_T*math.cos(ky_min*y)*
  --  math.exp(-(x-x0)^2/(2*sigma^2))*(1 + math.cos(kz_min*z))
  return 1e-2*rho_s/L_T*math.random(-1,1)
end

function fProfile(x,y,z,v,mu)
  return (2*math.pi*kineticTempProfile(x)*eV/kineticMass)^(-3/2)*
    math.exp(-kineticMass*v^2/(2*kineticTempProfile(x)*eV))*
    math.exp(-math.abs(mu)*bFieldProfile(x)/(kineticTempProfile(x)*eV))
end

-- initialize electron distribution function
initKineticF = Updater.EvalOnNodes5D {
   onGrid = grid_back_5d,
   basis = basis_5d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,v,mu,t)
		 return n0*fProfile(x,y,z,v,mu)
	 end
}
runUpdater(initKineticF, 0.0, 0.0, {}, {f})
f:sync()

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
fInitialPerturb = DataStruct.Field5D {
   onGrid = grid_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- to copy phi to a 5d field
copy3dTo5d = Updater.CopyNodalFields3D_5D {
  onGrid = grid_5d,
  targetBasis = basis_5d,
  sourceBasis = basis_3d,
}
runUpdater(copy3dTo5d, 0.0, 0.0, {nInitialPerturbDG}, {fInitialPerturb})

-- Magnetic Field (3D)
bField3d = DataStruct.Field3D {
   onGrid = grid_back_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bField3d
initb = Updater.EvalOnNodes3D {
  onGrid = grid_back_3d,
  basis = basis_3d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,t)
    return bFieldProfile(x)
  end
}
runUpdater(initb, 0.0, 0.0, {}, {bField3d})
bField3d:sync()

-- Jacobian Factor (5D)
jacobianField = DataStruct.Field5D {
   onGrid = grid_back_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
initJacobian = Updater.EvalOnNodes5D {
  onGrid = grid_back_5d,
  basis = basis_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return kineticMass*kineticMass*bFieldProfile(x)
  end
}
-- Fill out jacobian
runUpdater(initJacobian, 0.0, 0.0, {}, {jacobianField})
jacobianField:sync()

-- B_y^* field (5D)
bStarYField = DataStruct.Field5D {
   onGrid = grid_back_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarY
initbStarY = Updater.EvalOnNodes5D {
  onGrid = grid_back_5d,
  basis = basis_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return -kineticMass*vPara/(kineticCharge*x)
  end
}
runUpdater(initbStarY, 0.0, 0.0, {}, {bStarYField})
bStarYField:sync()

-- B_z^* field (5D)
bStarZField = DataStruct.Field5D {
   onGrid = grid_back_5d,
   numComponents = basis_5d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bStarZ
initbStarZ = Updater.EvalOnNodes5D {
  onGrid = grid_back_5d,
  basis = basis_5d,
  shareCommonNodes = false,
  evaluate = function (x,y,z,vPara,mu,t)
    return bFieldProfile(x)
  end
}
runUpdater(initbStarZ, 0.0, 0.0, {}, {bStarZField})
bStarZField:sync()

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
   onGrid = grid_back_5d,
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
hamilKE:sync()

-- to store number density
numDensityKinetic = DataStruct.Field3D {
   onGrid = grid_back_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
numDensityAdiabatic = numDensityKinetic:duplicate()
numDensityDelta = numDensityKinetic:duplicate()
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

-- to store the (fluctuating) electrostatic potential on spatial grid
phi3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to store electrostatic potential for addition to hamiltonian
phi5d = DataStruct.Field5D {
  onGrid = grid_5d,
  numComponents = basis_5d:numNodes(),
  ghost = {1, 1},
}

multiply5dCalc = Updater.FieldArithmeticUpdater5D {
  onGrid = grid_5d,
  basis = basis_5d,
  evaluate = function(fld1, fld2)
    return fld1*fld2
  end
}

grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Z_LOWER},
   upper = {X_UPPER, Z_UPPER},
   cells = {N_X, N_Z},
   decomposition = confDecomp2d,
}
-- create 2d basis functions
basis_2d = NodalFiniteElement2D.SerendipityElement {
   onGrid = grid_2d,
   polyOrder = polyOrder,
}
-- For zonal average
zonalAverageCalc = Updater.ZonalAverageCalc3D {
  onGrid = grid_3d,
  basis3d = basis_3d,
  basis2d = basis_2d,
}
-- to store the zonal-averaged potential in 2d
zonalAveragedPhi2d = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- to store the zonal-averaged potential in 3d
zonalAveragedPhi3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to store the zonal-averaged potential poisson solve result in 3d
zonalAveragedPhiResult3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to copy the (x,z) field back to (x,y,z)
copy2DTo3D = Updater.CopyNodalFields2D_3D {
   onGrid = grid_3d,
   sourceBasis = basis_2d,
   targetBasis = basis_3d,
   coordinateMap = {0,2},
}

function applyBcToTotalDistF(fld)
  -- Subtract out fluctuating component from f
  fFluctuating:copy(fld)
  fFluctuating:accumulate(-1.0, fBackground)
  fld:accumulate(-1.0, fFluctuating)
  -- Apply periodic boundary conditions to fluctuating component
  fFluctuating:sync()
  -- Add back to total field
  fld:accumulate(1.0, fFluctuating)
  fld:sync()
end

function calcHamiltonian(hamilKeIn, phi3dIn, hamilOut)
  -- accumulate static kinetic energy component
  hamilOut:copy(hamilKeIn)
  -- copy 3d potential to 5d field
  runUpdater(copy3dTo5d, 0.0, 0.0, {phi3dIn}, {phi5d})
  phi5d:sync()
  hamilOut:accumulate(kineticCharge, phi5d)
  hamilOut:sync()
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
  runUpdater(fieldEnergyCalc, curr, dt, {phi3d}, {fieldEnergy})
  -- Calc perturbed density
  numDensityDelta:copy(numDensityKinetic)
  numDensityDelta:accumulate(-1.0, numDensityAdiabatic)
end

-- create updater to solve Poisson equation
poissonSlvr = Updater.FemPoisson3D {
   onGrid = grid_3d,
   basis = basis_3d,
   periodicDirs = {0, 1, 2},
   sourceNodesShared = false, -- default true
   solutionNodesShared = false, -- default true
   writeStiffnessMatrix = false,
   modifierConstant = -kineticTemp/(adiabaticTemp*rho_s^2),
   isGyroKineticPoisson = true, -- solve 3d GK perp flag
}
-- create updater to solve Poisson equation for flux-surface averaged potential
fluxSurfaceAveragePotentialSlvr = Updater.FemPoisson3D {
   onGrid = grid_3d,
   basis = basis_3d,
   periodicDirs = {0, 1, 2},
   sourceNodesShared = false, -- default true
   solutionNodesShared = false, -- default true
   writeStiffnessMatrix = false,
   isGyroKineticPoisson = true, -- solve 3d GK perp flag
}
-- Stores right-hand side of poisson equation
poissonSource = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

function calcPotential(nKineticIn, nAdiabaticIn, phiOut)
  poissonSource:copy(nAdiabaticIn)
  poissonSource:accumulate(-1.0, nKineticIn)
  poissonSource:scale(-kineticTemp/(n0*rho_s^2)) -- kineticTemp is in eV

  -- compute flux-surface average of poissonSource
  --runUpdater(zonalAverageCalc, 0.0, 0.0, {poissonSource}, {zonalAveragedPhi2d})
  --zonalAveragedPhi2d:scale(1/(Y_UPPER-Y_LOWER))
  -- copy 2d result back to a 3d field for use in poisson solve
  --runUpdater(copy2DTo3D, 0.0, 0.0, {zonalAveragedPhi2d}, {zonalAveragedPhi3d})
  -- solve for flux-surface averaged potential
  --runUpdater(fluxSurfaceAveragePotentialSlvr, 0.0, 0.0, {zonalAveragedPhi3d}, {zonalAveragedPhiResult3d})
  -- accumulate flux-surface averaged potential to poissonSource to solve for phi
  --poissonSource:accumulate(1/(rho_s^2), zonalAveragedPhiResult3d)

  runUpdater(poissonSlvr, 0.0, 0.0, {poissonSource}, {phiOut})
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
  -- compute potential and apply boundary conditions
  calcPotential(numDensityKinetic, numDensityAdiabatic, phi3d)
  phi3d:sync()
  -- Compute hamiltonian
  calcHamiltonian(hamilKE, phi3d, hamil)

  -- RK stage 2
  local myStatus, myDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {f1, hamil}, {fNew})

  if (myStatus == false) then
    return myStatus, myDtSuggested
  end

  f1:combine(3.0/4.0, f, 1.0/4.0, fNew)
  applyBcToTotalDistF(f1)
  -- compute number density
  calcNumDensity(f1, numDensityKinetic)
  -- compute potential and apply boundary conditions
  calcPotential(numDensityKinetic, numDensityAdiabatic, phi3d)
  phi3d:sync()
  -- Compute hamiltonian
  calcHamiltonian(hamilKE, phi3d, hamil)

  -- RK stage 3
  local myStatus, myDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {f1, hamil}, {fNew})

  if (myStatus == false) then
    return myStatus, myDtSuggested
  end

  f1:combine(1.0/3.0, f, 2.0/3.0, fNew)
  applyBcToTotalDistF(f1)
  -- compute number density
  calcNumDensity(f1, numDensityKinetic)
  -- compute potential and apply boundary conditions
  calcPotential(numDensityKinetic, numDensityAdiabatic, phi3d)
  phi3d:sync()
  -- Compute hamiltonian
  calcHamiltonian(hamilKE, phi3d, hamil)

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
    fDup:copy(fNew)
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
      fNew:copy(fDup)
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
   --numDensityKinetic:write( string.format("n_%d.h5", frameNum), tCurr)
   fieldEnergy:write( string.format("fieldEnergy_%d.h5", frameNum), tCurr)
   numDensityDelta:write( string.format("nDelta_%d.h5", frameNum), tCurr)
   --phi3dSmoothed:write( string.format("phi_%d.h5", frameNum), tCurr)
   --vParaSquaredKinetic:write( string.format("vParaSq_%d.h5", frameNum), tCurr)
   --hamil:write( string.format("hamil_%d.h5", frameNum), tCurr)
   --f:write( string.format("f_%d.h5", frameNum), tCurr)
   phi3d:write( string.format("phi_%d.h5", frameNum), tCurr)
   f:write( string.format("f_%d.h5", frameNum), tCurr)
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
-- Store background f
fBackground:copy(f)

startFrame = 1

if Lucee.IsRestarting then
  f:read("f_" .. Lucee.RestartFrame .. ".h5")
  f:sync()

  startFrame = Lucee.RestartFrame + 1
  tCurr = tStart + tFrame*Lucee.RestartFrame
else
  -- Add perturbation to f
  runUpdater(multiply5dCalc, 0.0, 0.0, {f, fInitialPerturb}, {fNew})
  --fNew:sync()
  f:copy(fNew)
  -- Apply boundary conditions
  applyBcToTotalDistF(f)
end

-- Compute potential with perturbation added
calcNumDensity(f, numDensityKinetic)
calcPotential(numDensityKinetic, numDensityAdiabatic, phi3d)
phi3d:sync()
calcHamiltonian(hamilKE, phi3d, hamil)

-- Compute diagnostics for t = 0
calcDiagnostics(tCurr, 0.0)
writeFields(startFrame-1,tCurr)

tCurr = tStart

for frame = 1, nFrames do
  Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
  dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
  tCurr = tCurr+tFrame
  writeFields(frame, tCurr)
  Lucee.logInfo ("")
end
