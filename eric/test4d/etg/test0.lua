-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.1
-- parameters to control time-stepping
tStart = 0.0
tEnd = 1.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 1
tFrame = (tEnd-tStart)/nFrames -- time between frames

-- physical parameters
eV = Lucee.ElementaryCharge
kineticCharge = -Lucee.ElementaryCharge
ionCharge = Lucee.ElementaryCharge
kineticMass = Lucee.ElectronMass
adiabaticMass = 2*Lucee.ProtonMass
kineticTemp = 2072 -- [eV]
adiabaticTemp = 2072 -- [ev]
B0 = 1.91 -- [T]
R0 = 1.313 -- [m]
a = 0.4701 -- [m]
n0 = 4.992*10^(19) -- [1/m^3]
-- derived parameters
R = R0 + 0.5*a
vtKinetic = math.sqrt(kineticTemp*eV/kineticMass)
c_s = math.sqrt(kineticTemp*eV/kineticMass)
omega_s = math.abs(kineticCharge*B0/kineticMass)
rho_s = c_s/omega_s
deltaR = 32*rho_s
L_T = R/10
-- grid parameters: number of cells
N_X = 8
N_Y = 8
N_VPARA = 8
N_MU = 2
-- grid parameters: domain extent
X_LOWER = R
X_UPPER = R+deltaR
Y_LOWER = -deltaR/2
Y_UPPER = deltaR/2
VPARA_LOWER = -vtKinetic
VPARA_UPPER = vtKinetic
MU_LOWER = 0
MU_UPPER = 4*kineticMass*VPARA_UPPER*VPARA_UPPER/(2*B0)
-- full 4d phase space grid
grid_4d = Grid.RectCart4D {
   lower = {X_LOWER, Y_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_VPARA, N_MU},
}
-- 2d spatial grid
grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
}

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

-- create 
basis_4d = NodalFiniteElement4D.SerendipityElement {
   onGrid = grid_4d,
   polyOrder = polyOrder,
}

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
-- for RK time-stepping
f1 = f:duplicate()
-- updated solution
fNew = f:duplicate()
-- for use in time-stepping
fNewDup = fNew:duplicate()
-- to store background distfElc
backgroundf = f:duplicate()

function fProfile(x,y,v,mu)
  return n0*(2*math.pi*kineticTemp*eV/kineticMass)^(-3/2)*math.exp(-kineticMass*v^2/(2*kineticTemp*eV))*math.exp(-mu*bFieldProfile(x,y)/(kineticTemp*eV))
end

-- initialize electron distribution function
initKineticF = Updater.EvalOnNodes4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,v,mu,t)
		 return fProfile(x,y,v,mu)
	 end
}
runUpdater(initKineticF, 0.0, 0.0, {}, {f})

function bFieldProfile(x)
  return B0*R/x
end

-- Magnetic Field (2D)
bField2d = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize bfield2d
initbField2d = Updater.EvalOnNodes2D {
  onGrid = grid_4d,
  basis = basis_4d,
  shareCommonNodes = false,
  evaluate = function (x,y,t)
    return bFieldProfile(x)
  end
}
runUpdater(initbField2d, 0.0, 0.0, {}, {bField2d})

-- Jacobian Factor (4D)
jacobianField = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize Jacobian
initJacobian = Updater.EvalOnNodes4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   evaluate = function (x,y,vPara,mu,t)
      return kineticMass*kineticMass*bFieldProfile(x,y)
   end
}
runUpdater(initJacobian, 0.0, 0.0, {}, {jacobianField})

-- define equation to solve
advectionEqn = PoissonBracketEquation.AdvectionEquation4D {
}
gyroEqn = PoissonBracketEquation.GyroEquation4D {
  bParallelY = 
}

pbSlvr = Updater.PoissonBracket4D {
   onGrid = grid_4d,
   -- basis functions to use
   basis = basis_4d,
   -- CFL number
   cfl = cfl,
   -- equation to solve
   equation = advectionEqn,
}

-- Hamiltonian
hamil = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Static hamiltonian (KE only)
hamilKE = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}
-- Updater to initialize hamil
initHamilKE = Updater.EvalOnNodes4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   evaluate = function (x,y,vPara,mu,t)
      return 0.5*kineticMass*vPara*vPara + mu*bFieldProfile(x,y)
   end
}
runUpdater(initHamilKE, 0.0, 0.0, {}, {hamilKE})

-- to store number density
numDensityIon = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
numDensityElc = DataStruct.Field2D {
   onGrid = grid_2d,
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
-- to store the electrostatic potential on spatial grid
phi2d = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- to store electrostatic potential for addition to hamiltonian
phi4d = DataStruct.Field4d {
  onGrid = grid_2d,
  numComponents = basis_4d:numNodes(),
  ghost = {1, 1},
}

-- apply boundary conditions to a 4d field
function applyBc(fld)
  fld:applyPeriodicBc(0)
  fld:applyPeriodicBc(1)
  fld:applyCopyBc(2, "lower")
  fld:applyCopyBc(2, "upper")
  fld:applyCopyBc(3, "lower")
  fld:applyCopyBc(3, "upper")
end

function applyBcToFluctuatingPotential(fld)
  fld:applyPeriodicBc(0)
  fld:applyPeriodicBc(1)
end

function applyBcToBackgroundPotential(fld)
  fld:applyCopyBc(0, "lower")
  fld:applyCopyBc(0, "upper")
  fld:applyCopyBc(1, "lower")
  fld:applyCopyBc(1, "upper")
end

function applyBcToFluctuatingDistF(fld)
  fld:applyPeriodicBc(0)
  fld:applyPeriodicBc(1)
  fld:applyCopyBc(2, "lower")
  fld:applyCopyBc(2, "upper")
  fld:applyCopyBc(3, "lower")
  fld:applyCopyBc(3, "upper")
end

function applyBcToBackgroundDistF(fld)
  fld:applyCopyBc(0, "lower")
  fld:applyCopyBc(0, "upper")
  fld:applyCopyBc(1, "lower")
  fld:applyCopyBc(1, "upper")
  fld:applyCopyBc(2, "lower")
  fld:applyCopyBc(2, "upper")
  fld:applyCopyBc(3, "lower")
  fld:applyCopyBc(3, "upper")
end

function calcPotential()
end

function calcHamiltonian(hamilKeIn, phiIn, hamilOut)
  hamilOut:clear(0.0)
  -- copy potential to 4d field
  runUpdater(copy2dTo4d, 0.0, 0.0, {phiIn}, phi4d)
  hamilOut:accumulate(kineticCharge, phi4d)
  hamilOut:accumulate(1.0, hamilKeIn)
end

-- dynvector for field energy
fieldEnergy = DataStruct.DynVector { numComponents = 1, }
-- to compute field energy
fieldEnergyCalc = Updater.NormGrad2D {
   onGrid = grid_2d,
   basis = basis_2d,
}

function calcDiagnostics(curr, dt)
  runUpdater(fieldEnergyCalc, curr, dt, {phi2d}, {fieldEnergy})
  runUpdater(numDensityCalc, curr, dt, {f}, {numDensityElc})
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1
   local myStatus, myDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {f, hamil}, {f1})

   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   applyBc(f1)
   -- compute potential
   -- compute hamiltonian

   -- RK stage 2
   local myStatus, myDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {f1, hamil}, {fNew})

   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   f1:combine(3.0/4.0, f, 1.0/4.0, fNew)
   applyBc(f1)

   -- RK stage 3
   local myStatus, myDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {f1, hamil}, {fNew})

   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   f1:combine(1.0/3.0, f, 2.0/3.0, fNew)
   applyBc(f1)
   f:copy(f1)

   return myStatus, myDtSuggested
end

-- function to advance solution from tStart to tEnd
-- TODO: hamil duplicates, revert phi
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      fNewDup:copy(fNew)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	      myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	      -- time-step too large
	      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	      fNew:copy(fNewDup)
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
   f:write( string.format("f_%d.h5", frameNum), tCurr)
   numDensityElc:write( string.format("n_%d.h5", frameNum), tCurr)
   fieldEnergy:write( string.format("fieldEnergy_%d.h5", frameNum), tCurr)
end

applyBc(f)
fNew:copy(f)
-- compute phi
-- compute hamiltonian

-- write initial conditions
f:write("f_0.h5")

calcDiagnostics(0.0, 0.0)

tCurr = tStart
for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
