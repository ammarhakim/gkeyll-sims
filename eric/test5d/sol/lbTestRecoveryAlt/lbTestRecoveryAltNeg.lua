-- Input file for a 3X2V Simulation of LAPD
-- mpiexec -n 1 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i lbTestRecoveryAltNeg.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist

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
dtSuggested = 1e-7 -- initial time-step to use (will be adjusted)
nFrames = 1
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

-- physical parameters
eV        = Lucee.ElementaryCharge
elcCharge = -eV
ionCharge = eV
eps0      = Lucee.Epsilon0 -- permittivity of free space
ionMass   = 3.973*Lucee.ProtonMass 
elcMass   = ionMass/400 -- REDUCED ELECTRON MASS Lucee.ElectronMass
elcTemp   = 6 -- [eV]
ionTemp   = 1 -- [eV]
B0 = 0.0398  -- [T]
n0 = 2*10^18 -- [1/m^3]
V_bias = 0 -- [V]
-- derived parameters
vtElc   = math.sqrt(0.5*elcTemp*eV/elcMass) -- only used for grid
vtIon   = math.sqrt(ionTemp*eV/ionMass)
omega_i = math.abs(ionCharge*B0/ionMass)

c_s = math.sqrt(elcTemp*eV/ionMass)
rho_s = c_s/omega_i
epsilonDensity = 1e-5 -- density floor to enforce

deltaX  = 100*rho_s
deltaY  = 100*rho_s
R          = 40*rho_s -- [m]
L_parallel = 36*R -- [m]
-- grid parameters: number of cells
N_X = 1
N_Y = 1
N_Z = 1
N_VPARA = 10
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = -deltaX/2
X_UPPER = deltaX/2
Y_LOWER = -deltaY/2
Y_UPPER = deltaY/2
Z_LOWER = -L_parallel/2
Z_UPPER = L_parallel/2
-- source parameters
L_s = 0.5*rho_s
r_s = 20*rho_s
sourceAmplitude = 0.03*n0*c_s/R
-- total volume
totalVol = (X_UPPER-X_LOWER)*(Y_UPPER-Y_LOWER)*(Z_UPPER-Z_LOWER)

VPARA_UPPER_ELC = 4*vtElc
VPARA_LOWER_ELC = -VPARA_UPPER_ELC
MU_LOWER_ELC = 0
MU_UPPER_ELC = 0.75*elcMass*(VPARA_UPPER_ELC*VPARA_UPPER_ELC)/(2*B0)

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
grid_elc_5d = Grid.NonUniformRectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER_ELC, 0},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER_ELC, MU_UPPER_ELC},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
   mappings = {
     function(zeta)
       return zeta
     end,
     function(zeta)
       return zeta
     end,
     function(zeta)
       return zeta
     end,
     function(zeta)
       return zeta
     end,
     function(zeta)
       return zeta
     end,
   },
   decomposition = phaseDecomp,
}

-- 3d spatial grid
grid_3d = Grid.RectCart3D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER},
   cells = {N_X, N_Y, N_Z},
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
   num1DGaussPoints = 2,
}
-- create 3d basis functions
basis_3d = NodalFiniteElement3D.SerendipityElement {
   onGrid = grid_3d,
   polyOrder = polyOrder,
   num1DGaussPoints = 2,
}
-- create 2d basis functions for (x,y) grid
basis_2d = NodalFiniteElement2D.SerendipityElement {
   onGrid = grid_2d,
   polyOrder = polyOrder,
   num1DGaussPoints = 2,
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

-- Magnetic field
function bFieldProfile(x)
  return B0
end

hamilKeElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
hamilDerivKeElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
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

function initBackgroundProfile(x,y,z)
  return 1
end

function elcTempProfile(x,y,z)
  return elcTemp/3
end

function fElcProfile(x,y,z,v,mu)
  local u = 0.2*vtElc
  return math.exp(-elcMass*(v-u)^2/(2*elcTempProfile(x,y,z)*eV))*
    math.exp(-math.abs(mu)*bFieldProfile(x)/(elcTempProfile(x,y,z)*eV))
    --*(2*math.pi*elcTempProfile(x,y,z)*eV/elcMass)^(-3/2)
end
-- to copy phi to the 5d space the electron hamiltonian is stored on
copy3dTo5dElc = Updater.CopyNodalFields3D_5D {
  -- 5D phase-space grid 
  onGrid = grid_elc_5d,
  -- 5D phase-space basis functions
  targetBasis = basis_elc_5d,
  -- 3D spatial basis functions
  sourceBasis = basis_3d,
}
-- to multiply two 5d fields together
multiply5dCalcElc = Updater.FieldArithmeticUpdater5D {
  onGrid = grid_elc_5d,
  basis = basis_elc_5d,
  evaluate = function(fld1, fld2)
    return fld1*fld2
  end
}

-- initialize electron distribution function
initfElc = Updater.EvalOnNodes5D {
   onGrid = grid_elc_5d,
   basis = basis_elc_5d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,v,mu,t)
		 --return n0*fElcProfile(x,y,z,v,mu)
		 local targetTemp = 1*eV--elcTemp/2*eV
		 local v_m = math.sqrt(targetTemp*3/elcMass)
		 local uTar = 0.25*math.sqrt(targetTemp/elcMass)
		 local mu_m = 2*targetTemp/B0
		 --return math.cos(math.pi/(2*VPARA_UPPER_ELC)*v)
		 if (v <= uTar + v_m and v >= uTar - v_m) and math.abs(mu) <= mu_m then
		   return 1
     else
       return 0
     end
	 end
}
runUpdater(initfElc, 0.0, 0.0, {}, {fElc})
fElc:sync()

-- to scale initial condition
scaleInitDistfElc = Updater.SOLInitializeDensity {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/elcMass,
  evaluate = function(x,y,z,t)
    return n0
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
  evaluate = function(x,y,z,t)
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
  evaluate = function(x,y,z,vPara,mu,t)
    return elcMass*elcMass*bFieldProfile(x)
  end
}
-- Fill out jacobian
runUpdater(initJacobianElc, 0.0, 0.0, {}, {jacobianFieldElc})
jacobianFieldElc:sync()

-- B_y^* field for electrons (5D)
bStarYFieldElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
bStarYFieldElc:clear(0.0)
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
-- electron total tempearture
temperatureElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron parallel tempearture
temperatureParaElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- electron perp tempearture
temperaturePerpElc = DataStruct.Field3D {
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
kineticEnergyElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
totalEnergyElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

generalIntegralElcCalc = Updater.SOLGeneralIntegralAtNodeCalc {--SOLWeightedProjectionCalc {
   onGrid = grid_elc_5d,
   basis5d = basis_elc_5d,
   basis3d = basis_3d,
   basis2d = basis_2d,
   scaleFactor = 2*math.pi/elcMass,
}
-- to compute temperature for electrons
temperatureElcCalc = Updater.SOLTemperatureAtNodeCalc {
   onGrid = grid_3d,
   basis = basis_3d,
   speciesMass = elcMass,
}
-- to calculate the average of a 3d field
fieldInt3d = Updater.IntegrateGeneralField3D {
  onGrid = grid_3d,
  basis = basis_3d,
}
-- stores total ion/elc density. should only be computed once per time step
totalElcCounter = DataStruct.DynVector { numComponents = 1, }
totalMomentumCounter = DataStruct.DynVector { numComponents = 1, }
totalEnergyCounter = DataStruct.DynVector { numComponents = 1, }
-- stores total ion/elc energy. should only be computed once per time step
totalEnergyElcVec = DataStruct.DynVector { numComponents = 1, }

energyAtNodeAfterPositivity = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
paraEnergyAtNodeBeforePositivity = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
perpEnergyAtNodeBeforePositivity = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
energyAtNodeAfterMaxDrag = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

elcSourceTempHistory = DataStruct.DynVector { numComponents = 1, }

-- lenard berstein code
function getElcAlpha(t)
  -- see http://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node35.html
  --local n_e0 = avgElcDensity:lastInsertedData()/totalVol
  local logLambda = 6.6-0.5*math.log(n0/10^(20)) + 1.5*math.log(elcTemp)
  -- alpha needs to be multiplied by average n/(kTe)^(3/2) in C++ code
  -- enhanced by artificial factor of 10
  return logLambda*elcCharge^4/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(elcMass))/10
end

-- lenard berstein code
function getIonAlpha(t)
  -- see http://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node35.html
  --local n_i0 = avgElcDensity:lastInsertedData()/totalVol
  local logLambda = 6.6-0.5*math.log(n0/10^(20)) + 1.5*math.log(elcTemp)
  -- alpha needs to be multiplied by average n/(kTe)^(3/2) in C++ code
  return logLambda*elcCharge^4/(12*math.pi^(3/2)*eps0^2*math.sqrt(ionMass))
end

diffSlvrElc = Updater.LenardBernsteinDiff3D2VUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  cfl = cfl,
  onlyIncrement = true,
  speciesMass = elcMass,
  alpha = function(t)
    return getElcAlpha(t)
  end
}

dragSlvrElc = Updater.LenardBernsteinDrag3D2VUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  cfl = cfl,
  onlyIncrement = true,
  speciesMass = elcMass,
  alpha = function(t)
    return getElcAlpha(t)
  end
}

lenardBernsteinScaleElcUpdater = Updater.SOLLenardBernsteinScaleCell5D {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  cfl = cfl,
  speciesMass = elcMass,
  alpha = function(t)
    return getElcAlpha(t)
  end
}

positivityScaleElcUpdater = Updater.SOLPositivityScaleNodeUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
}

fElcCollisions = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
elcEnergyAtNodeDrag = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
elcEnergyAtNodeDiff = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
weightedFieldInt3d = Updater.SOLTotalIntegralCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
}

oneField3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
oneField3d:clear(1.0)

phiCurrentProductAtEdgeVec = DataStruct.DynVector { numComponents = 2, }
currentIntegratedAtEdgeVec = DataStruct.DynVector { numComponents = 2, }
beforePosEnergyVec = DataStruct.DynVector { numComponents = 1, }
afterPosEnergyVec = DataStruct.DynVector { numComponents = 1, }
beforePosMomentumVec = DataStruct.DynVector { numComponents = 1, }
afterPosMomentumVec = DataStruct.DynVector { numComponents = 1, }
beforeFloorIonEnergyVec = DataStruct.DynVector { numComponents = 1, }
beforeFloorElcEnergyVec = DataStruct.DynVector { numComponents = 1, }
afterFloorIonEnergyVec = DataStruct.DynVector { numComponents = 1, }
afterFloorElcEnergyVec = DataStruct.DynVector { numComponents = 1, }
beforeFloorElcNumberVec = DataStruct.DynVector { numComponents = 1, }
afterFloorElcNumberVec = DataStruct.DynVector { numComponents = 1, }
avgIonDensity = DataStruct.DynVector { numComponents = 1, }

-- to calculate elc energy integrated over a cell
elcEnergyAtCellCalc = Updater.SOLEnergyAtCellCalc {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  scaleFactor = 2*math.pi/elcMass^3
}
-- to get correct energy in every cell
positivityScaleCellElcUpdater = Updater.SOLPositivityScaleCellUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
}

setDesiredDensityElc = Updater.SOLDesiredDensity5D {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/elcMass,
}
kineticEnergyBeforeElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

kineticEnergyAfterElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

kineticEnergyStage1Elc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
kineticEnergyStage2Elc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
kineticEnergyStage3Elc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

function calcLenardBernsteinElcElc(tCurr, myDt, fInitial, fFinal)
  -- calculate inputs for collision operator
  -- compute electron moments
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fInitial, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fInitial, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fInitial, hamilKeElc, bField3d}, {kineticEnergyElc})
  -- compute electron temperature
  local negativeTemperature = runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc,
   kineticEnergyElc, mom1dir3Elc}, {temperatureElc})
  
  if (negativeTemperature == false) then
    endExecution = true
    Lucee.logInfo(string.format("-- calcLenardBernsteinElcElc: Ending execution because negative temperature. --"))
    return false, 0.0
  end

  local dragStatus, dragDtSuggested = runUpdater(dragSlvrElc, tCurr, myDt, {fInitial, mom1dir3Elc,
    numDensityElc, temperatureElc, numDensityElc}, {fElcCollisions})
  -- compute energy at each node in the cell for drag term only
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcCollisions, jacobianFieldElc, hamilKeElc}, {elcEnergyAtNodeDrag})
  -- accumulate drag term times dt to distribution function
  fFinal:accumulate(myDt, fElcCollisions)
  
  runUpdater(diffSlvrAltElc, tCurr, myDt, {fInitial, temperatureElc, temperatureElc, temperatureElc,
    --temperatureParaEqElc, temperaturePerpEqElc,
    bField3d, numDensityElc}, {fElcCollisions})

  -- compute energy integrated in a cell for diffusion term only
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcCollisions, jacobianFieldElc, hamilKeElc},
    {elcEnergyAtNodeDiff})

  -- accumulate correct fraction of diffusion term to distribution function
  local diffStatus, diffDtSuggested = runUpdater(lenardBernsteinScaleElcUpdater, tCurr, myDt, {fElcCollisions,
    elcEnergyAtNodeDrag, elcEnergyAtNodeDiff, temperatureElc, temperatureElc, temperatureElc,
    --temperatureParaEqElc, temperaturePerpEqElc,
    bField3d, numDensityElc}, {fFinal})

  if (diffStatus == false or dragStatus == false) then
    return false, math.min(dragDtSuggested, diffDtSuggested)
  end

  return true, math.min(dragDtSuggested, diffDtSuggested)
end

endExecution = false

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
  -- RK stage 1
  f1Elc:copy(fElc)
  local myStatus, myDtSuggested = calcLenardBernsteinElcElc(tCurr, myDt, fElc, f1Elc)
  if (myStatus == false) then
    return false, myDtSuggested
  end

  f1Elc:sync()

  -- RK stage 2
  fNewElc:copy(f1Elc)
  myStatus, myDtSuggested = calcLenardBernsteinElcElc(tCurr, myDt, f1Elc, fNewElc)
  if (myStatus == false) then
    return false, myDtSuggested
  end
  f1Elc:combine(3.0/4.0, fElc, 1.0/4.0, fNewElc)

  f1Elc:sync()

  -- RK stage 3
  fNewElc:copy(f1Elc)
  myStatus, myDtSuggested = calcLenardBernsteinElcElc(tCurr, myDt, f1Elc, fNewElc)

  if (myStatus == false) then
    return false, myDtSuggested
  end
  f1Elc:combine(1.0/3.0, fElc, 2.0/3.0, fNewElc)

  fElc:copy(f1Elc)

  -- calculate continuous time diagnostics
  calcContinuousDiagnostics(tCurr, myDt)

  return true, myDtSuggested
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
  local step = 1
  local tCurr = tStart
  local myDt = initDt
  local status, dtSuggested

  local startTime = os.clock()

  while tCurr<=tEnd do
    -- Store fields that might need to be reverted
    fDupElc:copy(fElc)

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
      end
    else
      calcDiagnostics(tCurr, myDt)
      tCurr = tCurr + myDt
      -- commented next line to maintain a uniform time step
      --myDt = dtSuggested
      step = step + 1

      if (tCurr >= tEnd) then
        break
      end
    end
  end

  totalTime = totalTime + os.clock()-startTime
  totalSteps = totalSteps + step

  return dtSuggested, tCurr
end

-- function containing diagnostics that must be called at the end of every frame
function calcDiagnostics(tCurr, myDt)
  --runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, oneFieldElc, bField3d}, {numDensityElc})
  --runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {totalElcCounter})
  --Lucee.logInfo(string.format("Avg number = %g",totalElcCounter:lastInsertedData()/totalVol))
end

-- function containing diagnostics that must be called at the end of every time step
function calcContinuousDiagnostics(tCurr, myDt)
  -- compute total guiding center density by integrating over whole domain
  -- note: includes a volume factor that should be removed in post-processing
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {totalElcCounter})

  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(fieldInt3d, tCurr, myDt, {mom1dir3Elc}, {totalMomentumCounter})

  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilKeElc, bField3d}, {kineticEnergyElc})
  runUpdater(fieldInt3d, tCurr, myDt, {kineticEnergyElc}, {totalEnergyCounter})

end

-- use to remember original density before positivity modification
numDensityTarget = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

positivityDragSliceMuElc = Updater.SOLPositivityDragSliceInMuUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
}

positivityDragSliceVParElc = Updater.SOLPositivityDragSliceUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
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

function applyElcPositivity(tCurr, myDt, distfElc)
  -- compute energy before positivity
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilKeElc}, {kineticEnergyBeforeElc})
  --runUpdater(generalIntegralElcCalc, tCurr, myDt, {distfElc, hamilKeElc, bField3d}, {kineticEnergyBeforeElc})
  runUpdater(fieldInt3d, tCurr, myDt, {kineticEnergyBeforeElc}, {beforePosEnergyVec})
  -- compute parallel energy before positivity
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilParaElc}, {paraEnergyAtNodeBeforePositivity})
  -- compute perp energy before positivity
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilPerpElc}, {perpEnergyAtNodeBeforePositivity})
  -- compute density before positivity
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {distfElc, oneFieldElc, bField3d}, {numDensityElc})
  -- store numDensityElc so we can restore it by scaling
  numDensityTarget:copy(numDensityElc)
  -- compute momentum before positivity
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {distfElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(fieldInt3d, tCurr, myDt, {mom1dir3Elc}, {beforePosMomentumVec})
  
  -- apply positivity to elcs
  distfElc:clearNegative()
  -- restore original density by scaling
  runUpdater(setDesiredDensityElc, tCurr, myDt, {numDensityTarget, bField3d}, {distfElc})

  -- need to sync ghost cells for drag update
  distfElc:sync()
  --- compute elc perp energy at cell after positivity is applied
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilPerpElc},
    {energyAtNodeAfterPositivity})
  -- take a maximum drag step in mu direction
  runUpdater(positivityDragSliceMuElc, tCurr, myDt, {distfElc, energyAtNodeAfterPositivity,
    perpEnergyAtNodeBeforePositivity}, {fElcPositivityDrag})
  -- compute energy at cell after maximum drag step
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcPositivityDrag, jacobianFieldElc, hamilPerpElc},
    {energyAtNodeAfterMaxDrag})
  -- Compute the difference between distf post-positivity and distf post-drag
  fElcPositivityDelta:combine(1.0, fElcPositivityDrag, -1.0, distfElc)
  -- add correct amount of fElcPositivityDelta to get correct energy at each cell
  runUpdater(positivityScaleCellElcUpdater, tCurr, myDt, {fElcPositivityDelta, perpEnergyAtNodeBeforePositivity,
    energyAtNodeAfterPositivity, energyAtNodeAfterMaxDrag}, {distfElc})

  -- need to sync ghost nodes for drag update
  distfElc:sync()
  -- compute elc parallel energy at nodes after positivity is applied
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilParaElc},
    {energyAtNodeAfterPositivity})
  -- take a maximum drag step in v-parallel direction
  runUpdater(positivityDragSliceVParElc, tCurr, myDt, {distfElc, energyAtNodeAfterPositivity,
    paraEnergyAtNodeBeforePositivity, mom1dir3Elc, numDensityElc}, {fElcPositivityDrag})
  -- compute energy in node after maximum drag step
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcPositivityDrag, jacobianFieldElc, hamilParaElc},
    {energyAtNodeAfterMaxDrag})
  -- Compute the difference between distf post-positivity and distf post-drag
  fElcPositivityDelta:combine(1.0, fElcPositivityDrag, -1.0, distfElc)
  -- add correct amount of fElcPositivityDelta to get correct energy at each node
  runUpdater(positivityScaleCellElcUpdater, tCurr, myDt, {fElcPositivityDelta, paraEnergyAtNodeBeforePositivity,
    energyAtNodeAfterPositivity, energyAtNodeAfterMaxDrag}, {distfElc})
  distfElc:sync()

  -- compute energy after positivity
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {distfElc, jacobianFieldElc, hamilKeElc}, {kineticEnergyAfterElc})
  --runUpdater(generalIntegralElcCalc, tCurr, myDt, {distfElc, hamilKeElc, bField3d}, {kineticEnergyAfterElc})
  runUpdater(fieldInt3d, tCurr, myDt, {kineticEnergyAfterElc}, {afterPosEnergyVec})
  -- compute momentum after positivity
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {distfElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(fieldInt3d, tCurr, myDt, {mom1dir3Elc}, {afterPosMomentumVec})

  -- compute increment in momentum
  local deltaU = (afterPosMomentumVec:lastInsertedData() - beforePosMomentumVec:lastInsertedData())/beforePosMomentumVec:lastInsertedData()
  Lucee.logInfo(string.format("Electron Momentum change = %g %%",100*deltaU  ))

  return (afterPosEnergyVec:lastInsertedData() - beforePosEnergyVec:lastInsertedData())--/myDt
end

-- write data to H5 file
function writeFields(frameNum, tCurr)
  fElc:write( string.format("fElc_%d.h5", frameNum), tCurr)
  totalElcCounter:write( string.format("totalNumberCounter_%d.h5", frameNum), tCurr)
  totalMomentumCounter:write( string.format("totalMomentumCounter_%d.h5", frameNum), tCurr)
  totalEnergyCounter:write( string.format("totalEnergyCounter_%d.h5", frameNum), tCurr)
end

runUpdater(scaleInitDistfElc, 0.0, 0.0, {bField3d}, {fElc})
fElc:sync()
 
tCurr = tStart
startFrame = 1
-- recalculate number density
--runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElc, oneFieldElc, bField3d}, {numDensityElc})

calcContinuousDiagnostics(tCurr, 0.0)
calcDiagnostics(tCurr, 0.0)
writeFields(startFrame-1,tCurr)

-- timing data
totalTime = 0
totalSteps = 0

diffSlvrAltElc = Updater.LenardBernsteinDiffAlternate3D2VUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  polyOrder = polyOrder,
  cfl = cfl,
  onlyIncrement = true,
  speciesMass = elcMass,
  alpha = function(t)
    return getElcAlpha(t)
  end
}
-- compute electron moments
--runUpdater(generalIntegralElcCalc, 0, 0, {fElc, oneFieldElc, bField3d}, {numDensityElc})
--runUpdater(generalIntegralElcCalc, 0, 0, {fElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
--runUpdater(generalIntegralElcCalc, 0, 0, {fElc, hamilKeElc, bField3d}, {kineticEnergyElc})
-- compute electron temperature
--runUpdater(temperatureElcCalc, 0, 0, {numDensityElc,
-- kineticEnergyElc, mom1dir3Elc}, {temperatureElc})

--runUpdater(fieldInt3d, 0, 0, {mom1dir3Elc}, {beforePosMomentumVec})
--Lucee.logInfo(string.format("Avg momentum = %g", beforePosMomentumVec:lastInsertedData()/totalVol))

--local stage1ElcEnergy = applyElcPositivity(0, 1e-8, fElc)
--Lucee.logInfo(string.format("Positivity Energy Change = %g", stage1ElcEnergy))
--Lucee.logInfo(string.format("Total Energy = %g", afterPosEnergyVec:lastInsertedData()))

--runUpdater(generalIntegralElcCalc, 0, 0, {fElc, oneFieldElc, bField3d}, {numDensityElc})
--runUpdater(fieldInt3d, 0.0, 0.0, {numDensityElc}, {totalElcCounter})
--Lucee.logInfo(string.format("Avg number = %g",totalElcCounter:lastInsertedData()/totalVol))

--fElc:write( string.format("fElc_%d.h5", 1), 0)

--f1Elc:clear(0.0)
--calcLenardBernsteinElcElc(0, 1e-8, fElc, f1Elc)
--fElcCollisions:copy(f1Elc)

--runUpdater(diffSlvrAltElc, 0, 0, {fElc, temperatureElc, temperatureElc, temperatureElc,
--bField3d, numDensityElc}, {fElcCollisions})
--runUpdater(dragSlvrElc, 0, 0, {fElc, mom1dir3Elc,
--    numDensityElc, temperatureElc, numDensityElc}, {fElcCollisions})

--fElcCollisions:write( string.format("fElcCollisions_%d.h5", 0), 0)

for frame = startFrame, nFrames do
  if endExecution == false then
    Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
    dtSuggested, tCurr = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
    writeFields(frame, tCurr)
    -- print out timing data
    Lucee.logInfo (string.format("Total time for %g steps was %g s", totalSteps, totalTime))
    Lucee.logInfo (string.format("Time per time step = %g s", totalTime/totalSteps))

    Lucee.logInfo ("")
  else
    Lucee.logInfo (string.format("-- Not advancing solution from %g", tCurr))
 end
end
