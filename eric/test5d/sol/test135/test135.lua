-- Input file for a 3X2V Simulation of NSTX
-- mpiexec -n 1 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i test135.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist
filename = "test135"

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
tEnd = 1000e-6
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 4000
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

-- physical parameters
eV        = Lucee.ElementaryCharge
elcCharge = -eV
ionCharge = eV
eps0      = Lucee.Epsilon0 -- permittivity of free space
ionMass   = 2.014*Lucee.ProtonMass -- (Deuterium)
elcMass   = ionMass/400 -- REDUCED ELECTRON MASS Lucee.ElectronMass
elcTemp   = 40 -- [eV]
ionTemp   = 40 -- [eV]
B0 = 0.5  -- [T]
R0 = 0.85 -- [m], Major Radius
a0 = 0.5
n0 = 7*10^18 -- [1/m^3]
V_bias = 0 -- [V]
P_SOL = 5.4e6/5 -- [W]
xSource = -0.04 + R0 + a0 -- [m], source start coordinate
lambdaSource = 0.02 -- [m], characteristic length scale of density and temperature
-- derived parameters
vtElc   = math.sqrt(elcTemp*eV/elcMass) -- only used for grid
vtIon   = math.sqrt(ionTemp*eV/ionMass)
omega_i = math.abs(ionCharge*B0*(R0/(R0+a0))/ionMass)
c_s = math.sqrt(elcTemp*eV/ionMass)
rho_s = c_s/omega_i
epsilonDensity = 1e-7 -- density floor to enforce

deltaX  = 50*rho_s
deltaY  = 200*rho_s
qR = 8 -- [m, from Myra 2016]
L_parallel = 2*2*math.pi/3*qR
-- grid parameters: number of cells
N_X = 18
N_Y = 72
N_Z = 10
N_VPARA = 10
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = (R0+a0) - deltaX/2
X_UPPER = (R0+a0) + deltaX/2
Y_LOWER = -deltaY/2
Y_UPPER = deltaY/2
Z_LOWER = -L_parallel/2
Z_UPPER = L_parallel/2
-- source parameters
L_s = 0.5*rho_s
r_s = 20*rho_s
sourceAmplitude = 1--n0*c_s/L_parallel
-- total volume
totalVol = (X_UPPER-X_LOWER)*(Y_UPPER-Y_LOWER)*(Z_UPPER-Z_LOWER)

VPARA_UPPER_ELC = 4*vtElc
VPARA_LOWER_ELC = -VPARA_UPPER_ELC
MU_LOWER_ELC = 0
MU_UPPER_ELC = 0.75*elcMass*(VPARA_UPPER_ELC*VPARA_UPPER_ELC)/(2*B0*(R0/(R0+a0)))

VPARA_UPPER_ION = 4*vtIon
VPARA_LOWER_ION = -VPARA_UPPER_ION
MU_LOWER_ION = 0
MU_UPPER_ION = 0.75*ionMass*(VPARA_UPPER_ION*VPARA_UPPER_ION)/(2*B0*(R0/(R0+a0)))

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

-- 3d spatial grid
grid_3d = Grid.RectCart3D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER},
   cells = {N_X, N_Y, N_Z},
   periodicDirs = {1}, -- periodic in y only
   decomposition = confDecomp,
}

-- create 3d basis functions
basis_3d = NodalFiniteElement3D.SerendipityElement {
   onGrid = grid_3d,
   polyOrder = polyOrder,
   num1DGaussPoints = 3,
}
numDensityElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

S_n = 3.92032e+22
-- electron density profile
function initialElectronDensityProfile(x,y,z)
  -- Peak temperature of the source goes into sound speed
  local nPeak = 4*math.sqrt(3)/3*math.sqrt(ionMass/(60*eV))*S_n/2*(L_parallel/4)
  local L_s = L_parallel/2
  if math.abs(z) <= L_s/2 then
    return nPeak*(1 + math.sqrt(1-(z/L_s/2)^2))/2
  else
    return nPeak/2
  end
end

initNumDensityElc = Updater.EvalOnNodes3D {
   onGrid = grid_3d,
   basis = basis_3d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		 return initialElectronDensityProfile(x,y,z)
	 end
}
runUpdater(initNumDensityElc, 0.0, 0.0, {}, {numDensityElc})
numDensityElc:write( string.format("nElc_%d.h5", 0), 0)


function initialIonDensityProfile(x)
  local dxOrig = 50*rho_s/18 -- Original grid spacing in x on which Greg's notes are based
  local elcDens = n0--initialElectronDensityProfile(x)
  --n0 -- Assume electron density is flat for now
  local w = 2*dxOrig
  local phi_2 = 200 -- 4*Te
  local phi_5 = 40 -- 4*Te
  local w3 = dxOrig -- Width of region 3
  local phi_4 = (phi_2*3*lambdaSource + phi_5*w3)/(3*lambdaSource + w3)

  if (x <= R0 + a0 - 7*dxOrig) then
    -- cubic falloff region
    return elcDens - 6*n0*rho_s^2/elcTemp*(x-( R0+a0-7*dxOrig ))/w^3*phi_2
  elseif (x <= R0 + a0 - 3*dxOrig) then
    -- constant region
    return elcDens
  elseif (x <= R0 + a0 - 2*dxOrig) then
    return elcDens + n0/elcTemp*6*rho_s^2/w3^3*(x-(R0+a0-3*dxOrig))*(phi_2-phi_4)
  elseif (x <= R0 + a0 + 7*dxOrig) then
    -- exponential region
    return elcDens - n0/elcTemp*(rho_s/lambdaSource)^2*math.exp(-(x-(R0+a0-2*dxOrig))/lambdaSource)*(phi_4-phi_5)
  else
    -- cubic falloff region
    return elcDens + 6*n0*rho_s^2/elcTemp*(x-( R0+a0+7*dxOrig ))/w^3*phi_5
  end
end

numDensityIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
numDensityDelta = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

initNumDensityIon = Updater.EvalOnNodes3D {
   onGrid = grid_3d,
   basis = basis_3d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		 return initialIonDensityProfile(x)
	 end
}
runUpdater(initNumDensityIon, 0.0, 0.0, {}, {numDensityIon})
numDensityElc:clear(n0)

numDensityElc:sync()
numDensityIon:sync()

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
   sourceNodesShared = false, -- default true
   solutionNodesShared = false, -- default true
   writeStiffnessMatrix = false,
   modifierConstant = 0,--kPerpTimesRhoS^2,
   laplacianWeight = rho_s^2,
   isGyroKineticPoisson = true,
   -- boundary conditions to apply
   periodicDirs = {1},
   bcLeft = { T ="D", V = 0.0},
   bcRight = { T ="D", V = 0.0}, -- dirichlet in x
   bcBack = { T ="N", V = 0.0},
   bcFront = { T ="N", V = 0.0}, -- no bc in z
}
function calcPotential(nIon, nElc, phiOut)
  -- calculate average ion density
  --local n_i0 = 3.5e17--avgIonDensity:lastInsertedData()/totalVol--2.325e17
  local n_i0 = n0--avgIonDensity:lastInsertedData()/totalVol--2.325e17
  -- solve for potential
  numDensityDelta:combine(1.0, nElc, -1.0, nIon)
  numDensityDelta:scale(elcTemp/(n_i0))

  local statusPhi = runUpdater(fieldSolver, 0.0, 0.0, {numDensityDelta}, {phiOut})

  if statusPhi == false then
    Lucee.logInfo(string.format("-- phi3d solve failed --"))
    return false
  end

  return true
end

calcPotential(numDensityIon, numDensityElc, phi3d)

numDensityDelta:combine(1.0, numDensityIon, -1.0, numDensityElc)
numDensityDelta:write( string.format("nDelta_%d.h5", 0), 0)

phi3d:write( string.format("phi_%d.h5", 0), 0)
