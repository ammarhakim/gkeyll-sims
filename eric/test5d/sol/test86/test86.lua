-- Input file for a 3X2V Simulation of LAPD
-- mpiexec -n 4 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i test86.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist
filename = "test86"

-- phase-space decomposition (5d)
phaseDecomp = DecompRegionCalc5D.CartProd { cuts = {18, 18, 1, 1, 1} }
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
tEnd = 2000e-6
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 2000
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
vtElc   = math.sqrt(elcTemp*eV/elcMass) -- only used for grid
vtIon   = math.sqrt(ionTemp*eV/ionMass)
omega_i = math.abs(ionCharge*B0/ionMass)

c_s = math.sqrt(elcTemp*eV/ionMass)
rho_s = c_s/omega_i

deltaX  = 100*rho_s
deltaY  = 100*rho_s
R          = 40*rho_s -- [m]
L_parallel = 36*R -- [m]
-- grid parameters: number of cells
N_X = 36
N_Y = 36
N_Z = 10
N_VPARA = 8
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
MU_UPPER_ELC = elcMass*(VPARA_UPPER_ELC*VPARA_UPPER_ELC)/(2*B0)

VPARA_UPPER_ION = 4*vtIon
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
   decomposition = phaseDecomp,
}
grid_ion_5d = Grid.RectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER_ION, MU_LOWER_ION},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER_ION, MU_UPPER_ION},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
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
basis_ion_5d = NodalFiniteElement5D.SerendipityElement {
   onGrid = grid_ion_5d,
   polyOrder = polyOrder,
   num1DGaussPoints = 2,
}
-- create 3d basis functions
basis_3d = NodalFiniteElement3D.SerendipityElement {
   onGrid = grid_3d,
   polyOrder = polyOrder,
   num1DGaussPoints = 2,
}
-- create 2d basis functions
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
-- distribution function after positivity drag
fIonPositivityDrag = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
   --writeGhost = {1,1},
}
-- delta positivity term
fIonPositivityDelta = DataStruct.Field5D {
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
hamilParaElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
hamilDerivKeElc = DataStruct.Field5D {
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
hamilParaIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
hamilDerivKeIon = DataStruct.Field5D {
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
  return B0
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
-- Compute parallel velocity derivative of hamiltonian
calcHamilKeDerivElc = Updater.SOLDerivativeCalc5D {
  onGrid = grid_elc_5d,
  basis = basis_elc_5d,
  scaleFactor = 1/elcMass,
  dir = 3,
}
runUpdater(calcHamilKeDerivElc, 0.0, 0.0, {hamilKeElc}, {hamilDerivKeElc})
-- compute parallel energy
initHamilParaElc = Updater.EvalOnNodes5D {
   onGrid = grid_elc_5d,
   basis = basis_elc_5d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,vPara,mu,t)
      return 0.5*elcMass*vPara*vPara
   end
}
runUpdater(initHamilParaElc, 0.0, 0.0, {}, {hamilParaElc})

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
-- Compute parallel velocity derivative of hamiltonian
calcHamilKeDerivIon = Updater.SOLDerivativeCalc5D {
  onGrid = grid_ion_5d,
  basis = basis_ion_5d,
  scaleFactor = 1/ionMass,
  dir = 3,
}
runUpdater(calcHamilKeDerivIon, 0.0, 0.0, {hamilKeIon}, {hamilDerivKeIon})
-- compute parallel energy
initHamilParaIon = Updater.EvalOnNodes5D {
   onGrid = grid_ion_5d,
   basis = basis_ion_5d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,vPara,mu,t)
      return 0.5*ionMass*vPara*vPara
   end
}
runUpdater(initHamilParaIon, 0.0, 0.0, {}, {hamilParaIon})

-- this function approximates a steady-state LAPD profile that is "1" at center
-- and 1/20 at the edge
function initBackgroundProfile(x,y,z)
  local r = math.sqrt(x^2 + y^2)
  local a = deltaX/2
  local edgeVal = 1/20
  if r <= a then
    return (1-edgeVal)*(1-r^2/a^2)^3 + edgeVal
  else
    return edgeVal
  end
end

function elcTempProfile(x,y,z)
  return elcTemp*initBackgroundProfile(x,y,z)
end

function fElcProfile(x,y,z,v,mu)
  return math.exp(-elcMass*v^2/(2*elcTempProfile(x,y,z)*eV))*
    math.exp(-math.abs(mu)*bFieldProfile(x)/(elcTempProfile(x,y,z)*eV))
    --*(2*math.pi*elcTempProfile(x,y,z)*eV/elcMass)^(-3/2)
end

function ionTempProfile(x,y,z)
  return ionTemp
end

function fIonProfile(x,y,z,v,mu)
  return math.exp(-ionMass*v^2/(2*ionTempProfile(x,y,z)*eV))*
    math.exp(-math.abs(mu)*bFieldProfile(x)/(ionTempProfile(x,y,z)*eV))
    --*(2*math.pi*ionTempProfile(x,y,z)*eV/ionMass)^(-3/2)
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
initfElc = Updater.EvalOnNodes5D {
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
initfIon = Updater.EvalOnNodes5D {
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

-- lua implementation of tanh function since it is removed in Lua 5.3
function tanh(x)
  if x == 0 then return 0.0 end
  local neg = false
  if x < 0 then x = -x; neg = true end
  if x < 0.54930614433405 then
    local y = x * x
    x = x + x * y *
        ((-0.96437492777225469787e0  * y +
          -0.99225929672236083313e2) * y +
          -0.16134119023996228053e4) /
        (((0.10000000000000000000e1  * y +
           0.11274474380534949335e3) * y +
           0.22337720718962312926e4) * y +
           0.48402357071988688686e4)
  else
    x = math.exp(x)
    x = 1.0 - 2.0 / (x * x + 1.0)
  end
  if neg then x = -x end
  return x
end

-- "1" for region of source, "0" outside
function sourceTempProfile(x,y,z)
  --local r = math.sqrt(x^2 + y^2)
  --return 0.5*(1 - tanh( (r-r_s)/L_s )) + 0.01
  return initBackgroundProfile(x,y,z)
end

-- "1" for region of source, "0" outside
function sourceDensityProfile(x,y,z)
  local r = math.sqrt(x^2 + y^2)
  return 0.5*(1 - tanh( (r-r_s)/L_s ))
end

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
		 return math.exp(-ionMass*v^2/(2*0.85*ionTemp*eV))*
      math.exp(-math.abs(mu)*bFieldProfile(x)/(0.85*ionTemp*eV))
	 end
}
runUpdater(initSourceIon, 0.0, 0.0, {}, {ionSource})

-- to scale source term
scaleInitSourceIon = Updater.SOLInitializeDensity {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/ionMass,
  evaluate = function(x,y,z,t)
    return sourceAmplitude*sourceDensityProfile(x,y,z)
  end
}

-- to scale source term
scaleInitDistfIon = Updater.SOLInitializeDensity {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/ionMass,
  evaluate = function(x,y,z,t)
    return 0.5*n0*initBackgroundProfile(x,y,z)
  end
}

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
   evaluate = function(x,y,z,v,mu,t)
     local elcSourceTemp = 6
		 return math.exp(-elcMass*v^2/(2*elcSourceTemp*sourceTempProfile(x,y,z)*eV))*
      math.exp(-math.abs(mu)*bFieldProfile(x)/(elcSourceTemp*sourceTempProfile(x,y,z)*eV))
	 end
}
runUpdater(initSourceElc, 0.0, 0.0, {}, {elcSource})

-- to scale source term
scaleInitSourceElc = Updater.SOLInitializeDensity {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/elcMass,
  evaluate = function(x,y,z,t)
    return sourceAmplitude*sourceDensityProfile(x,y,z)
  end
}
-- to scale source term
scaleInitDistfElc = Updater.SOLInitializeDensity {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/elcMass,
  evaluate = function(x,y,z,t)
    return 0.5*n0*initBackgroundProfile(x,y,z)
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
  evaluate = function(x,y,z,vPara,mu,t)
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
bStarYFieldElc:clear(0.0)
-- B_y^* field for ions (5D)
bStarYFieldIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}

bStarYFieldIon:clear(0.0)
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
-- to compute edge current
totalCurrent = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to compute edge current
edgeCurrent = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to compute ion density at domain boundaries in z
edgeDensityIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to compute ion flux at domain boundaries in z
edgeFluxIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to compute ion kinetic energy at domain boundaries in z
edgeKeIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to compute ion temperature at domain boundaries in z
edgeTempIon = DataStruct.Field3D {
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
-- to compute elc density at domain boundaries in z
edgeDensityElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to compute electron flux at domain boundaries in z
edgeFluxElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to compute electron kinetic energy at domain boundaries in z
edgeKeElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- to compute electron temperature at domain boundaries in z
edgeTempElc = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
edgeMomentCalcElc = Updater.MomentsAtEdges5D {
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  onGrid = grid_elc_5d,
  polyOrder = polyOrder,
  integrateGhosts = true,
}
-- to compute sheath potential on lower and upper surfaces in z (reflect electrons)
biasedSheathElcCalc = Updater.LAPDSheath5D {
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  onGrid = grid_elc_5d,
  polyOrder = polyOrder,
  speciesMass = elcMass,
  speciesCharge = -eV,
}
-- to compute sheath potential on lower and upper surfaces in z (reflect ions)
biasedSheathIonCalc = Updater.LAPDSheath5D {
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  onGrid = grid_ion_5d,
  polyOrder = polyOrder,
  speciesMass = ionMass,
  speciesCharge = eV,
}

-- stores result of potential solve
phi3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- potential used for sheath bc
phiSheath = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- stores E_parallel
eParallel3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
-- To calculate x and y derivatives of the potential
xDeriv3dCalc = Updater.SOLDerivativeCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
  dir = 0,
}
yDeriv3dCalc = Updater.SOLDerivativeCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
  dir = 1,
}
zDeriv3dCalc = Updater.SOLDerivativeCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
  dir = 2,
}
-- To store x and y derivatives of the potential
xDerivPhi3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
xDerivPhiSq3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
yDerivPhi3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
yDerivPhiSq3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
secondHamil3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
secondHamil5dElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
secondHamil5dIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
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
   bcLeft = { T ="D", V = 0.0},
   bcRight = { T ="D", V = 0.0}, -- dirichlet in x
   bcBottom = { T ="D", V = 0.0},
   bcTop = { T ="D", V = 0.0}, -- dirichlet in y
   bcBack = { T ="N", V = 0.0},
   bcFront = { T ="N", V = 0.0}, -- no bc in z
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
pbSlvrIon = Updater.PoissonBracketSimple5D {
   onGrid = grid_ion_5d,
   basis = basis_ion_5d,
   cfl = cfl,
   equation = gyroEqnIon,
   -- let solver know about additional jacobian factor
   jacobianFactor = ionMass^2*B0,
   updateDirections = {0,1,2,3},
   zeroFluxDirections = {0,1,3,4},
   fluxType = "upwind",
}
pbSlvrElc = Updater.PoissonBracketSimple5D {
   onGrid = grid_elc_5d,
   basis = basis_elc_5d,
   cfl = cfl,
   equation = gyroEqnElc,
   -- let solver know about additional jacobian factor
   jacobianFactor = elcMass^2*B0,
   updateDirections = {0,1,2,3},
   zeroFluxDirections = {0,1,3,4},
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
-- electron total pressure
pressureElc = DataStruct.Field3D {
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
-- used to calculate number density without B weight
oneFieldElc = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
oneFieldElc:clear(1.0)
oneFieldIon = DataStruct.Field5D {
   onGrid = grid_ion_5d,
   numComponents = basis_ion_5d:numNodes(),
   ghost = {1, 1},
}
oneFieldIon:clear(1.0)
-- used to store number density without B weight
kineticEnergyIon = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
kineticEnergyElc = DataStruct.Field3D {
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
-- ion parallel tempearture
temperatureParaIon = DataStruct.Field3D {
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

-- to compute pressure for electrons
pressureElcCalc = Updater.SOLPressureAtNodeCalc {
   onGrid = grid_3d,
   basis = basis_3d,
   speciesMass = elcMass,
}

generalIntegralElcCalc = Updater.SOLGeneralIntegralAtNodeCalc {--SOLWeightedProjectionCalc {
   onGrid = grid_elc_5d,
   basis5d = basis_elc_5d,
   basis3d = basis_3d,
   basis2d = basis_2d,
   scaleFactor = 2*math.pi/elcMass,
}

generalIntegralIonCalc = Updater.SOLGeneralIntegralAtNodeCalc {--SOLWeightedProjectionCalc {
   onGrid = grid_ion_5d,
   basis5d = basis_ion_5d,
   basis3d = basis_3d,
   basis2d = basis_2d,
   scaleFactor = 2*math.pi/ionMass,
}
-- to compute temperature for electrons
temperatureElcCalc = Updater.SOLTemperatureAtNodeCalc {
   onGrid = grid_3d,
   basis = basis_3d,
   speciesMass = elcMass,
}
-- to compute temperature for electrons
temperatureIonCalc = Updater.SOLTemperatureAtNodeCalc {
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
avgElcDensity = DataStruct.DynVector { numComponents = 1, }
totalIonCounter = DataStruct.DynVector { numComponents = 1, }
totalElcCounter = DataStruct.DynVector { numComponents = 1, }

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
-- diagnostic for flux across edge
totalEdgeFluxElcCalc = Updater.SOLFluxAcrossEdgeCalc {
  onGrid = grid_elc_5d,
  basis = basis_elc_5d,
  scaleFactor = 2*math.pi/elcMass^3,
  integrateGhosts = true,
}
totalEdgeFluxElc = DataStruct.DynVector { numComponents = 2, }
-- diagnostic for flux across edge
totalEdgeFluxIonCalc = Updater.SOLFluxAcrossEdgeCalc {
  onGrid = grid_ion_5d,
  basis = basis_ion_5d,
  scaleFactor = 2*math.pi/ionMass^3,
  integrateGhosts = true,
}
totalEdgeFluxIon = DataStruct.DynVector { numComponents = 2, }

-- heat flux at edge diagnostic
heatFluxAtEdgeElcCalc = Updater.SOLHeatFluxCalc {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/elcMass, -- B will be an input
}
heatFluxLowerElc = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
heatFluxUpperElc = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}

heatFluxAtEdgeIonCalc = Updater.SOLHeatFluxCalc {
  onGrid = grid_ion_5d,
  basis5d = basis_ion_5d,
  basis3d = basis_3d,
  basis2d = basis_2d,
  scaleFactor = 2*math.pi/ionMass, -- B will be an input
}
heatFluxLowerIon = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
heatFluxUpperIon = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}

-- diagnostic to look at phi3d in a single cell at all times
recordFieldInCellPhiCalc = Updater.RecordFieldInCell3D {
  onGrid = grid_3d,
  cellIndex = {18,18,2},
}
fieldInCellPhi = DataStruct.DynVector { numComponents = basis_3d:numNodes(), }
fieldInCellNumElc = DataStruct.DynVector { numComponents = basis_3d:numNodes(), }

-- total energy vs time dynvectors
elcEnergyHistory = DataStruct.DynVector { numComponents = 1, }

-- to calculate elc energy integrated over a cell
elcEnergyAtCellCalc = Updater.SOLEnergyAtCellCalc {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  scaleFactor = 2*math.pi/elcMass^3
}
elcEnergyAtCellAfterPositivity = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

ionSourceTempHistory = DataStruct.DynVector { numComponents = 1, }
elcSourceTempHistory = DataStruct.DynVector { numComponents = 1, }

-- lenard berstein code
function getElcAlpha(t)
  -- see http://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node35.html
  local n_e0 = avgElcDensity:lastInsertedData()/totalVol
  local logLambda = 6.6-0.5*math.log(n_e0/10^(20)) + 1.5*math.log(elcTemp)
  -- alpha needs to be multiplied by average n/(kTe)^(3/2) in C++ code
  return logLambda*elcCharge^4/(6*math.sqrt(2)*math.pi^(3/2)*eps0^2*math.sqrt(elcMass))
end

diffSlvrElc = Updater.LenardBernsteinDiff5DUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  cfl = cfl,
  onlyIncrement = true,
  speciesMass = elcMass,
  alpha = function(t)
    return getElcAlpha(t)
  end
}
diffSlvrElcIon = Updater.LenardBernsteinDiff5DUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  cfl = cfl,
  onlyIncrement = true,
  speciesMass = elcMass,
  alpha = function(t)
    -- collision frequency will be 1/1.96 that for elc-elc collisions
    return getElcAlpha(t)/1.96
  end
}
dragSlvrElc = Updater.LenardBernsteinDrag5DUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  cfl = cfl,
  onlyIncrement = true,
  speciesMass = elcMass,
  alpha = function(t)
    return getElcAlpha(t)
  end
}

lenardBernsteinScaleElcUpdater = Updater.SOLPositivityScaleCellUpdater {
  onGrid = grid_elc_5d,
  basis5d = basis_elc_5d,
  basis3d = basis_3d,
  positivityChecks = false,
}

fElcCollisions = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
gVParElcCollisions = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
gMuElcCollisions = DataStruct.Field5D {
   onGrid = grid_elc_5d,
   numComponents = basis_elc_5d:numNodes(),
   ghost = {1, 1},
}
elcEnergyAtCellDrag = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
elcEnergyAtCellDiff = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
weightedFieldInt3d = Updater.SOLTotalIntegralCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
}
ohmicHeatingVec = DataStruct.DynVector { numComponents = 1, }
-- stores product of E_parallel times J by projection
ohmicHeatingField = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

initPhiOffset = Updater.EvalOnNodes3D {
  onGrid = grid_3d,
  basis = basis_3d,
  shareCommonNodes = false,
  evaluate = function(x,y,z,t)
    local r = math.sqrt(x^2 + y^2)
    if r >= r_s then
      return V_bias
    else
      return 0
    end
  end
}
-- offset potential for sheath bc
phiOffset = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
runUpdater(initPhiOffset, 0, 0, {}, {phiOffset})
phiOffset:sync()
phiOffset:write( string.format("phiOffset_%d.h5", 0), 0)

oneField3d = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
oneField3d:clear(1.0)

phiCurrentProductAtEdgeCalc = Updater.SOLFluxAcrossEdgeCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
}

currentIntegratedAtEdgeCalc = Updater.SOLFluxAcrossEdgeCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
}

phiCurrentProductAtEdgeVec = DataStruct.DynVector { numComponents = 2, }
currentIntegratedAtEdgeVec = DataStruct.DynVector { numComponents = 2, }

function applyPositivity(tCurr, myDt, distfElc, distfIon)  -- apply positivity to ions
  local ionStatus = runUpdater(positivityIonUpdater, tCurr, myDt, {bField3d}, {distfIon})
  -- apply positivity to elcs
  local elcStatus = runUpdater(positivityElcUpdater, tCurr, myDt, {bField3d}, {distfElc})
  
  return true
end

function calcPotential(nIon, nElc, phiOut)
  -- calculate average ion density
  runUpdater(fieldInt3d, 0.0, 0.0, {nIon}, {avgIonDensity})
  local n_i0 = 2.325e17--avgIonDensity:lastInsertedData()/totalVol--2.325e17
  -- solve for potential
  numDensityDelta:combine(1.0, nElc, -1.0, nIon)
  numDensityDelta:scale(elcTemp/(n_i0))

  local statusPhi = runUpdater(fieldSolver, 0.0, 0.0, {numDensityDelta}, {phiOut})
  if statusPhi == false then
    Lucee.logInfo(string.format("-- phi3d solve failed --"))
    Lucee.logInfo(string.format("-- n_i0 = %g --", n_i0))
    return false
  end

  return true
end

function calcHamiltonianElc(phiIn, hamilOut)
  -- copy phi3d into a 5d field
  runUpdater(copy3dTo5dElc, 0.0, 0.0, {phiIn}, {hamilOut})
  hamilOut:scale(elcCharge)
  -- compute second order hamiltonian correction
  --calcSecondOrderHamiltonian(phiIn, secondHamil3d)
  --runUpdater(copy3dTo5dElc, 0.0, 0.0, {secondHamil3d}, {secondHamil5dElc})
  --hamilOut:accumulate(-0.5*elcMass/B0^2,secondHamil5dElc)
  -- add the kinetic energy component
  hamilOut:accumulate(1.0, hamilKeElc)
  hamilOut:sync()
end

productProjection3dCalc = Updater.ASquaredProjection3D {
  onGrid = grid_3d,
  basis = basis_3d,
}

function calcSecondOrderHamiltonian(phiIn, hamil3dOut)
  runUpdater(xDeriv3dCalc, 0.0, 0.0, {phiIn}, {xDerivPhi3d})
  runUpdater(yDeriv3dCalc, 0.0, 0.0, {phiIn}, {yDerivPhi3d})
  -- square each field by projection
  runUpdater(productProjection3dCalc, 0.0, 0.0, {xDerivPhi3d,xDerivPhi3d}, {xDerivPhiSq3d})
  runUpdater(productProjection3dCalc, 0.0, 0.0, {yDerivPhi3d,yDerivPhi3d}, {yDerivPhiSq3d})
  hamil3dOut:combine(1,xDerivPhiSq3d, 1,yDerivPhiSq3d)
end

function calcHamiltonianIon(phiIn, hamilOut)
  -- copy phi3d into a 5d field
  runUpdater(copy3dTo5dIon, 0.0, 0.0, {phiIn}, {hamilOut})
  hamilOut:scale(ionCharge)
  -- compute second order hamiltonian correction
  --calcSecondOrderHamiltonian(phiIn, secondHamil3d)
  --runUpdater(copy3dTo5dIon, 0.0, 0.0, {secondHamil3d}, {secondHamil5dIon})
  --hamilOut:accumulate(-0.5*ionMass/B0^2,secondHamil5dIon)
  -- add the kinetic energy component
  hamilOut:accumulate(1.0, hamilKeIon)
  hamilOut:sync()
end

function applyBcToDistF(tCurr, myDt, fElc, fIon, phiIn)
  fIon:sync()
  fElc:sync()
  -- apply zero inflow boundary conditions in z to ions (not sure if needed)
  runUpdater(bcLowerZIon, tCurr, myDt, {}, {fIon})
  runUpdater(bcUpperZIon, tCurr, myDt, {}, {fIon})
  -- apply zero inflow boundary conditions in z to electrons (not sure if needed)
  runUpdater(bcLowerZElc, tCurr, myDt, {}, {fElc})
  runUpdater(bcUpperZElc, tCurr, myDt, {}, {fElc})
  -- compute potential for sheath bc
  phiSheath:combine(1.0, phiIn, -1.0, phiOffset)
  -- apply biased sheath boundary conditions to electron distribution function
  sheathStatus = runUpdater(biasedSheathElcCalc, tCurr, myDt, {phiSheath, hamilDerivKeElc}, {fElc})
  -- apply biased sheath boundary conditions to ion distribution function
  sheathStatus = runUpdater(biasedSheathIonCalc, tCurr, myDt, {phiSheath, hamilDerivKeIon}, {fIon})
  -- sync cells in distribution function (not sure if needed)
  fElc:sync()
  fIon:sync()
end

function calcLenardBernsteinElcElc(tCurr, myDt, fInitial, fFinal)
  -- calculate inputs for collision operator
  -- compute electron moments
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fInitial, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityElc}, {avgElcDensity})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fInitial, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fInitial, hamilKeElc, bField3d}, {kineticEnergyElc})
  -- compute electron temperature
  negativeTemperature = runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc,
    kineticEnergyElc, mom1dir3Elc}, {temperatureElc})

  if (negativeTemperature == false) then
    endExecution = true
    Lucee.logInfo(string.format("-- Ending execution because negative temperature. --"))
    return false, 0.0
  end

  local dragStatus, dragDtSuggested = runUpdater(dragSlvrElc, tCurr, myDt, {fInitial, mom1dir3Elc,
    numDensityElc, temperatureElc, numDensityElc}, {fElcCollisions})
  -- compute energy integrated in a cell for drag term only
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcCollisions, jacobianFieldElc, hamilKeElc},
    {elcEnergyAtCellDrag, elcEnergyHistory})
  -- accumulate drag term times dt to distribution function
  fFinal:accumulate(myDt, fElcCollisions)
  
  local diffStatus, diffDtSuggested = runUpdater(diffSlvrElc, tCurr, myDt, {fInitial, temperatureElc,
    bField3d, numDensityElc}, {fElcCollisions, gVParElcCollisions, gMuElcCollisions})

  if (diffStatus == false or dragStatus == false) then
    return false, math.min(dragDtSuggested, diffDtSuggested)
  end

  -- compute energy integrated in a cell for diffusion term only
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcCollisions, jacobianFieldElc, hamilKeElc},
    {elcEnergyAtCellDiff, elcEnergyHistory})

  elcEnergyAtCellDiff:scale(-1.0/myDt)
  -- need a 3d field with all 0's as input to scaler. doesn't matter which field
  elcEnergyAtCellAfterPositivity:clear(0.0)
  -- accumulate correct fraction of diffusion term to distribution function
  runUpdater(lenardBernsteinScaleElcUpdater, tCurr, myDt, {fElcCollisions, elcEnergyAtCellDrag,
    elcEnergyAtCellAfterPositivity, elcEnergyAtCellDiff}, {fFinal})

  if (diffStatus == false or dragStatus == false) then
    return false, math.min(dragDtSuggested, diffDtSuggested)
  end

  return true, math.min(dragDtSuggested, diffDtSuggested)
end

function calcLenardBernsteinElcIon(tCurr, myDt, fInitialElc, fInitialIon, fFinal)
  -- calculate inputs for collision operator
  -- compute electron moments
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fInitialElc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityElc}, {avgElcDensity})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fInitialElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fInitialElc, hamilKeElc, bField3d}, {kineticEnergyElc})
  -- compute electron temperature
  negativeTemperature = runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc,
    kineticEnergyElc, mom1dir3Elc}, {temperatureElc})
  -- compute number density and first parallel velocity moment for ions
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fInitialIon, oneFieldIon, bField3d}, {numDensityIon})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fInitialIon, hamilDerivKeIon, bField3d}, {mom1dir3Ion})

  if (negativeTemperature == false) then
    endExecution = true
    Lucee.logInfo(string.format("-- Ending execution because negative temperature. --"))
    return false, 0.0
  end

  -- compute drag term
  local dragStatus, dragDtSuggested = runUpdater(dragSlvrElc, tCurr, myDt, {fInitialElc, mom1dir3Ion,
    numDensityIon, temperatureElc, numDensityElc}, {fElcCollisions})
  -- compute energy integrated in a cell for drag term only
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcCollisions, jacobianFieldElc, hamilKeElc},
    {elcEnergyAtCellDrag, elcEnergyHistory})
  -- accumulate drag term times dt to distribution function
  fFinal:accumulate(myDt, fElcCollisions)
  
  local diffStatus, diffDtSuggested = runUpdater(diffSlvrElcIon, tCurr, myDt, {fInitialElc, temperatureElc,
    bField3d, numDensityElc}, {fElcCollisions, gVParElcCollisions, gMuElcCollisions})

  if (diffStatus == false or dragStatus == false) then
    return false, math.min(dragDtSuggested, diffDtSuggested)
  end

  -- compute energy integrated in a cell for diffusion term only
  runUpdater(elcEnergyAtCellCalc, tCurr, myDt, {fElcCollisions, jacobianFieldElc, hamilKeElc},
    {elcEnergyAtCellDiff, elcEnergyHistory})

  elcEnergyAtCellDiff:scale(-1.0/myDt)
  -- need a 3d field with all 0's as input to scaler. doesn't matter which field
  elcEnergyAtCellAfterPositivity:clear(0.0)
  -- accumulate correct fraction of diffusion term to distribution function
  runUpdater(lenardBernsteinScaleElcUpdater, tCurr, myDt, {fElcCollisions, elcEnergyAtCellDrag,
    elcEnergyAtCellAfterPositivity, elcEnergyAtCellDiff}, {fFinal})

  if (diffStatus == false or dragStatus == false) then
    return false, math.min(dragDtSuggested, diffDtSuggested)
  end

  return true, math.min(dragDtSuggested, diffDtSuggested)
end

endExecution = false

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
  -- RK stage 1
  local myStatusElc, myDtSuggestedElc = runUpdater(pbSlvrElc, tCurr, myDt, {fElc, hamilElc}, {f1Elc})
  local myStatusLBElc, myDtSuggestedLBElc = calcLenardBernsteinElcElc(tCurr, myDt, fElc, f1Elc)
  local myStatusLBElcIon, myDtSuggestedLBElcIon = calcLenardBernsteinElcIon(tCurr, myDt, fElc, fIon, f1Elc)
  local myStatusIon, myDtSuggestedIon = runUpdater(pbSlvrIon, tCurr, myDt, {fIon, hamilIon}, {f1Ion})

  local netStatus = (myStatusElc == false or myStatusIon == false) or (myStatusLBElc == false or myStatusLBElcIon == false)
  local netDtSuggested = math.min(math.min(math.min(myDtSuggestedElc, myDtSuggestedIon), myDtSuggestedLBElc), myDtSuggestedLBElcIon)
  if (netStatus == true) then
    return false, netDtSuggested
  end
  -- add on particle source
  f1Elc:accumulate(myDt, elcSource)
  f1Ion:accumulate(myDt, ionSource)
  local positivityStatus = applyPositivity(tCurr, myDt, f1Elc, f1Ion)
  if (positivityStatus == false) then
    fElc:copy(f1Elc)
    fIon:copy(f1Ion)
    return false, netDtSuggested
  end
  -- calculate density of each species for potential solve
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {f1Elc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {f1Ion, oneFieldIon, bField3d}, {numDensityIon})
  -- solve for potential
  local potentialStatus = calcPotential(numDensityIon, numDensityElc, phi3d)
  if (potentialStatus == false) then
    endExecution = true
    calcContinuousDiagnostics(tCurr, myDt)
    return false, netDtSuggested
  end
  -- calculate new hamiltonians
  calcHamiltonianElc(phi3d, hamilElc)
  calcHamiltonianIon(phi3d, hamilIon)
  -- apply boundary conditions to distribution functions
  applyBcToDistF(tCurr, myDt, f1Elc, f1Ion, phi3d)

  -- RK stage 2
  myStatusElc, myDtSuggestedElc = runUpdater(pbSlvrElc, tCurr, myDt, {f1Elc, hamilElc}, {fNewElc})
  myStatusLBElc, myDtSuggestedLBElc = calcLenardBernsteinElcElc(tCurr, myDt, f1Elc, fNewElc)
  myStatusLBElcIon, myDtSuggestedLBElcIon = calcLenardBernsteinElcIon(tCurr, myDt, f1Elc, f1Ion, fNewElc)
  myStatusIon, myDtSuggestedIon = runUpdater(pbSlvrIon, tCurr, myDt, {f1Ion, hamilIon}, {fNewIon})

  netStatus = (myStatusElc == false or myStatusIon == false) or (myStatusLBElc == false or myStatusLBElcIon == false)
  netDtSuggested = math.min(math.min(math.min(myDtSuggestedElc, myDtSuggestedIon), myDtSuggestedLBElc), myDtSuggestedLBElcIon)

  if (netStatus == true) then
    return false, netDtSuggested
  end
  -- add on particle source
  fNewElc:accumulate(myDt, elcSource)
  fNewIon:accumulate(myDt, ionSource)
  positivityStatus = applyPositivity(tCurr, myDt, fNewElc, fNewIon)
  if (positivityStatus == false) then
    fElc:copy(fNewElc)
    fIon:copy(fNewIon)
    return false, netDtSuggested
  end
  
  f1Elc:combine(3.0/4.0, fElc, 1.0/4.0, fNewElc)
  f1Ion:combine(3.0/4.0, fIon, 1.0/4.0, fNewIon)
  -- calculate density of each species for potential solve
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {f1Elc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {f1Ion, oneFieldIon, bField3d}, {numDensityIon})
  -- solve for potential
  potentialStatus = calcPotential(numDensityIon, numDensityElc, phi3d)
  if (potentialStatus == false) then
    endExecution = true
    calcContinuousDiagnostics(tCurr, myDt)
    return false, netDtSuggested
  end
  -- calculate new hamiltonians
  calcHamiltonianElc(phi3d, hamilElc)
  calcHamiltonianIon(phi3d, hamilIon)
  -- apply boundary conditions to distribution functions
  applyBcToDistF(tCurr, myDt, f1Elc, f1Ion, phi3d)

  -- RK stage 3
  myStatusElc, myDtSuggestedElc = runUpdater(pbSlvrElc, tCurr, myDt, {f1Elc, hamilElc}, {fNewElc})
  myStatusLBElc, myDtSuggestedLBElc = calcLenardBernsteinElcElc(tCurr, myDt, f1Elc, fNewElc)
  myStatusLBElcIon, myDtSuggestedLBElcIon = calcLenardBernsteinElcIon(tCurr, myDt, f1Elc, f1Ion, fNewElc)
  myStatusIon, myDtSuggestedIon = runUpdater(pbSlvrIon, tCurr, myDt, {f1Ion, hamilIon}, {fNewIon})

  netStatus = (myStatusElc == false or myStatusIon == false) or (myStatusLBElc == false or myStatusLBElcIon == false)
  netDtSuggested = math.min(math.min(math.min(myDtSuggestedElc, myDtSuggestedIon), myDtSuggestedLBElc), myDtSuggestedLBElcIon)

  if (netStatus == true) then
    return false, netDtSuggested
  end
  -- add on particle source
  fNewElc:accumulate(myDt, elcSource)
  fNewIon:accumulate(myDt, ionSource)
  positivityStatus = applyPositivity(tCurr, myDt, fNewElc, fNewIon)
  if (positivityStatus == false) then
    fElc:copy(fNewElc)
    fIon:copy(fNewIon)
    return false, netDtSuggested
  end
  f1Elc:combine(1.0/3.0, fElc, 2.0/3.0, fNewElc)
  f1Ion:combine(1.0/3.0, fIon, 2.0/3.0, fNewIon)
  -- calculate density of each species for potential solve
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {f1Elc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {f1Ion, oneFieldIon, bField3d}, {numDensityIon})
  -- solve for potential
  potentialStatus = calcPotential(numDensityIon, numDensityElc, phi3d)
  if (potentialStatus == false) then
    endExecution = true
    calcContinuousDiagnostics(tCurr, myDt)
    return false, netDtSuggested
  end
  -- calculate new hamiltonians
  calcHamiltonianElc(phi3d, hamilElc)
  calcHamiltonianIon(phi3d, hamilIon)
  -- apply boundary conditions to distribution functions
  applyBcToDistF(tCurr, myDt, f1Elc, f1Ion, phi3d)

  fElc:copy(f1Elc)
  fIon:copy(f1Ion)

  -- calculate continuous time diagnostics
  calcContinuousDiagnostics(tCurr, myDt)

  return true, netDtSuggested
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

-- function containing diagnostics that must be called at the end of every frame
function calcDiagnostics(tCurr, myDt)
  -- compute electron moments
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilKeElc, bField3d}, {kineticEnergyElc})
  -- compute electron temperature
  negativeTemperature = runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc,
    kineticEnergyElc, mom1dir3Elc}, {temperatureElc})
  -- convert to ev
  temperatureElc:scale(1/eV)

  -- compute parallel temperature
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilParaElc, bField3d}, {kineticEnergyElc})
  kineticEnergyElc:scale(2) -- get rid of 1/2 factor in hamilParaElc
  runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc, kineticEnergyElc, mom1dir3Elc},
    {temperatureParaElc})
  temperatureParaElc:scale(1/eV)

  -- compute electron pressure
  runUpdater(pressureElcCalc, tCurr, myDt, {numDensityElc, kineticEnergyElc, mom1dir3Elc}, {pressureElc})

  -- compute ion moments
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, oneFieldIon, bField3d}, {numDensityIon})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, hamilDerivKeIon, bField3d}, {mom1dir3Ion})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, hamilKeIon, bField3d}, {kineticEnergyIon})
  -- compute ion temperature
  negativeTemperature = runUpdater(temperatureIonCalc, tCurr, myDt, {numDensityIon,
    kineticEnergyIon, mom1dir3Ion}, {temperatureIon})
  -- convert to ev
  temperatureIon:scale(1/eV)

  -- compute parallel temperature for ions
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, hamilParaIon, bField3d}, {kineticEnergyIon})
  kineticEnergyIon:scale(2) -- get rid of 1/2 factor in hamilParaIon
  runUpdater(temperatureIonCalc, tCurr, myDt, {numDensityIon, kineticEnergyIon, mom1dir3Ion},
    {temperatureParaIon})
  temperatureParaIon:scale(1/eV)

  -- calculate quantites on edge
  runUpdater(edgeMomentCalcElc, tCurr, myDt, {fElc, oneFieldElc}, {edgeDensityElc})
  edgeDensityElc:scale(2*math.pi*B0/elcMass)
  runUpdater(edgeMomentCalcElc, tCurr, myDt, {fElc, hamilDerivKeElc}, {edgeFluxElc})
  edgeFluxElc:scale(2*math.pi*B0/elcMass)
  runUpdater(edgeMomentCalcElc, tCurr, myDt, {fElc, hamilKeElc}, {edgeKeElc})
  edgeKeElc:scale(2*math.pi*B0/elcMass)
  runUpdater(temperatureElcCalc, tCurr, myDt, {edgeDensityElc,
    edgeKeElc, edgeFluxElc}, {edgeTempElc})
  edgeTempElc:scale(1/eV)

  -- calculate quantites on edge
  runUpdater(edgeMomentCalcIon, tCurr, myDt, {fIon, oneFieldIon}, {edgeDensityIon})
  edgeDensityIon:scale(2*math.pi*B0/ionMass)
  runUpdater(edgeMomentCalcIon, tCurr, myDt, {fIon, hamilDerivKeIon}, {edgeFluxIon})
  edgeFluxIon:scale(2*math.pi*B0/ionMass)
  runUpdater(edgeMomentCalcIon, tCurr, myDt, {fIon, hamilKeIon}, {edgeKeIon})
  edgeKeIon:scale(2*math.pi*B0/ionMass)
  runUpdater(temperatureIonCalc, tCurr, myDt, {edgeDensityIon,
    edgeKeIon, edgeFluxIon}, {edgeTempIon})
  edgeTempIon:scale(1/eV)

  -- compute net edge current
  edgeCurrent:combine(elcCharge, edgeFluxElc, ionCharge, edgeFluxIon)

  -- compute heat flux at edge
  runUpdater(heatFluxAtEdgeElcCalc, tCurr, myDt, {fElc,hamilDerivKeElc,hamilKeElc,bField3d},
    {heatFluxLowerElc, heatFluxUpperElc})
  runUpdater(heatFluxAtEdgeIonCalc, tCurr, myDt, {fIon,hamilDerivKeIon,hamilKeIon,bField3d},
    {heatFluxLowerIon, heatFluxUpperIon})

  -- compute ohmic heating field
  runUpdater(productProjection3dCalc, tCurr, myDt, {totalCurrent,eParallel3d}, {ohmicHeatingField})
end

-- function to compute source temperature
function calcSourceDiagnostics(tCurr, myDt, elcSource, ionSource)
  -- compute electron moments
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {elcSource, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {elcSource, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralElcCalc, tCurr, myDt, {elcSource, hamilKeElc, bField3d}, {kineticEnergyElc})
  -- compute electron temperature
  runUpdater(temperatureElcCalc, tCurr, myDt, {numDensityElc,
    kineticEnergyElc, mom1dir3Elc}, {temperatureElc})
  -- convert to ev
  temperatureElc:scale(1/eV)

  -- compute electron moments
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {ionSource, oneFieldIon, bField3d}, {numDensityIon})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {ionSource, hamilDerivKeIon, bField3d}, {mom1dir3Ion})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {ionSource, hamilKeIon, bField3d}, {kineticEnergyIon})
  -- compute electron temperature
  runUpdater(temperatureIonCalc, tCurr, myDt, {numDensityIon,
    kineticEnergyIon, mom1dir3Ion}, {temperatureIon})
  -- convert to ev
  temperatureIon:scale(1/eV)

  numDensityElc:write( string.format("elcSource_%d.h5", 0), 0)
  numDensityIon:write( string.format("ionSource_%d.h5", 0), 0)

  -- compute average ion temperature
  local volume = (X_UPPER-X_LOWER)*(Y_UPPER-Y_LOWER)*(Z_UPPER-Z_LOWER)
  runUpdater(fieldInt3d, tCurr, myDt, {temperatureElc}, {elcSourceTempHistory})
  temperatureElc:write( string.format("elcSourceTemp_%d.h5", 0), 0)

  Lucee.logInfo(string.format( 'Average source elc temp = %g eV',elcSourceTempHistory:lastInsertedData()/volume ))
  runUpdater(fieldInt3d, tCurr, myDt, {temperatureIon}, {ionSourceTempHistory})
  temperatureIon:write( string.format("ionSourceTemp_%d.h5", 0), 0)
  Lucee.logInfo(string.format( 'Average source ion temp = %g eV',ionSourceTempHistory:lastInsertedData()/volume ))

  -- compute average source kinetic energy
  runUpdater(fieldInt3d, tCurr, myDt, {kineticEnergyElc}, {elcSourceTempHistory})
  Lucee.logInfo(string.format( 'Average source elc energy = %g W',elcSourceTempHistory:lastInsertedData() ))

  runUpdater(fieldInt3d, tCurr, myDt, {kineticEnergyIon}, {ionSourceTempHistory})
  Lucee.logInfo(string.format( 'Average source ion energy = %g W',ionSourceTempHistory:lastInsertedData() ))
end

-- function containing diagnostics that must be called at the end of every time step
function calcContinuousDiagnostics(tCurr, myDt)
  -- store phi in a cell
  runUpdater(recordFieldInCellPhiCalc, tCurr, myDt, {phi3d}, {fieldInCellPhi})
  runUpdater(recordFieldInCellPhiCalc, tCurr, myDt, {numDensityElc}, {fieldInCellNumElc})
  -- compute total guiding center density by integrating over whole domain
  -- note: includes a volume factor that should be removed in post-processing
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityElc}, {totalElcCounter})
  runUpdater(fieldInt3d, tCurr, myDt, {numDensityIon}, {totalIonCounter})

  -- Compute outward flux
  runUpdater(totalEdgeFluxIonCalc, tCurr, myDt, {fIon, jacobianFieldIon, hamilDerivKeIon},
    {totalEdgeFluxIon})
  runUpdater(totalEdgeFluxElcCalc, tCurr, myDt, {fElc, jacobianFieldElc, hamilDerivKeElc},
    {totalEdgeFluxElc})

  runUpdater(generalIntegralElcCalc, tCurr, myDt, {fElc, hamilDerivKeElc, bField3d}, {mom1dir3Elc})
  runUpdater(generalIntegralIonCalc, tCurr, myDt, {fIon, hamilDerivKeIon, bField3d}, {mom1dir3Ion})
  -- compute ohmic heating over entire box
  runUpdater(zDeriv3dCalc, tCurr, myDt, {phi3d},{eParallel3d})
  eParallel3d:scale(-1)
  totalCurrent:combine(elcCharge, mom1dir3Elc, ionCharge, mom1dir3Ion)
  runUpdater(weightedFieldInt3d, tCurr, myDt, {totalCurrent,eParallel3d}, {ohmicHeatingVec})

  -- compute net edge current
  runUpdater(edgeMomentCalcElc, tCurr, myDt, {fElc, hamilDerivKeElc}, {edgeFluxElc})
  edgeFluxElc:scale(2*math.pi*B0/elcMass)
  runUpdater(edgeMomentCalcIon, tCurr, myDt, {fIon, hamilDerivKeIon}, {edgeFluxIon})
  edgeFluxIon:scale(2*math.pi*B0/ionMass)
  edgeCurrent:combine(elcCharge, edgeFluxElc, ionCharge, edgeFluxIon)
  -- compute phi vs current integrated at edges
  runUpdater(phiCurrentProductAtEdgeCalc, tCurr, myDt, {edgeCurrent, phi3d}, {phiCurrentProductAtEdgeVec})
  -- compute current integrated at edges
  runUpdater(currentIntegratedAtEdgeCalc, tCurr, myDt, {edgeCurrent, oneField3d}, {currentIntegratedAtEdgeVec})

  return true
end

-- create output file directories
function createDirs()
  if Lucee.getRank() == 0 then
    os.execute("mkdir density")
    os.execute("mkdir temperature")
    os.execute("mkdir pressure")
    os.execute("mkdir phi")
    os.execute("mkdir edge")
    os.execute("mkdir dynvector")
    os.execute("mkdir heatflux")
    os.execute("mkdir ohmicheating")
  end
end

-- write data to H5 file
function writeFields(frameNum, tCurr)
  numDensityElc:write( string.format("nElc_%d.h5", frameNum), tCurr)
  temperatureElc:write( string.format("tElc_%d.h5", frameNum), tCurr)
  --temperatureParaElc:write( string.format("tParaElc_%d.h5", frameNum), tCurr)
  
  numDensityIon:write( string.format("nIon_%d.h5", frameNum), tCurr)
  temperatureIon:write( string.format("tIon_%d.h5", frameNum), tCurr)
  --temperatureParaIon:write( string.format("tParaIon_%d.h5", frameNum), tCurr)
  
  pressureElc:write( string.format("pressureElc_%d.h5", frameNum), tCurr)
  
  phi3d:write( string.format("phi_%d.h5", frameNum), tCurr)

  fElc:write( string.format("fElc_%d.h5", frameNum), tCurr)
  fIon:write( string.format("fIon_%d.h5", frameNum), tCurr)

  ohmicHeatingVec:write( string.format("ohmicHeatingVec_%d.h5", frameNum), tCurr)
  ohmicHeatingField:write( string.format("ohmicHeating_%d.h5", frameNum), tCurr)
  
  totalIonCounter:write( string.format("totalIonCounter_%d.h5", frameNum), tCurr)
  totalElcCounter:write( string.format("totalElcCounter_%d.h5", frameNum), tCurr)

  -- edge fields
  edgeCurrent:write( string.format("edgeCurrent_%d.h5", frameNum), tCurr)

  edgeDensityIon:write( string.format("edgeDensityIon_%d.h5", frameNum), tCurr)
  edgeFluxIon:write( string.format("edgeFluxIon_%d.h5", frameNum), tCurr)
  edgeTempIon:write( string.format("edgeTempIon_%d.h5", frameNum), tCurr)

  edgeDensityElc:write( string.format("edgeDensityElc_%d.h5", frameNum), tCurr)
  edgeFluxElc:write( string.format("edgeFluxElc_%d.h5", frameNum), tCurr)
  edgeTempElc:write( string.format("edgeTempElc_%d.h5", frameNum), tCurr)

  fieldInCellPhi:write( string.format("phiInCell_%d.h5", frameNum), tCurr)
  fieldInCellNumElc:write( string.format("nElcInCell_%d.h5", frameNum), tCurr)

  heatFluxUpperElc:write( string.format("heatFluxUpperElc_%d.h5", frameNum), tCurr)
  heatFluxUpperIon:write( string.format("heatFluxUpperIon_%d.h5", frameNum), tCurr)

  phiCurrentProductAtEdgeVec:write( string.format("phiCurrentProductAtEdgeVec_%d.h5", frameNum), tCurr)
  currentIntegratedAtEdgeVec:write( string.format("currentIntegratedAtEdgeVec_%d.h5", frameNum), tCurr)

  if Lucee.getRank() == 0 then
    os.execute("mv ".. filename .. "_" .. string.format("nElc_%d.h5", frameNum) .. " density")
    os.execute("mv ".. filename .. "_" .. string.format("tElc_%d.h5", frameNum) .. " temperature")
    os.execute("mv ".. filename .. "_" .. string.format("nIon_%d.h5", frameNum) .. " density")
    os.execute("mv ".. filename .. "_" .. string.format("tIon_%d.h5", frameNum) .. " temperature")
    os.execute("mv ".. filename .. "_" .. string.format("pressureElc_%d.h5", frameNum) .. " pressure")
    os.execute("mv ".. filename .. "_" .. string.format("phi_%d.h5", frameNum) .. " phi")
    os.execute("mv ".. filename .. "_" .. string.format("ohmicHeatingVec_%d.h5", frameNum) .. " ohmicheating")
    os.execute("mv ".. filename .. "_" .. string.format("ohmicHeating_%d.h5", frameNum) .. " ohmicheating")
    os.execute("mv ".. filename .. "_" .. string.format("totalIonCounter_%d.h5", frameNum) .. " dynvector")
    os.execute("mv ".. filename .. "_" .. string.format("totalElcCounter_%d.h5", frameNum) .. " dynvector")
    os.execute("mv ".. filename .. "_" .. string.format("edgeCurrent_%d.h5", frameNum) .. " edge")
    os.execute("mv ".. filename .. "_" .. string.format("edgeDensityIon_%d.h5", frameNum) .. " edge")
    os.execute("mv ".. filename .. "_" .. string.format("edgeFluxIon_%d.h5", frameNum) .. " edge")
    os.execute("mv ".. filename .. "_" .. string.format("edgeTempIon_%d.h5", frameNum) .. " edge")
    os.execute("mv ".. filename .. "_" .. string.format("edgeDensityElc_%d.h5", frameNum) .. " edge")
    os.execute("mv ".. filename .. "_" .. string.format("edgeFluxElc_%d.h5", frameNum) .. " edge")
    os.execute("mv ".. filename .. "_" .. string.format("edgeTempElc_%d.h5", frameNum) .. " edge")
    os.execute("mv ".. filename .. "_" .. string.format("phiInCell_%d.h5", frameNum) .. " dynvector")
    os.execute("mv ".. filename .. "_" .. string.format("nElcInCell_%d.h5", frameNum) .. " dynvector")
    os.execute("mv ".. filename .. "_" .. string.format("heatFluxUpperElc_%d.h5", frameNum) .. " heatflux")
    os.execute("mv ".. filename .. "_" .. string.format("heatFluxUpperIon_%d.h5", frameNum) .. " heatflux")
    os.execute("mv ".. filename .. "_" .. string.format("phiCurrentProductAtEdgeVec_%d.h5", frameNum) .. " dynvector")
    os.execute("mv ".. filename .. "_" .. string.format("currentIntegratedAtEdgeVec_%d.h5", frameNum) .. " dynvector")
  end
end

-- apply positivity preservation to the sources
applyPositivity(0.0, 0.0, elcSource, ionSource)
-- scale source terms to have correct density
runUpdater(scaleInitSourceIon, 0.0, 0.0, {bField3d}, {ionSource})
runUpdater(scaleInitSourceElc, 0.0, 0.0, {bField3d}, {elcSource})
-- ensure that source terms are in exact balance
runUpdater(generalIntegralElcCalc, 0.0, 0.0, {elcSource, oneFieldElc, bField3d}, {numDensityElc})
-- store electron density temporarily in avgIonDensity structure
runUpdater(fieldInt3d, 0.0, 0.0, {numDensityElc}, {avgIonDensity})
totalSourceElc = avgIonDensity:lastInsertedData()
runUpdater(generalIntegralIonCalc, 0.0, 0.0, {ionSource, oneFieldIon, bField3d}, {numDensityIon})
-- store ion density temporarily in avgIonDensity structure
runUpdater(fieldInt3d, 0.0, 0.0, {numDensityIon}, {avgIonDensity})
totalSourceIon = avgIonDensity:lastInsertedData()
-- scale ions to match electrons
ionSource:scale(totalSourceElc/totalSourceIon)
--print( string.format("ionSource scaled by %g (%g/%g)",totalSourceElc/totalSourceIon,totalSourceElc,totalSourceIon) )
ionSource:sync()
elcSource:sync()

-- write out some source fields
calcSourceDiagnostics(0.0, 0.0, elcSource, ionSource)

-- create output file directories
createDirs()

if Lucee.IsRestarting then
  fElc:read("fElc_" .. Lucee.RestartFrame .. ".h5")
  fElc:sync()
  fIon:read("fIon_" .. Lucee.RestartFrame .. ".h5")
  fIon:sync()

  startFrame = Lucee.RestartFrame + 1
  tCurr = tStart + tFrame*Lucee.RestartFrame
  runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralIonCalc, 0.0, 0.0, {fIon, oneFieldIon, bField3d}, {numDensityIon})
  -- solve for potential
  calcPotential(numDensityIon, numDensityElc, phi3d)
  -- calculate new hamiltonians
  calcHamiltonianElc(phi3d, hamilElc)
  calcHamiltonianIon(phi3d, hamilIon)
  -- apply boundary conditions to distribution functions
  applyBcToDistF(0.0, 0.0, fElc, fIon, phi3d)
else
  -- Make sure distribution functions are positive
  applyPositivity(0.0, 0.0, fElc, fIon)
  runUpdater(scaleInitDistfElc, 0.0, 0.0, {bField3d}, {fElc})
  runUpdater(scaleInitDistfIon, 0.0, 0.0, {bField3d}, {fIon})
  -- multiply some small random density perturbations
  runUpdater(multiply5dCalcIon, 0.0, 0.0, {fIon, fIonInitialPerturb}, {fDupIon})
  fIon:copy(fDupIon)
  fIon:sync()
  runUpdater(multiply5dCalcElc, 0.0, 0.0, {fElc, fElcInitialPerturb}, {fDupElc})
  fElc:copy(fDupElc)
  fElc:sync()
  -- make sure total ions and electrons are equal
  runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElc, oneFieldElc, bField3d}, {numDensityElc})
  -- store electron density temporarily in avgIonDensity structure
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityElc}, {avgIonDensity})
  totalInitElc = avgIonDensity:lastInsertedData()
  runUpdater(generalIntegralIonCalc, 0.0, 0.0, {fIon, oneFieldIon, bField3d}, {numDensityIon})
  -- store ion density temporarily in avgIonDensity structure
  runUpdater(fieldInt3d, 0.0, 0.0, {numDensityIon}, {avgIonDensity})
  totalInitIon = avgIonDensity:lastInsertedData()
  -- scale ions to match electrons
  fIon:scale(totalInitElc/totalInitIon)
  fIon:sync()
  fElc:sync()
  
  tCurr = tStart
  startFrame = 1
  -- recalculate number density
  runUpdater(generalIntegralElcCalc, 0.0, 0.0, {fElc, oneFieldElc, bField3d}, {numDensityElc})
  runUpdater(generalIntegralIonCalc, 0.0, 0.0, {fIon, oneFieldIon, bField3d}, {numDensityIon})
  -- solve for potential
  calcPotential(numDensityIon, numDensityElc, phi3d)
  -- calculate new hamiltonians
  calcHamiltonianElc(phi3d, hamilElc)
  calcHamiltonianIon(phi3d, hamilIon)
  -- apply boundary conditions to distribution functions
  applyBcToDistF(0.0, 0.0, fElc, fIon, phi3d)
end

calcDiagnostics(tCurr, 0.0)
calcContinuousDiagnostics(tCurr, 0.0)
writeFields(startFrame-1,tCurr)

for frame = startFrame, nFrames do
  if endExecution == false then
    Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
    dtSuggested, tCurr = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
    writeFields(frame, tCurr)
    Lucee.logInfo ("")
  else
    --Lucee.logInfo (string.format("-- Not advancing solution from %g", tCurr))
  end
end
