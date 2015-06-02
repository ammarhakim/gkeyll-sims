-- Input file for SOL problem with kinetic ions and electrons
-- Not sure what the math.pi's are all about, but this file was used for paper figure
-- 5-31-2015: Used to test calculation of cutoff reflection boundary conditions.

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.1

-- parameters to control time-stepping
tStart = 0.0
tEnd = 350e-6
nFrames = 5

-- physical constants
-- eletron mass (kg)
elcMass = Lucee.ElectronMass
-- electron volt to joules
eV = Lucee.ElementaryCharge
-- Deuterium ion mass (kg)
ionMass = 2.014*Lucee.ProtonMass
-- signed ion charge
ionCharge = Lucee.ElementaryCharge
-- signed electron charge
elcCharge = -Lucee.ElementaryCharge

-- physical input parameters
kPerpTimesRho = 0.2
-- pedestal density (1/m^3)
nPed = 5e19
-- pedestal (source) temperature (eV)
tPed = 1500
-- Fixed value of Te that must be independent of time (eV)
Te0 = 75
-- ELM pulse duration (seconds)
tELM = 200e-6
-- Parallel length (m)
lParallel = 40
-- Source length (m)
lSource = 25
-- Particle source proportionality factor
A = 1.2
-- Magnetic field (Tesla)
B0 = 2

mu0 = Lucee.Mu0

-- Derived parameters
-- Ion cyclotron frequency
Omega_i = Lucee.ElementaryCharge*B0/ionMass
-- Ion sound speed
c_i = math.sqrt(Te0*eV/ionMass)
-- Sound-based gyroradius
rho_i = c_i/Omega_i
-- k_perpendicular
kPerp = kPerpTimesRho/rho_i
-- thermal velocities
vtElc = math.sqrt(tPed*eV/elcMass)
vtIon = math.sqrt(tPed*eV/ionMass)
-- Pedestal sound speed (m/s)
cPed = math.sqrt(2*tPed*eV/ionMass)
-- Particle source
Sn   = A*nPed*cPed/lSource
-- number of cells
N_Z, N_VPARA, N_MU = 8, 32, 16
-- domain extents
Z_LOWER, Z_UPPER = -lParallel, lParallel
-- electron velocity extents
VPARA_UPPER_ELC = math.min(4, 2.5*math.sqrt(N_VPARA/4))*vtElc
VPARA_LOWER_ELC = -VPARA_UPPER_ELC
MU_LOWER_ELC = 0
MU_UPPER_ELC = math.min(8, 4*math.sqrt(N_MU/2))*elcMass*vtElc*vtElc/B0
-- ion velocity extents
VPARA_UPPER_ION = math.min(4, 2.5*math.sqrt(N_VPARA/4))*vtIon
VPARA_LOWER_ION = -VPARA_UPPER_ION
MU_LOWER_ION = 0
MU_UPPER_ION = math.min(8, 4*math.sqrt(N_MU/2))*ionMass*vtIon*vtIon/B0

-- phase space grid for electrons
gridElc = Grid.RectCart3D {
   lower = {Z_LOWER, VPARA_LOWER_ELC, MU_LOWER_ELC},
   upper = {Z_UPPER, VPARA_UPPER_ELC, MU_UPPER_ELC},
   cells = {N_Z, N_VPARA, N_MU},
}
-- phase space grid for ions
gridIon = Grid.RectCart3D {
   lower = {Z_LOWER, VPARA_LOWER_ION, MU_LOWER_ION},
   upper = {Z_UPPER, VPARA_UPPER_ION, MU_UPPER_ION},
   cells = {N_Z, N_VPARA, N_MU},
}
-- spatial grid for both species
grid_1d = Grid.RectCart1D {
   lower = {Z_LOWER},
   upper = {Z_UPPER},
   cells = {N_Z},
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

-- create FEM nodal basis
basisElc = NodalFiniteElement3D.SerendipityElement {
   onGrid = gridElc,
   polyOrder = polyOrder,
}
basisIon = NodalFiniteElement3D.SerendipityElement {
   onGrid = gridIon,
   polyOrder = polyOrder,
}

-- spatial FEM nodal basis
basis_1d = NodalFiniteElement1D.LagrangeTensor {
   onGrid = grid_1d,
   polyOrder = polyOrder,
}

-- distribution function for electrons
distfElc = DataStruct.Field3D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
-- distribution function for ions
distfIon = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}
-- unit density ion distribution function
distfIonUnit = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}

-- Electron and ion number densitites
numDensityElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
numDensityIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
numDensityIon3dDg = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}
-- Moments for electrons and ions
mom0Elc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom0Ion = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom1ParaElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom1ParaIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom2ParaElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom2ParaIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom1MuElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom1MuIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
tPerpElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
tPerpIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- Updater to compute electron number density
mom0CalcElc = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = gridElc,
   basis3d = basisElc,
   basis1d = basis_1d,
   moment = 0,
}
mom0CalcIon = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = gridIon,
   basis3d = basisIon,
   basis1d = basis_1d,
   moment = 0,
}

-- Updater to compute parallel velocity moments
mom1ParaCalcElc = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = gridElc,
   basis3d = basisElc,
   basis1d = basis_1d,
   moment = 1,
   momentDirection = 1,
}
mom1ParaCalcIon = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = gridIon,
   basis3d = basisIon,
   basis1d = basis_1d,
   moment = 1,
   momentDirection = 1,
}
mom2ParaCalcElc = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = gridElc,
   basis3d = basisElc,
   basis1d = basis_1d,
   moment = 2,
   momentDirection = 1,
}
mom2ParaCalcIon = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = gridIon,
   basis3d = basisIon,
   basis1d = basis_1d,
   moment = 2,
   momentDirection = 1,
}

-- Updater to compute <mu> moment
mom1MuCalcElc = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = gridElc,
   basis3d = basisElc,
   basis1d = basis_1d,
   moment = 1,
   momentDirection = 2,
}
mom1MuCalcIon = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = gridIon,
   basis3d = basisIon,
   basis1d = basis_1d,
   moment = 1,
   momentDirection = 2,
}

function maxwellian(mass, vt, vPara, mu)
  local vtPerp = math.sqrt(tPed*eV/mass)
  return 1/(2*Lucee.Pi*vtPerp^2*math.sqrt(2*Lucee.Pi*vt^2))*math.exp(-vPara^2/(2*vt^2))*
    math.exp(-mu*B0/(tPed*eV))
end

-- Maxwellian (Right half only)
function maxwellianRight(mass, vt, vPara, mu)
  local vtPerp = math.sqrt(tPed*eV/mass)
  if vPara >=0 then
    return 1/(2*Lucee.Pi*vtPerp^2*math.sqrt(2*Lucee.Pi*vt^2))*math.exp(-vPara^2/(2*vt^2))*
      math.exp(-mu*B0/(tPed*eV))
  else
    return 0
  end
end

-- Maxwellian (Left half only)
function maxwellianLeft(mass, vt, vPara, mu)
  local vtPerp = math.sqrt(tPed*eV/mass)
  if vPara <= 0 then
    return 1/(2*Lucee.Pi*vtPerp^2*math.sqrt(2*Lucee.Pi*vt^2))*math.exp(-vPara^2/(2*vt^2))*
      math.exp(-mu*B0/(tPed*eV))
  else
    return 0
  end
end

-- Return initial Te(x) in eV
function initialElectronTemp(x)
  return Te0
end

-- Return initial ne(x) in 1/m^3
function initialElcDens(x)
  local backgroundDens = 0.7 + 0.3*(1 - math.abs(x/lParallel))
  if math.abs(x) < lSource/2 then
       backgroundDens = backgroundDens + 0.5*math.cos(math.pi*x/lSource)
     end
  backgroundDens = backgroundDens*1e19
  return backgroundDens
end

-- Return initial Ti(x) in eV
function initialIonTemp(x)
  local backgroundTemp = 100 + 45*(1 - math.abs(x/lParallel))
  if math.abs(x) < lSource/2 then
       backgroundTemp = backgroundTemp + 30*math.cos(math.pi*x/lSource)
     end
  return backgroundTemp
end

-- updater to initialize distribution function
initDistfElc = Updater.ProjectOnNodalBasis3D {
   onGrid = gridElc,
   -- basis functions to use
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(z,vPara,mu,t)
     local vTe = math.sqrt(initialElectronTemp(z)*eV/elcMass)
     return initialElcDens(z)*maxwellian(elcMass, vTe, vPara, mu)
	 end
}
runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})
-- Calculate initial electron density field
runUpdater(mom0CalcElc, 0.0, 0.0, {distfElc}, {numDensityElc})
numDensityElc:scale(2*math.pi*B0/elcMass)

-- updater to initialize distribution function
initDistfIon = Updater.ProjectOnNodalBasis3D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(z,vPara,mu,t)
     local backgroundTemp = initialIonTemp(z)
     local nHat = 2
     local vTi = math.sqrt(2*math.pi/(math.pi-1)*backgroundTemp*eV/ionMass)

     if z > lSource/2 then
       return nHat*maxwellianRight(ionMass, vTi, vPara, mu)
     elseif z < -lSource/2 then
       return nHat*maxwellianLeft(ionMass, vTi, vPara, mu)
     else
       -- Must be between the source boundaries, use a linear combo.
       return ((lSource/2 + z)*nHat*maxwellianRight(ionMass, vTi, vPara, mu) + 
          (lSource/2 - z)*nHat*maxwellianLeft(ionMass, vTi, vPara, mu))/lSource
     end
	 end
}
runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIonUnit})
-- Calculate initial ion density field
initIonDensityCalc = Updater.SOLIonDensityInitialization {
  onGrid = grid_1d,
  basis = basis_1d,
  kPerpTimesRho = kPerpTimesRho,
}
runUpdater(initIonDensityCalc, 0.0, 0.0, {numDensityElc}, {numDensityIon})

-- Copy numDensityIon to a 3d field
copyTo3DIonDg = Updater.NodalCopy1DTo3DFieldUpdater {
  onGrid = gridIon,
  basis1d = basis_1d,
  basis3d = basisIon,
  polyOrder = polyOrder,
}
runUpdater(copyTo3DIonDg, 0.0, 0.0, {numDensityIon}, {numDensityIon3dDg})

-- Multiply with gaussian to get total distfIon
multiply3dFields = Updater.FieldArithmeticUpdater3D {
  onGrid = gridIon,
  basis = basisIon,
  evaluate = function(An,Bn,t)
    return An*Bn
  end,
}
runUpdater(multiply3dFields, 0.0, 0.0, {numDensityIon3dDg, distfIonUnit}, {distfIon})

-- extra fields for performing RK update
distfNewElc = DataStruct.Field3D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distf1Elc = DataStruct.Field3D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distfNewIon = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}
distf1Ion = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}

-- Electron Hamiltonian
hamilElc = DataStruct.Field3D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
-- Ion Hamiltonian
hamilIon = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}

-- Fields to store the kinetic energy part of the ion and electron hamiltonians
hamilElcKeDg = DataStruct.Field3D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
hamilIonKeDg = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}

-- updater to initialize electron kinetic energy term in Hamiltonian
initHamilElcKe = Updater.EvalOnNodes3D {
   onGrid = gridElc,
   basis = basisElc,
   shareCommonNodes = false,
   evaluate = function (z,vPara,mu,t)
      return elcMass*vPara^2/2 + mu*B0
   end
}
runUpdater(initHamilElcKe, 0.0, 0.0, {}, {hamilElcKeDg})

-- updater to initialize ion kinetic energy term in Hamiltonian
initHamilIonKe = Updater.EvalOnNodes3D {
   onGrid = gridIon,
   basis = basisIon,
   shareCommonNodes = false,
   evaluate = function (z,vPara,mu,t)
      return ionMass*vPara^2/2 + mu*B0
   end
}
runUpdater(initHamilIonKe, 0.0, 0.0, {}, {hamilIonKeDg})

-- particle source (ELECTRONS)
particleSourceElc = DataStruct.Field3D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
particleSourceElc:clear(0.0)
-- updater to fill out particle source
particleSourceUpdaterElc = Updater.ProjectOnNodalBasis3D {
   onGrid = gridElc,
   -- basis functions to use
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(z,vPara,mu,t)
    if math.abs(z) < lSource/2 then
      if t < tELM then
        return Sn*math.cos(math.pi*z/lSource)*maxwellian(elcMass, math.sqrt(tPed*eV/elcMass), vPara, mu)
      else
        return Sn/9*math.cos(math.pi*z/lSource)*maxwellian(elcMass, math.sqrt(210*eV/elcMass), vPara, mu)
      end
    else
      return 0
    end
	 end
}
runUpdater(particleSourceUpdaterElc, 0.0, 0.0, {}, {particleSourceElc})

-- particle source (IONS)
particleSourceIon = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
particleSourceIon:clear(0.0)
-- updater to fill out particle source
particleSourceUpdaterIon = Updater.ProjectOnNodalBasis3D {
   onGrid = gridIon,
   -- basis functions to use
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(z,vPara,mu,t)
    if math.abs(z) < lSource/2 then
      if t < tELM then
        return Sn*math.cos(math.pi*z/lSource)*maxwellian(ionMass, math.sqrt(tPed*eV/ionMass), vPara, mu)
      else
        return Sn/9*math.cos(math.pi*z/lSource)*maxwellian(ionMass, math.sqrt(260*eV/ionMass), vPara, mu)
      end
    else
      return 0
    end
	 end
}
runUpdater(particleSourceUpdaterIon, 0.0, 0.0, {}, {particleSourceIon})

-- keeps track of whether to update particlesources again
postELM = false

-- define electron equation to solve
pbElcEqn = PoissonBracketEquation.SOL3D {
  speciesMass = elcMass,
}
pbSlvrElc = Updater.PoissonBracketSimple3D {
  onGrid = gridElc,
  basis = basisElc,
  cfl = cfl,
  equation = pbElcEqn,
  updateDirections = {0,1},
  zeroFluxDirections = {1,2},
  onlyIncrement = false,
  fluxType = "upwind",
}
-- define ion equation to solve
pbIonEqn = PoissonBracketEquation.SOL3D {
  speciesMass = ionMass,
}
pbSlvrIon = Updater.PoissonBracketSimple3D {
  onGrid = gridIon,
  basis = basisIon,
  cfl = cfl,
  equation = pbIonEqn,
  updateDirections = {0,1},
  zeroFluxDirections = {1,2},
  onlyIncrement = false,
  fluxType = "upwind",
}

-- Output dynvector for reflectingBc (always size 2)
cutoffVelocities = DataStruct.DynVector { numComponents = 2, }
-- Temporary fields so we don't get intermediate rk3 values
cutoffVelocities1 = DataStruct.DynVector { numComponents = 2, }
cutoffVelocities2 = DataStruct.DynVector { numComponents = 2, }

-- compute moments from distribution function
function calcMoments(curr, dt, distfElcIn, distfIonIn)
  -- moment 0
  runUpdater(mom0CalcElc, curr, dt, {distfElcIn}, {numDensityElc})
  numDensityElc:scale(2*math.pi*B0/elcMass)
  runUpdater(mom0CalcIon, curr, dt, {distfIonIn}, {numDensityIon})
  numDensityIon:scale(2*math.pi*B0/ionMass)
  -- moment 1 para
  runUpdater(mom1ParaCalcElc, curr, dt, {distfElcIn}, {mom1ParaElc})
  mom1ParaElc:scale(2*math.pi*B0/elcMass)
  runUpdater(mom1ParaCalcIon, curr, dt, {distfIonIn}, {mom1ParaIon})
  mom1ParaIon:scale(2*math.pi*B0/ionMass)
  -- moment 2 para
  runUpdater(mom2ParaCalcElc, curr, dt, {distfElcIn}, {mom2ParaElc})
  mom2ParaElc:scale(2*math.pi*B0/elcMass)
  runUpdater(mom2ParaCalcIon, curr, dt, {distfIonIn}, {mom2ParaIon})
  mom2ParaIon:scale(2*math.pi*B0/ionMass)
  -- moment 1 mu, scaled to be equal to <v_\perp^2>
  runUpdater(mom1MuCalcElc, curr, dt, {distfElcIn}, {mom1MuElc})
  mom1MuElc:scale(2*math.pi*B0/elcMass*2*B0/elcMass)
  runUpdater(mom1MuCalcIon, curr, dt, {distfIonIn}, {mom1MuIon})
  mom1MuIon:scale(2*math.pi*B0/ionMass*2*B0/ionMass)
end

-- field to store continuous potential in 1D
phi1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {1, 1},
}
-- continuous potential after it has been modified
phi1dAfterBc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {1, 1},
}

phi1dDg = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- updater to compute phi electrostatically
--electrostaticPhiCalc = Updater.ElectrostaticContPhiUpdater {
--   onGrid = grid_1d,
--   basis = basis_1d,
--   kPerpTimesRho = kPerpTimesRho,
--   Te0 = Te0,
--}

-- updater to ensure phi at boundary has value of the sheath potential
setPhiAtBoundaryCalc = Updater.SetPhiAtBoundaryUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
  elcMass = elcMass,
  elcCharge = eV,
}

-- Dynvectors to store 0-3rd moments at left and right edges
momentsAtEdgesElc = DataStruct.DynVector { numComponents = 8, }
momentsAtEdgesIon = DataStruct.DynVector { numComponents = 8, }

momentsAtEdgesElcCalc = Updater.MomentsAtEdges3DUpdater {
  onGrid = gridElc,
  basis = basisElc,
  scaleFactor = 2*math.pi*B0/elcMass;
}

momentsAtEdgesIonCalc = Updater.MomentsAtEdges3DUpdater {
  onGrid = gridIon,
  basis = basisIon,
  scaleFactor = 2*math.pi*B0/ionMass
}

-- dynvector for heat flux at edge
heatFluxAtEdge = DataStruct.DynVector { numComponents = 6, }
-- dynvector for sheath power transmission coefficients
sheathCoefficients = DataStruct.DynVector { numComponents = 7, }

-- to compute total particle energy
heatFluxAtEdgeCalc = Updater.KineticHeatFluxAtEdge3DUpdater {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   ionMass = ionMass,
   electronMass = elcMass,
   -- Enable calculation of sheath coefficients
   computeSheathCoefficient = false,
}

-- updater to move a 1d discontinuous field to a 1d continuous field
--contFromDisContCalc = Updater.ContFromDisCont1D {
--   onGrid = grid_1d,
--   basis = basis_1d,
--}

-- compute hamiltonian for electrons
function calcHamiltonianElc(curr, dt, phiDgIn, hamilOut)
   -- clear out fields (is this needed?)
   hamilOut:clear(0.0)

   -- Accumulate q*Phi contribution to hamiltonian
   runUpdater(copyTo3DElcDg, curr, dt, {phiIn}, {hamilOut})
   hamilOut:scale(elcCharge)

   -- Accumulate the KE hamiltonian to the full hamiltonian
   hamilOut:accumulate(1.0, hamilElcKeDg)
end

-- to compute second order hamiltonian term
--mhdHamiltonianCalc = Updater.MHDHamiltonianUpdater {
--  onGrid = grid_1d,
--  basis = basis_1d,
--}
-- store result of mhdHamiltonianCalc
--mhdHamiltonian1dDg = DataStruct.Field1D {
--   onGrid = grid_1d,
--   numComponents = basis_1d:numNodes(),
--   ghost = {1, 1},
--}
-- result of mhdHamiltonianCalc as continuous field
--mhdHamiltonian1d = DataStruct.Field1D {
--   onGrid = grid_1d,
--  numComponents = basis_1d:numExclusiveNodes(),
--   ghost = {1, 1},
--}
--mhdHamiltonian2d = DataStruct.Field3D {
--   onGrid = gridIon,
--   numComponents = basisIon:numExclusiveNodes(),
--   ghost = {1, 1},
--}
--mhdHamiltonian2dDg = DataStruct.Field2D {
--   onGrid = gridIon,
--   numComponents = basisIon:numNodes(),
--   ghost = {1, 1},
--}
-- stores mhd hamiltonian contribution to total energy
mhdHamiltonianEnergy = DataStruct.DynVector { numComponents = 1, }

-- compute hamiltonian for electrons
function calcHamiltonianIon(curr, dt, phiDgIn, hamilOut)
   -- clear out fields (is this needed?)
   hamilOut:clear(0.0)

   -- Accumulate q*Phi contribution to hamiltonian
   runUpdater(copyTo3DIonDg, curr, dt, {phiDgIn}, {hamilOut})
   hamilOut:scale(ionCharge)

   -- Accumulate the KE hamiltonian to the full hamiltonian
   hamilOut:accumulate(1.0, hamilIonKeDg)

   -- compute second order hamiltonian term
   --runUpdater(copyCToD1d, curr, dt, {phiIn}, {phi1dDg})
   --runUpdater(mhdHamiltonianCalc, curr, dt, {phi1dDg, numDensityIon}, {mhdHamiltonian1dDg})
   --runUpdater(contFromDisContCalc, curr, dt, {mhdHamiltonian1dDg}, {mhdHamiltonian1d})
   --runUpdater(copyTo2DIon, curr, dt, {mhdHamiltonian1d}, {mhdHamiltonian2d})
   --mhdHamiltonian2d:scale(-0.5*kPerpTimesRho*kPerpTimesRho*Lucee.ElementaryCharge/Te0)
   --hamilOut:accumulate(1.0, mhdHamiltonian2d)
end

-- Outflow BCs
function getRepTbl(pOrder, val)
   if pOrder == 1 then
      return {val, val, val, val, val, val, val, val}
   elseif pOrder == 2 then
      return {val, val, val, val, val, val, val, val,
              val, val, val, val,
              val, val, val, val, val, val, val, val}
   end
end

function getCountTbl(pOrder, val)
   if pOrder == 1 then
      return {0, 1, 2, 3, 4, 5, 6, 7}
   elseif pOrder == 2 then
      return {0, 1, 2, 3, 4, 5, 6, 7,
              8, 9, 10, 11,
              12, 13, 14, 15, 16, 17, 18, 19}
   end
end

bcConst = BoundaryCondition.Const { 
   components = getCountTbl(polyOrder),
   values = getRepTbl(polyOrder, 0.0),
}

function makeBcObjIon()
   local bcLower = Updater.Bc3D {
      onGrid = gridIon,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "lower",
   }
   local bcUpper = Updater.Bc3D {
      onGrid = gridIon,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "upper",
   }
   return bcLower, bcUpper
end

-- make objects to apply BCs
bcLowerIon, bcUpperIon = makeBcObjIon()

-- updater to apply boundary condition on distribution function
reflectingBc = Updater.SOL3DElectrostaticDistFuncReflectionBc {
   onGrid = gridElc,
   basis = basisElc,
   basis1d = basis_1d,
   edge = "both",
   -- physical parameters
   scaleFactor = 2*math.pi*B0/elcMass,
}


-- apply boundary conditions
function applyBc(curr, dt, fldElc, fldIon, cutoffV)
  runUpdater(bcLowerIon, curr, dt, {}, {fldIon})
  runUpdater(bcUpperIon, curr, dt, {}, {fldIon})
  -- Use reflecting BC for the electrons
  runUpdater(reflectingBc, curr, dt, {mom0Ion, mom1ParaIon}, {fldElc, cutoffV})
  -- (zero flux bc's applied via poisson updater)
end

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   -- compute moments at edges
   runUpdater(momentsAtEdgesElcCalc, curr, dt, {distfElc, hamilKeElcDg}, {momentsAtEdgesElc})
   runUpdater(momentsAtEdgesIonCalc, curr, dt, {distfIon, hamilKeIonDg}, {momentsAtEdgesIon})
   -- modify phi so that is is equal to phi_s at edge
   -- TODO: change it so it takes phi1dDg as input instead of phi1d
   runUpdater(setPhiAtBoundaryCalc, curr, dt, {phi1d, cutoffVelocities}, {phi1dAfterBc})
   runUpdater(copyCToD1d, curr, dt, {phi1dAfterBc}, {phi1dDg})
   -- compute n*T_perp for both species
   tPerpElc:copy(mom1MuElc)
   tPerpElc:scale(0.5*elcMass/eV)
   tPerpIon:copy(mom1MuIon)
   tPerpIon:scale(0.5*ionMass/eV)
   -- compute heat flux at edges
   runUpdater(heatFluxAtEdgeCalc, curr, dt, {phi1dDg, momentsAtEdgesElc,
    momentsAtEdgesIon, numDensityElc, numDensityIon, tPerpElc, tPerpIon},
    {heatFluxAtEdge, sheathCoefficients})
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local statusElc, dtSuggestedElc
   local statusIon, dtSuggestedIon

   if (postELM == false) then
     if tCurr + myDt > tELM then
       postELM = true
       runUpdater(particleSourceUpdaterElc, tCurr, myDt, {}, {particleSourceElc})
       runUpdater(particleSourceUpdaterIon, tCurr, myDt, {}, {particleSourceIon})
     end
   end
   
   -- RK stage 1
   statusElc, dtSuggestedElc = runUpdater(pbSlvrElc, tCurr, myDt, {distfElc, hamilElc}, {distf1Elc})
   statusIon, dtSuggestedIon = runUpdater(pbSlvrIon, tCurr, myDt, {distfIon, hamilIon}, {distf1Ion})
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   
   distf1Elc:accumulate(myDt, particleSourceElc)
   distf1Ion:accumulate(myDt, particleSourceIon)
   
   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion, cutoffVelocities1)
   
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon}, {phi1d})
   -- Copy the continuous phi back onto a dg field
   runUpdater(copyCToD1d, tCurr, myDt, {phi1d}, {phi1dDg})

   -- Compute Hamiltonian
   calcHamiltonianElc(tCurr, myDt, phi1dDg, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1dDg, hamilIon)

   -- RK stage 2
   statusElc, dtSuggestedElc = runUpdater(pbSlvrElc, tCurr, myDt, {distf1Elc, hamilElc}, {distfNewElc})
   statusIon, dtSuggestedIon = runUpdater(pbSlvrIon, tCurr, myDt, {distf1Ion, hamilIon}, {distfNewIon})
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   
   distfNewElc:accumulate(myDt, particleSourceElc)
   distfNewIon:accumulate(myDt, particleSourceIon)
   
   distf1Elc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   distf1Ion:combine(3.0/4.0, distfIon, 1.0/4.0, distfNewIon)
   
   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion, cutoffVelocities2)
   
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon}, {phi1d})
   -- Copy the continuous phi back onto a dg field
   runUpdater(copyCToD1d, tCurr, myDt, {phi1d}, {phi1dDg})

   -- Compute Hamiltonian
   calcHamiltonianElc(tCurr, myDt, phi1dDg, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1dDg, hamilIon)

   -- RK stage 3
   statusElc, dtSuggestedElc = runUpdater(pbSlvrElc, tCurr, myDt, {distf1Elc, hamilElc}, {distfNewElc})
   statusIon, dtSuggestedIon = runUpdater(pbSlvrIon, tCurr, myDt, {distf1Ion, hamilIon}, {distfNewIon})
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end

   distfNewElc:accumulate(myDt, particleSourceElc)
   distfNewIon:accumulate(myDt, particleSourceIon)

   distf1Elc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   distf1Ion:combine(1.0/3.0, distfIon, 2.0/3.0, distfNewIon)
   
   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion, cutoffVelocities)

   distfElc:copy(distf1Elc)
   distfIon:copy(distf1Ion)

   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon}, {phi1d})
   -- Copy the continuous phi back onto a dg field
   runUpdater(copyCToD1d, tCurr, myDt, {phi1d}, {phi1dDg})

   -- Compute Hamiltonian
   calcHamiltonianElc(tCurr, myDt, phi1dDg, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1dDg, hamilIon)

   return true, math.min(dtSuggestedElc, dtSuggestedIon)
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      distfDupElc:copy(distfElc)
      distfDupIon:copy(distfIon)
      hamilDupElc:copy(hamilElc)
      hamilDupIon:copy(hamilIon)

      if (tCurr+myDt > tEnd) then
	      myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	      distfElc:copy(distfDupElc)
	      distfIon:copy(distfDupIon)
	      hamilElc:copy(hamilDupElc)
	      hamilIon:copy(hamilDupIon)
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

-- write data to H5 files
function writeFields(frameNum, tCurr)
   numDensityElc:write( string.format("numDensityElc_%d.h5", frameNum), tCurr)
   numDensityIon:write( string.format("numDensityIon_%d.h5", frameNum), tCurr)
   --distfElc:write( string.format("distfElc_%d.h5", frameNum), tCurr)
   --distfIon:write( string.format("distfIon_%d.h5", frameNum), tCurr)
   phi1dDg:write( string.format("phi_%d.h5", frameNum), tCurr)
   heatFluxAtEdge:write( string.format("heatFluxAtEdge_%d.h5", frameNum), tCurr)
   --cutoffVelocities:write( string.format("cutoffV_%d.h5", frameNum), tCurr)
   --sheathCoefficients:write( string.format("sheathCoefficients_%d.h5", frameNum) ,tCurr)
end


-- unit density ion distribution function
distfIonUnitTest = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}
distfIonUnitTestHigh = DataStruct.Field3D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}

-- unit density ion distribution function
distfElcUnitTest = DataStruct.Field3D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}

initDistfIonTest = Updater.ProjectOnNodalBasis3D {
   onGrid = gridIon,
   basis = basisIon,
   shareCommonNodes = false,
   evaluate = function(z,vPara,mu,t)
     local backgroundTemp = 100 -- ev
     local nHat = 2
     local vTi = math.sqrt(backgroundTemp*eV/ionMass)

     if z > lSource/2 then
       return nHat*maxwellianRight(ionMass, vTi, vPara, mu)
     elseif z < -lSource/2 then
       return nHat*maxwellianLeft(ionMass, vTi, vPara, mu)
     else
       -- Must be between the source boundaries, use a linear combo.
       return ((lSource/2 + z)*nHat*maxwellianRight(ionMass, vTi, vPara, mu) + 
          (lSource/2 - z)*nHat*maxwellianLeft(ionMass, vTi, vPara, mu))/lSource
     end
	 end
}
runUpdater(initDistfIonTest, 0.0, 0.0, {}, {distfIonUnitTest})

initDistfIonTestHigh = Updater.ProjectOnNodalBasis3D {
   onGrid = gridIon,
   basis = basisIon,
   shareCommonNodes = false,
   evaluate = function(z,vPara,mu,t)
     local backgroundTemp = 500 -- ev
     local nHat = 2
     local vTi = math.sqrt(backgroundTemp*eV/ionMass)

     if z > lSource/2 then
       return nHat*maxwellianRight(ionMass, vTi, vPara, mu)
     elseif z < -lSource/2 then
       return nHat*maxwellianLeft(ionMass, vTi, vPara, mu)
     else
       -- Must be between the source boundaries, use a linear combo.
       return ((lSource/2 + z)*nHat*maxwellianRight(ionMass, vTi, vPara, mu) + 
          (lSource/2 - z)*nHat*maxwellianLeft(ionMass, vTi, vPara, mu))/lSource
     end
	 end
}
runUpdater(initDistfIonTestHigh, 0.0, 0.0, {}, {distfIonUnitTestHigh})

initDistfElcTest = Updater.ProjectOnNodalBasis3D {
   onGrid = gridElc,
   basis = basisElc,
   shareCommonNodes = false,
   evaluate = function(z,vPara,mu,t)
     local backgroundTemp = 500 -- ev
     local vTe = math.sqrt(backgroundTemp*eV/elcMass)
     return maxwellian(elcMass, vTe, vPara, mu)
	 end
}
runUpdater(initDistfElcTest, 0.0, 0.0, {}, {distfElcUnitTest})


--applyBc(0.0, 0.0, distfElcUnitTest, distfIonUnitTest, cutoffVelocities)
runUpdater(bcLowerIon, 0.0, 0.0, {}, {distfIonUnitTest})
runUpdater(bcUpperIon, 0.0, 0.0, {}, {distfIonUnitTest})
calcMoments(0.0, 0.0, distfElcUnitTest, distfIonUnitTest)
numDensityElc:write( string.format("numDensityElc_%d.h5", 0), 0)
-- Use reflecting BC for the electrons
--runUpdater(reflectingBc, 0.0, 0.0, {mom0Ion, mom1ParaIon}, {distfElcUnitTest, cutoffVelocities})
--applyBc(0.0, 0.0, distfIonUnitTest, distfIonUnitTest, cutoffVelocities)
--
--
reflectingBcIonTest = Updater.SOL3DElectrostaticDistFuncReflectionBc {
   onGrid = gridIon,
   basis = basisIon,
   basis1d = basis_1d,
   edge = "both",
   -- physical parameters
   scaleFactor = 2*math.pi*B0/ionMass,
}

runUpdater(reflectingBc, 0.0, 0.0, {mom1ParaIon}, {distfElcUnitTest, cutoffVelocities})

