-- Input file for SOL problem with kinetic ions and electrons
-- Used to plot moment at edge before BC is applied...
-- Modified so that ghost cells are all = {2, 2}

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.0001

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

-- Pedestal sound speed (m/s)
cPed = math.sqrt(2*tPed*eV/ionMass)
-- Particle source
Sn   = A*nPed*cPed/lSource

-- domain extents
XL, XU = -lParallel, lParallel
-- number of cells
NX, NP = 8, 32
-- compute max thermal speed to set velocity space extents
vtElc = math.sqrt(tPed*eV/elcMass)
PL_ELC, PU_ELC = -6.0*elcMass*vtElc, 6.0*elcMass*vtElc

vtIon = math.sqrt(tPed*eV/ionMass)
PL_ION, PU_ION = -6.0*ionMass*vtIon, 6.0*ionMass*vtIon

-- parameters to control time-stepping
tStart = 0.0
tEnd = 1.6e-6
nFrames = 1

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

-- Determine number of global nodes per cell for use in creating CG
-- fields. Note that this looks a bit odd as this not the number of
-- *local* nodes but the number of nodes in each cell to give the
-- correct number of global nodes in fields.
if (polyOrder == 1) then
   numCgNodesPerCell = 1
   numCgNodesPerCell_1d = 1
elseif (polyOrder == 2) then
   numCgNodesPerCell = 3
   numCgNodesPerCell_1d = 2
elseif (polyOrder == 3) then
   numCgNodesPerCell = 5
   numCgNodesPerCell_1d = 3
end

-- phase space grid for electrons
gridElc = Grid.RectCart2D {
   lower = {XL, PL_ELC},
   upper = {XU, PU_ELC},
   cells = {NX, NP},
}
gridIon = Grid.RectCart2D {
   lower = {XL, PL_ION},
   upper = {XU, PU_ION},
   cells = {NX, NP},
}

-- create FEM nodal basis
basisElc = NodalFiniteElement2D.SerendipityElement {
   onGrid = gridElc,
   polyOrder = polyOrder,
}
basisIon = NodalFiniteElement2D.SerendipityElement {
   onGrid = gridIon,
   polyOrder = polyOrder,
}

-- spatial grid
grid_1d = Grid.RectCart1D {
   lower = {XL},
   upper = {XU},
   cells = {NX},
}
-- spatial FEM nodal basis
basis_1d = NodalFiniteElement1D.LagrangeTensor {
   onGrid = grid_1d,
   polyOrder = polyOrder,
}

-- distribution function for electrons
distfElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {2, 2},
}
-- distribution function for ions
distfIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}
-- unit density ion distribution function
distfIonUnit = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}

-- Electron and ion number densitites
numDensityElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}
numDensityIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}
numDensityIon2dDg = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}
-- Moments for electrons and ions
mom0Elc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}
mom0Ion = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}
mom1Elc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}
mom1Ion = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}
mom2Elc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}
mom2Ion = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}

-- Updater to compute electron number density
mom0CalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = gridElc,
   basis2d = basisElc,
   basis1d = basis_1d,
   moment = 0,
}
mom0CalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = gridIon,
   basis2d = basisIon,
   basis1d = basis_1d,
   moment = 0,
}

-- Updater to compute electron momentum
mom1CalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = gridElc,
   basis2d = basisElc,
   basis1d = basis_1d,
   moment = 1,
}
mom1CalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = gridIon,
   basis2d = basisIon,
   basis1d = basis_1d,
   moment = 1,
}

-- Updater to compute electron energy
mom2CalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = gridElc,
   basis2d = basisElc,
   basis1d = basis_1d,
   moment = 2,
}
mom2CalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = gridIon,
   basis2d = basisIon,
   basis1d = basis_1d,
   moment = 2,
}

-- Maxwellian with number density 'n0', species mass 'mass' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, mass, vt, p)
   return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-p^2/(2*mass^2*vt^2))
end

-- Maxwellian (Right half only)
function maxwellianRight(n0, mass, vt, p)
  if p>=0 then
    return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-p^2/(2*mass^2*vt^2))
  else
    return 0
  end
end

-- Maxwellian (Left half only)
function maxwellianLeft(n0, mass, vt, p)
  if p <= 0 then
    return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-p^2/(2*mass^2*vt^2))
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
initDistfElc = Updater.ProjectOnNodalBasis2D {
   onGrid = gridElc,
   -- basis functions to use
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
     local vTe = math.sqrt(initialElectronTemp(x)*eV/elcMass)
     return maxwellian(initialElcDens(x), elcMass, vTe, y)
	 end
}
runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})
-- Calculate initial electron density field
runUpdater(mom0CalcElc, 0.0, 0.0, {distfElc}, {numDensityElc})
numDensityElc:scale(1/elcMass)

-- updater to initialize distribution function
initDistfIon = Updater.ProjectOnNodalBasis2D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
     local backgroundTemp = initialIonTemp(x)
     local nHat = 2
     local vTi = math.sqrt(2*math.pi/(math.pi-1)*backgroundTemp*eV/ionMass)

     if x > lSource/2 then
       return maxwellianRight(nHat, ionMass, vTi, y)
     elseif x < -lSource/2 then
       return maxwellianLeft(nHat, ionMass, vTi, y)
     else
       -- Must be between the source boundaries, use a linear combo.
       return ((lSource/2 + x)*maxwellianRight(nHat, ionMass, vTi, y) + (lSource/2 - x)*maxwellianLeft(nHat, ionMass, vTi, y))/lSource
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

-- Copy numDensityIon to a 2d field
copyTo2DIonDg = Updater.NodalCopyFaceToInteriorUpdater {
  onGrid = gridIon,
  basis1d = basis_1d,
  basis2d = basisIon,
  dir = 1,
  shareCommonNodes = false,
}
runUpdater(copyTo2DIonDg, 0.0, 0.0, {numDensityIon}, {numDensityIon2dDg})

-- Multiply with gaussian to get total distfIon
multiply2dFields = Updater.FieldArithmeticUpdater2D {
  onGrid = gridIon,
  basis = basisIon,
  evaluate = function(An,Bn,t)
    return An*Bn
  end,
}
runUpdater(multiply2dFields, 0.0, 0.0, {numDensityIon2dDg, distfIonUnit}, {distfIon})

-- extra fields for performing RK update
distfNewElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {2, 2},
}
distf1Elc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {2, 2},
}
distfNewIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}
distf1Ion = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}

-- Electron Hamiltonian
hamilElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {2, 2},
}
-- Ion Hamiltonian
hamilIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {2, 2},
}
-- DG versions of the hamiltonians
hamilElcDg = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {2, 2},
}
hamilIonDg = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}

-- create updater to copy a continuous field to a discontinuous field
copyCToDElc2D = Updater.CopyContToDisCont2D {
   onGrid = gridElc,
   basis = basisElc,
}
-- create updater to copy a continuous field to a discontinuous field
copyCToDIon2D = Updater.CopyContToDisCont2D {
   onGrid = gridIon,
   basis = basisIon,
}

-- create field to store p^2 term in Hamiltonian
hamilPSquaredElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {2, 2},
}
hamilPSquaredIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {2, 2},
}

-- updater to initialize electron kinetic energy term in Hamiltonian
initHamilPSquaredElc = Updater.EvalOnNodes2D {
   onGrid = gridElc,
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = true, -- Hamiltonian is continuous
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local p = y
      return p^2
   end
}
runUpdater(initHamilPSquaredElc, 0.0, 0.0, {}, {hamilPSquaredElc})

-- updater to initialize ion kinetic energy term in Hamiltonian
initHamilPSquaredIon = Updater.EvalOnNodes2D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = true, -- Hamiltonian is continuous
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local p = y
      return p^2
   end
}
runUpdater(initHamilPSquaredIon, 0.0, 0.0, {}, {hamilPSquaredIon})

hamilPerpElcDg = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {2, 2},
}
hamilPerpIonDg = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}

-- updater to initialize perpendicular kinetic energy term in the hamiltonian
initHamilPerpElc = Updater.EvalOnNodes2D {
   onGrid = gridElc,
   basis = basisElc,
   shareCommonNodes = false,
   evaluate = function (x,y,z,t)
      return tPed*eV
   end
}
runUpdater(initHamilPerpElc, 0.0, 0.0, {}, {hamilPerpElcDg})
-- updater to initialize perpendicular kinetic energy term in the hamiltonian
initHamilPerpIon = Updater.EvalOnNodes2D {
   onGrid = gridIon,
   basis = basisIon,
   shareCommonNodes = false,
   evaluate = function (x,y,z,t)
      return tPed*eV
   end
}
runUpdater(initHamilPerpIon, 0.0, 0.0, {}, {hamilPerpIonDg})

-- particle source (ELECTRONS)
particleSourceElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {2, 2},
}
-- clear out contents
particleSourceElc:clear(0.0)
-- updater to fill out particle source
particleSourceUpdaterElc = Updater.ProjectOnNodalBasis2D {
   onGrid = gridElc,
   -- basis functions to use
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
    if math.abs(x) < lSource/2 then
      if t < tELM then
        return maxwellian(Sn*math.cos(math.pi*x/lSource), elcMass, math.sqrt(tPed*eV/elcMass), y)
      else
        return maxwellian(Sn/9*math.cos(math.pi*x/lSource), elcMass, math.sqrt(210*eV/elcMass), y)
      end
    else
      return 0
    end
	 end
}
particleSourceUpdaterElc:setOut( {particleSourceElc} )
-- initialize
particleSourceUpdaterElc:advance(0.0) -- time is irrelevant
-- particle source (IONS)
particleSourceIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}
-- clear out contents
particleSourceIon:clear(0.0)
-- updater to fill out particle source
particleSourceUpdaterIon = Updater.ProjectOnNodalBasis2D {
   onGrid = gridIon,
   -- basis functions to use
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
    if math.abs(x) < lSource/2 then
      if t < tELM then
        return maxwellian(Sn*math.cos(math.pi*x/lSource), ionMass, math.sqrt(tPed*eV/ionMass), y)
      else
        return maxwellian(Sn/9*math.cos(math.pi*x/lSource), ionMass, math.sqrt(260*eV/ionMass), y)
      end
    else
      return 0
    end
	 end
}
particleSourceUpdaterIon:setOut( {particleSourceIon} )
-- initialize
particleSourceUpdaterIon:advance(0.0)

-- keeps track of whether to update particlesources again
postELM = false

-- Updater for electron Vlasov equation
vlasovSolverElc = Updater.PoissonBracket {
   onGrid = gridElc,
   basis = basisElc,
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
}
vlasovSolverIon = Updater.PoissonBracket {
   onGrid = gridIon,
   basis = basisIon,
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
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
  runUpdater(mom0CalcElc, curr, dt, {distfElcIn}, {mom0Elc})
  runUpdater(mom0CalcIon, curr, dt, {distfIonIn}, {mom0Ion})
  -- moment 1
  runUpdater(mom1CalcElc, curr, dt, {distfElcIn}, {mom1Elc})
  runUpdater(mom1CalcIon, curr, dt, {distfIonIn}, {mom1Ion})
  -- moment 2
  runUpdater(mom2CalcElc, curr, dt, {distfElcIn}, {mom2Elc})
  runUpdater(mom2CalcIon, curr, dt, {distfIonIn}, {mom2Ion})
  
  -- Use moment 0 fields to calculate number densities
  numDensityElc:clear(0.0)
  numDensityElc:accumulate(1/elcMass, mom0Elc)
  numDensityIon:clear(0.0)
  numDensityIon:accumulate(1/ionMass, mom0Ion)
end

-- field to store continuous vector potential
aParallel1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {2, 2},
}

aParallel1dDg = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}

-- field to store continuous potential in 1D
phi1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {2, 2},
}
-- continuous potential after it has been modified
phi1dAfterBc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {2, 2},
}

phi1dDg = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}

-- updater to compute phi electrostatically
electrostaticPhiCalc = Updater.ElectrostaticContPhiUpdater {
  -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   kPerpTimesRho = kPerpTimesRho,
   Te0 = Te0,
}

setPhiAtBoundaryCalc = Updater.SetPhiAtBoundaryUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
  elcMass = elcMass,
  elcCharge = eV,
}

-- updater to copy 1D field to 2D field
copyTo2DElc = Updater.Copy1DTo2DNodalField {
   -- grid for updater
   onGrid = gridElc,
}
copyTo2DIon = Updater.Copy1DTo2DNodalField {
   -- grid for updater
   onGrid = gridIon,
}

copyTo2DElcDg = Updater.NodalCopyFaceToInteriorUpdater {
  onGrid = gridElc,
  basis1d = basis_1d,
  basis2d = basisElc,
  dir = 1,
  shareCommonNodes = false,
}

-- function to copy 1D field to 2D field
function copyCont1DTo2D(copier, curr, dt, phi1, phi2)
   phi2:clear(0.0)
   return runUpdater(copier, curr, dt, {phi1}, {phi2})
end

-- Dynvectors to store 0-3rd moments at left and right edges
momentsAtEdgesElc = DataStruct.DynVector { numComponents = 8, }
momentsAtEdgesIon = DataStruct.DynVector { numComponents = 8, }
-- Dynvectors to store 0-3rd moments at left and right edges BEFORE APPLYBC
momentsAtEdgesElcBeforeBc = DataStruct.DynVector { numComponents = 8, }
momentsAtEdgesIonBeforeBc = DataStruct.DynVector { numComponents = 8, }

momentsAtEdgesElcCalc = Updater.ElectromagneticMomentsAtEdgesUpdater {
  onGrid = gridElc,
  basis = basisElc,
  speciesMass = elcMass,
  speciesCharge = elcCharge,
}

momentsAtEdgesIonCalc = Updater.ElectromagneticMomentsAtEdgesUpdater {
  onGrid = gridIon,
  basis = basisIon,
  speciesMass = ionMass,
  speciesCharge = ionCharge,
}

-- Fields to store the KE part of the hamiltonians
hamilKeElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {2, 2},
}
hamilKeIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {2, 2},
}
-- DG versions of the KE part of the hamiltonians
hamilKeElcDg = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {2, 2},
}
hamilKeIonDg = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}

-- dynvector for heat flux at edge
heatFluxAtEdge = DataStruct.DynVector { numComponents = 6, }
-- dynvector for sheath power transmission coefficients
sheathCoefficients = DataStruct.DynVector { numComponents = 7, }

-- to compute total particle energy
heatFluxAtEdgeCalc = Updater.KineticHeatFluxAtEdgeUpdater {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   ionMass = ionMass,
   electronMass = elcMass,
   -- Perpendicular temperature of ions and electrons
   tPerp = tPed,
   -- Enable calculation of sheath coefficients
   computeSheathCoefficient = true,
}

-- dynvector for energy computed using discrete hamiltonian
hamilElcEnergy = DataStruct.DynVector { numComponents = 1, }
hamilSrcElcEnergy = DataStruct.DynVector { numComponents = 1, }
-- updater to compute energy using discrete hamiltonian
hamilElcEnergyCalc = Updater.KineticEnergyUpdater {
  onGrid = gridElc,
  basis = basisElc,
}
-- dynvector for energy computed using discrete hamiltonian
hamilIonEnergy = DataStruct.DynVector { numComponents = 1, }
hamilSrcIonEnergy = DataStruct.DynVector { numComponents = 1, }
-- updater to compute energy using discrete hamiltonian
hamilIonEnergyCalc = Updater.KineticEnergyUpdater {
  onGrid = gridIon,
  basis = basisIon,
}

-- dynvector for total energy in system
totalEnergy = DataStruct.DynVector { numComponents = 3, }
-- updater to compute total energy
totalEnergyCalc = Updater.KineticTotalEnergyUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
  ionMass = ionMass,
  electronMass = elcMass,
}

-- updater to move a discontinuous field to a continuous field
contFromDisContCalc = Updater.ContFromDisCont1D {
   onGrid = grid_1d,
   basis = basis_1d,
}
-- compute continuous 1d field from discontinuous 1d field
function calcContinuous1dField(curr, dt, disContIn, contOut)
   contOut:clear(0.0)
   return runUpdater(contFromDisContCalc, curr, dt, {disContIn}, {contOut})
end

-- Updaters to compute A(x)*p on ion and electron grids
aTimesPElcCalc = Updater.ATimesPUpdater {
  onGrid = gridElc,
  basis = basisElc,
}
aTimesPIonCalc = Updater.ATimesPUpdater {
  onGrid = gridIon,
  basis = basisIon,
}
-- CG fields to store A(x) on 2D grids
aParallel2dElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {2, 2},
}
aParallel2dIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {2, 2},
}
-- CG fields to store A(x)^2 on 2D grids
aSquared2dElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {2, 2},
}
aSquared2dIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {2, 2},
}
-- CG fields to store A(x)*p on 2D grids
aTimesPElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {2, 2},
}
aTimesPIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {2, 2},
}

-- Updater to compute A(x)*A(x) on same basis functions as A(x) (DG field)
aSquaredCalc = Updater.ASquaredProjectionUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
}
-- field to store DG projection of A(x)^2
aSquared1dDg = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}
-- field to store continuous projection of A(x)^2
aSquared1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {2, 2},
}

-- compute hamiltonian for electrons
function calcHamiltonianElc(curr, dt, phiIn, aParallel1dIn, aSquared1dIn, hamilKeOut, hamilOut)
   -- clear out fields (is this needed?)
   hamilOut:clear(0.0)
   hamilKeOut:clear(0.0)
   aSquared2dElc:clear(0.0)
   aParallel2dElc:clear(0.0)

   -- Accumulate q*Phi contribution to hamiltonian
   copyCont1DTo2D(copyTo2DElc, curr, dt, phiIn, hamilOut)
   hamilOut:scale(elcCharge)

   -- Accumulate p^2/2m to KE hamiltonian
   hamilKeOut:accumulate(0.5/elcMass, hamilPSquaredElc)

   -- Accumulate the KE hamiltonian to the full hamiltonian
   hamilOut:accumulate(1.0, hamilKeOut)
end

-- to compute second order hamiltonian term
mhdHamiltonianCalc = Updater.MHDHamiltonianUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
}
-- store result of mhdHamiltonianCalc
mhdHamiltonian1dDg = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {2, 2},
}
-- result of mhdHamiltonianCalc as continuous field
mhdHamiltonian1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {2, 2},
}
mhdHamiltonian2d = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {2, 2},
}
mhdHamiltonian2dDg = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {2, 2},
}
-- stores mhd hamiltonian contribution to total energy
mhdHamiltonianEnergy = DataStruct.DynVector { numComponents = 1, }

-- compute hamiltonian for ions
function calcHamiltonianIon(curr, dt, phiIn, aParallel1dIn, aSquared1dIn, hamilKeOut, hamilOut)
   -- clear out fields (is this needed?)
   hamilOut:clear(0.0)
   hamilKeOut:clear(0.0)
   aSquared2dIon:clear(0.0)
   aParallel2dIon:clear(0.0)

   -- Accumulate q*Phi contribution to hamiltonian
   copyCont1DTo2D(copyTo2DIon, curr, dt, phiIn, hamilOut)
   hamilOut:scale(ionCharge)

   -- Accumulate p^2/2m to hamiltonian
   hamilKeOut:accumulate(0.5/ionMass, hamilPSquaredIon)

   -- Accumulate the KE hamiltonian to the full hamiltonian
   hamilOut:accumulate(1.0, hamilKeOut)

   -- compute second order hamiltonian term
   runUpdater(copyCToD, curr, dt, {phiIn}, {phi1dDg})
   runUpdater(mhdHamiltonianCalc, curr, dt, {phi1dDg, numDensityIon}, {mhdHamiltonian1dDg})
   runUpdater(contFromDisContCalc, curr, dt, {mhdHamiltonian1dDg}, {mhdHamiltonian1d})
   runUpdater(copyTo2DIon, curr, dt, {mhdHamiltonian1d}, {mhdHamiltonian2d})
   mhdHamiltonian2d:scale(-0.5*kPerpTimesRho*kPerpTimesRho*Lucee.ElementaryCharge/Te0)
   hamilOut:accumulate(1.0, mhdHamiltonian2d)
end

-- A HACK
function getRepTbl(pOrder, val)
   if pOrder == 1 then
      return {val, val, val, val}
   elseif pOrder == 2 then
      return {val, val, val, val, val, val, val, val}
   elseif pOrder == 3 then
      return {val, val, val, val, val, val, val, val, val, val, val, val}
   end
end
function getCountTbl(pOrder, val)
   if pOrder == 1 then
      return {0, 1, 2, 3}
   elseif pOrder == 2 then
      return {0, 1, 2, 3, 4, 5, 6, 7}
   elseif pOrder == 3 then
      return {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
   end
end

bcConst = BoundaryCondition.Const { 
   components = getCountTbl(polyOrder),
   values = getRepTbl(polyOrder, 0.0),
}

function makeBcObjIon()
   local bcLower = Updater.Bc2D {
      onGrid = gridIon,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "lower",
   }
   local bcUpper = Updater.Bc2D {
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
reflectingBc = Updater.ElectromagneticDistFuncReflectionBc {
   onGrid = gridElc,
   basis = basisElc,
   edge = "both",
   -- physical parameters
   elcCharge = elcCharge,
   ionCharge = ionCharge,
   elcMass = elcMass,
   ionMass = ionMass,
}

-- apply boundary conditions
function applyBc(curr, dt, fldElc, fldIon, aParallel1dDgIn, cutoffV)
   for i,bc in ipairs({bcLowerIon, bcUpperIon}) do
      runUpdater(bc, curr, dt, {}, {fldIon})
   end
   -- Use reflecting BC for the electrons
   runUpdater(reflectingBc, curr, dt, {mom0Ion, mom1Ion, aParallel1dDgIn}, {fldElc, cutoffV})

   for i,fld in ipairs({fldElc, fldIon}) do
      fld:applyCopyBc(1, "lower")
      fld:applyCopyBc(1, "upper")
   end
end

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   -- take the KE part of the hamiltonian and copy it to a DG field
   runUpdater(copyCToDElc2D, 0.0, 0.0, {hamilKeElc}, {hamilKeElcDg})
   runUpdater(copyCToDIon2D, 0.0, 0.0, {hamilKeIon}, {hamilKeIonDg})
   -- compute moments at edges
   runUpdater(momentsAtEdgesElcCalc, curr, dt, {distfElc, hamilKeElcDg, aParallel1dDg}, {momentsAtEdgesElc})
   runUpdater(momentsAtEdgesIonCalc, curr, dt, {distfIon, hamilKeIonDg, aParallel1dDg}, {momentsAtEdgesIon})
   -- modify phi so that is is equal to phi_s at edge
   runUpdater(setPhiAtBoundaryCalc, curr, dt, {phi1d, cutoffVelocities}, {phi1dAfterBc})
   calcDiscontinuousField(curr, dt, phi1dAfterBc, phi1dDg)
   -- compute heat flux at edges
   runUpdater(heatFluxAtEdgeCalc, curr, dt, {phi1dDg, momentsAtEdgesElc, momentsAtEdgesIon}, {heatFluxAtEdge, sheathCoefficients})
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local statusElc, dtSuggestedElc
   local statusIon, dtSuggestedIon

   if (postELM == false) then
     if tCurr + myDt > tELM then
       postELM = true
       particleSourceUpdaterElc:advance(tCurr + myDt)
       particleSourceUpdaterIon:advance(tCurr + myDt)
     end
   end
   
   -- RK stage 1
   statusElc, dtSuggestedElc = runUpdater(vlasovSolverElc, tCurr, myDt, {distfElc, hamilElc}, {distf1Elc})
   statusIon, dtSuggestedIon = runUpdater(vlasovSolverIon, tCurr, myDt, {distfIon, hamilIon}, {distf1Ion})
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   
   distf1Elc:accumulate(myDt, particleSourceElc)
   distf1Ion:accumulate(myDt, particleSourceIon)
   
   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion, aParallel1dDg, cutoffVelocities1)
   
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon}, {phi1d})

   -- Compute Hamiltonian
   calcHamiltonianElc(tCurr, myDt, phi1d, aParallel1d, aSquared1d, hamilKeElc, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, aParallel1d, aSquared1d, hamilKeIon, hamilIon)

   -- RK stage 2
   statusElc, dtSuggestedElc = runUpdater(vlasovSolverElc, tCurr, myDt, {distf1Elc, hamilElc}, {distfNewElc})
   statusIon, dtSuggestedIon = runUpdater(vlasovSolverIon, tCurr, myDt, {distf1Ion, hamilIon}, {distfNewIon})
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   
   distfNewElc:accumulate(myDt, particleSourceElc)
   distfNewIon:accumulate(myDt, particleSourceIon)
   
   distf1Elc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   distf1Ion:combine(3.0/4.0, distfIon, 1.0/4.0, distfNewIon)
   
   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion, aParallel1dDg, cutoffVelocities2)
   
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon}, {phi1d})

   -- Compute HamiltonianmomentsAtEdgesElcBeforeBc
   calcHamiltonianElc(tCurr, myDt, phi1d, aParallel1d, aSquared1d, hamilKeElc, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, aParallel1d, aSquared1d, hamilKeIon, hamilIon)

   -- RK stage 3
   statusElc, dtSuggestedElc = runUpdater(vlasovSolverElc, tCurr, myDt, {distf1Elc, hamilElc}, {distfNewElc})
   statusIon, dtSuggestedIon = runUpdater(vlasovSolverIon, tCurr, myDt, {distf1Ion, hamilIon}, {distfNewIon})
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end

   distfNewElc:accumulate(myDt, particleSourceElc)
   distfNewIon:accumulate(myDt, particleSourceIon)

   distf1Elc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   distf1Ion:combine(1.0/3.0, distfIon, 2.0/3.0, distfNewIon)
   
   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)
   runUpdater(momentsAtEdgesElcCalc, tCurr, myDt, {distf1Elc, hamilKeElcDg, aParallel1dDg}, {momentsAtEdgesElcBeforeBc})
   runUpdater(momentsAtEdgesIonCalc, tCurr, myDt, {distf1Ion, hamilKeIonDg, aParallel1dDg}, {momentsAtEdgesIonBeforeBc})
   applyBc(tCurr, myDt, distf1Elc, distf1Ion, aParallel1dDg, cutoffVelocities)

   distfElc:copy(distf1Elc)
   distfIon:copy(distf1Ion)

   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon}, {phi1d})

   -- Compute Hamiltonian
   calcHamiltonianElc(tCurr, myDt, phi1d, aParallel1d, aSquared1d, hamilKeElc, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, aParallel1d, aSquared1d, hamilKeIon, hamilIon)

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

-- create updater to copy a continuous 1D field to a 1D discontinuous field
copyCToD = Updater.CopyContToDisCont1D {
   onGrid = grid_1d,
   basis = basis_1d,
}

-- Copy a continuous field to a discontinuous field
function calcDiscontinuousField(tCurr, dt, cgIn, dgOut)
   runUpdater(copyCToD, tCurr, dt, {cgIn}, {dgOut})
end

-- write data to H5 files
function writeFields(frameNum, tCurr)
   numDensityElc:write( string.format("numDensityElc_%d.h5", frameNum), tCurr)
   numDensityIon:write( string.format("numDensityIon_%d.h5", frameNum), tCurr)
   distfElc:write( string.format("distfElc_%d.h5", frameNum), tCurr)
   distfIon:write( string.format("distfIon_%d.h5", frameNum), tCurr)
   phi1dDg:write( string.format("phi_%d.h5", frameNum), tCurr)
   heatFluxAtEdge:write( string.format("heatFluxAtEdge_%d.h5", frameNum), tCurr)
   cutoffVelocities:write( string.format("cutoffV_%d.h5", frameNum), tCurr)
   sheathCoefficients:write( string.format("sheathCoefficients_%d.h5", frameNum) ,tCurr)
   momentsAtEdgesElcBeforeBc:write( string.format("momentsAtEdgesElcBeforeBc_%d.h5", frameNum), tCurr)
   momentsAtEdgesIonBeforeBc:write( string.format("momentsAtEdgesIonBeforeBc_%d.h5", frameNum), tCurr)
end

calcMoments(0.0, 0.0, distfElc, distfIon)
-- Compute A-parallel
aParallel1d:clear(0.0)
aSquared1d:clear(0.0)

applyBc(0.0, 0.0, distfElc, distfIon, aParallel1dDg, cutoffVelocities)

-- calculate initial phi
runUpdater(electrostaticPhiCalc, 0.0, 0.0, {numDensityElc, numDensityIon}, {phi1d})
-- calculate initial Hamiltonian
calcHamiltonianElc(0.0, 0.0, phi1d, aParallel1d, aSquared1d, hamilKeElc, hamilElc)
calcHamiltonianIon(0.0, 0.0, phi1d, aParallel1d, aSquared1d, hamilKeIon, hamilIon)
-- compute initial diagnostics
calcDiagnostics(0.0, 0.0)
-- write out initial conditions
writeFields(0, 0.0)
-- make a duplicate in case we need it
distfDupElc = distfElc:duplicate()
distfDupIon = distfIon:duplicate()
hamilDupElc = hamilElc:duplicate()
hamilDupIon = hamilIon:duplicate()

-- parameters to control time-stepping
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
tFrame = (tEnd-tStart)/nFrames
tCurr = tStart

for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame, tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
