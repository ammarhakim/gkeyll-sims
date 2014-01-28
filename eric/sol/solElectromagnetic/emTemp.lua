-- Input file for SOL problem with kinetic ions and electrons
-- This test has positivty preservation turned off
-- PolyOrder 2 and kPerp = 0.2

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.1

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

kPerpTimesRho = 0.1
-- pedestal density (1/m^3)
nPed = 5e19
-- pedestal (source) temperature (eV)
tPed = 1500
-- Fixed value of Te that must be independent of time (eV)
Te0 = 250
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
tEnd = 350e-6
nFrames = 5

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
   ghost = {1, 1},
}
-- distribution function for ions
distfIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
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

-- Return initial ne(x) in 1/m^3
function initialElectronDens(x)
  local backgroundDens = 0.7 + 0.3*(1 - math.abs(x/lParallel))
  if math.abs(x) < lSource/2 then
       backgroundDens = backgroundDens + 0.5*math.cos(math.pi*x/lSource)
     end
  backgroundDens = backgroundDens*1e19
  return backgroundDens
end

-- Return initial Te(x) in eV
function initialElectronTemp(x)
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
     -- Average background density
     local n0 = (0.7 + 0.3/2 + 0.5*lSource/(lParallel*math.pi))*1e19
     local nHat = (initialElectronDens(x) + kPerpTimesRho^2*n0)/(1 + kPerpTimesRho^2)
     local vTe = math.sqrt(initialElectronTemp(x)*eV/elcMass)
     return maxwellian(nHat, elcMass, vTe, y)
	 end
}
runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})

-- Return initial ne(x) in 1/m^3
function initialIonDens(x)
  local backgroundDens = 0.7 + 0.3*(1 - math.abs(x/lParallel))
  if math.abs(x) < lSource/2 then
       backgroundDens = backgroundDens + 0.5*math.cos(math.pi*x/lSource)
     end
  backgroundDens = backgroundDens*1e19
  return backgroundDens
end

-- Return initial Te(x) in eV
function initialIonTemp(x)
  local backgroundTemp = 100 + 45*(1 - math.abs(x/lParallel))
  if math.abs(x) < lSource/2 then
       backgroundTemp = backgroundTemp + 30*math.cos(math.pi*x/lSource)
     end
  return backgroundTemp
end

-- updater to initialize distribution function
initDistfIon = Updater.ProjectOnNodalBasis2D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
     local backgroundTemp = initialIonTemp(x)
     local backgroundDens = initialIonDens(x)
     local nHat = 2*backgroundDens
     local vTi = math.sqrt(backgroundTemp*eV/ionMass)

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
runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})

-- extra fields for performing RK update
distfNewElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distf1Elc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distfNewIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}
distf1Ion = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}

-- Electron Hamiltonian
hamilElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- Ion Hamiltonian
hamilIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- DG versions of the hamiltonians
hamilElcDg = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
hamilIonDg = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
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
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
hamilPSquaredIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
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
   ghost = {1, 1},
}
hamilPerpIonDg = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
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
   ghost = {1, 1},
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
   ghost = {1, 1},
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
mom1Elc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom1Ion = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom2Elc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
mom2Ion = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
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
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
}

aParallel1dDg = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- Updater to compute aParallel
electromagneticACalc = Updater.ElectromagneticAUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
  kPerp = kPerp,
  elcMass = elcMass,
  ionMass = ionMass,
  elcCharge = elcCharge,
  ionCharge = ionCharge,
  mu0 = mu0,
}

-- field to store continuous potential in 1D
phi1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
}

phi1dDg = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- updater to compute phi electrostatically
electrostaticPhiCalc = Updater.ElectrostaticPhiUpdater {
  -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   kPerpTimesRho = kPerpTimesRho,
   Te0 = Te0,
   useCutoffVelocities = true,
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
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
hamilKeIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- DG versions of the KE part of the hamiltonians
hamilKeElcDg = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
hamilKeIonDg = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}

-- dynvector for heat flux at edge
heatFluxAtEdge = DataStruct.DynVector { numComponents = 6, }

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
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
aParallel2dIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- CG fields to store A(x)^2 on 2D grids
aSquared2dElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
aSquared2dIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- CG fields to store A(x)*p on 2D grids
aTimesPElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
aTimesPIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
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
   ghost = {1, 1},
}
-- field to store continuous projection of A(x)^2
aSquared1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
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

   -- Compute projection of A(x)^2 onto 2D grid
   copyCont1DTo2D(copyTo2DElc, curr, dt, aSquared1dIn, aSquared2dElc)
   -- Accumulate q^2/2m*A^2 to KE hamiltonian
   hamilKeOut:accumulate(0.5*elcCharge^2/elcMass, aSquared2dElc)

   -- Compute projection of aParallel onto 2D grid
   copyCont1DTo2D(copyTo2DElc, curr, dt, aParallel1dIn, aParallel2dElc)
   -- Compute A(x)*p on a 2D grid
   runUpdater(aTimesPElcCalc, curr, dt, {aParallel2dElc}, {aTimesPElc})
   -- Accumulate -q/m*p*A to KE hamiltonian
   hamilKeOut:accumulate(-elcCharge/elcMass, aTimesPElc)

   -- Accumulate the KE hamiltonian to the full hamiltonian
   hamilOut:accumulate(1.0, hamilKeOut)
end
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

   -- Compute projection of A(x)^2 onto 2D grid
   copyCont1DTo2D(copyTo2DIon, curr, dt, aSquared1dIn, aSquared2dIon)
   -- Accumulate q^2/2m*A^2 to hamiltonian
   hamilKeOut:accumulate(0.5*ionCharge^2/ionMass, aSquared2dIon)

   -- Compute projection of aParallel onto 2D grid
   copyCont1DTo2D(copyTo2DIon, curr, dt, aParallel1dIn, aParallel2dIon)
   -- Compute A(x)*p on a 2D grid
   runUpdater(aTimesPIonCalc, curr, dt, {aParallel2dIon}, {aTimesPIon})
   -- Accumulate -q/m*p*A to hamiltonian
   hamilKeOut:accumulate(-ionCharge/ionMass, aTimesPIon)

   -- Accumulate the KE hamiltonian to the full hamiltonian
   hamilOut:accumulate(1.0, hamilKeOut)
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
   calcDiscontinuousField(0.0, 0.0, phi1d, phi1dDg)

   -- take the KE part of the hamiltonian and copy it to a DG field
   runUpdater(copyCToDElc2D, 0.0, 0.0, {hamilKeElc}, {hamilKeElcDg})
   runUpdater(copyCToDIon2D, 0.0, 0.0, {hamilKeIon}, {hamilKeIonDg})
   -- compute moments at edges
   runUpdater(momentsAtEdgesElcCalc, curr, dt, {distfElc, hamilKeElcDg, aParallel1dDg}, {momentsAtEdgesElc})
   runUpdater(momentsAtEdgesIonCalc, curr, dt, {distfIon, hamilKeIonDg, aParallel1dDg}, {momentsAtEdgesIon})
   -- compute heat flux at edges
   runUpdater(heatFluxAtEdgeCalc, curr, dt, {phi1dDg, momentsAtEdgesElc, momentsAtEdgesIon}, {heatFluxAtEdge})

   -- copy hamiltonian to DG field
   --runUpdater(copyCToDElc2D, curr, dt, {hamilElc}, {hamilElcDg})
   --runUpdater(copyCToDIon2D, curr, dt, {hamilIon}, {hamilIonDg})
   -- accumulate perpendicular energy terms
   --hamilElcDg:accumulate(1.0, hamilPerpElcDg)
   --hamilIonDg:accumulate(1.0, hamilPerpIonDg)
   -- compute energy using discrete hamiltonian
   --runUpdater(hamilElcEnergyCalc, curr, dt, {distfElc, hamilElcDg}, {hamilElcEnergy})
   --runUpdater(hamilIonEnergyCalc, curr, dt, {distfIon, hamilIonDg}, {hamilIonEnergy})
   -- need to divide by mass to get correct terms!!
   --
   -- compute source energy using discrete hamiltonian
   --runUpdater(hamilElcEnergyCalc, curr, dt, {particleSourceElc, hamilElcDg}, {hamilSrcElcEnergy})
   --runUpdater(hamilIonEnergyCalc, curr, dt, {particleSourceIon, hamilIonDg}, {hamilSrcIonEnergy})
   
   --runUpdater(totalEnergyCalc, curr, dt, {heatFluxAtEdge, hamilElcEnergy, 
   -- hamilIonEnergy, hamilSrcElcEnergy, hamilSrcIonEnergy}, {totalEnergy})
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
   -- Compute A-parallel
   runUpdater(electromagneticACalc, tCurr, myDt, {mom0Elc, mom0Ion, mom1Elc, mom1Ion}, {aParallel1dDg})
   calcContinuous1dField(tCurr, myDt, aParallel1dDg, aParallel1d)
   -- Copy continuous field back to discontinuous field
   calcDiscontinuousField(tCurr, myDt, aParallel1d, aParallel1dDg)

   applyBc(tCurr, myDt, distf1Elc, distf1Ion, aParallel1dDg, cutoffVelocities1)
   
   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon, cutoffVelocities1}, {phi1dDg})
   calcContinuous1dField(tCurr, myDt, phi1dDg, phi1d)
   -- Compute A-parallel
   runUpdater(electromagneticACalc, tCurr, myDt, {mom0Elc, mom0Ion, mom1Elc, mom1Ion}, {aParallel1dDg})
   calcContinuous1dField(tCurr, myDt, aParallel1dDg, aParallel1d)
   -- Copy continuous field back to discontinuous field
   calcDiscontinuousField(tCurr, myDt, aParallel1d, aParallel1dDg)
   -- Compute projection of A(x)^2
   runUpdater(aSquaredCalc, tCurr, myDt, {aParallel1dDg}, {aSquared1dDg})
   calcContinuous1dField(tCurr, myDt, aSquared1dDg, aSquared1d)
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

   -- Compute A-parallel
   runUpdater(electromagneticACalc, tCurr, myDt, {mom0Elc, mom0Ion, mom1Elc, mom1Ion}, {aParallel1dDg})
   calcContinuous1dField(tCurr, myDt, aParallel1dDg, aParallel1d)
   -- Copy continuous field back to discontinuous field
   calcDiscontinuousField(tCurr, myDt, aParallel1d, aParallel1dDg)

   applyBc(tCurr, myDt, distf1Elc, distf1Ion, aParallel1dDg, cutoffVelocities2)
   
   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon, cutoffVelocities2}, {phi1dDg})
   calcContinuous1dField(tCurr, myDt, phi1dDg, phi1d)
   -- Compute A-parallel
   runUpdater(electromagneticACalc, tCurr, myDt, {mom0Elc, mom0Ion, mom1Elc, mom1Ion}, {aParallel1dDg})
   calcContinuous1dField(tCurr, myDt, aParallel1dDg, aParallel1d)
   -- Copy continuous field back to discontinuous field
   calcDiscontinuousField(tCurr, myDt, aParallel1d, aParallel1dDg)
   -- Compute projection of A(x)^2
   runUpdater(aSquaredCalc, tCurr, myDt, {aParallel1dDg}, {aSquared1dDg})
   calcContinuous1dField(tCurr, myDt, aSquared1dDg, aSquared1d)
   -- Compute Hamiltonian
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
   -- Compute A-parallel
   runUpdater(electromagneticACalc, tCurr, myDt, {mom0Elc, mom0Ion, mom1Elc, mom1Ion}, {aParallel1dDg})
   calcContinuous1dField(tCurr, myDt, aParallel1dDg, aParallel1d)
   -- Copy continuous field back to discontinuous field
   calcDiscontinuousField(tCurr, myDt, aParallel1d, aParallel1dDg)

   applyBc(tCurr, myDt, distf1Elc, distf1Ion, aParallel1dDg, cutoffVelocities)

   distfElc:copy(distf1Elc)
   distfIon:copy(distf1Ion)

   calcMoments(tCurr, myDt, distfElc, distfIon)
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon, cutoffVelocities}, {phi1dDg})
   calcContinuous1dField(tCurr, myDt, phi1dDg, phi1d)
   -- Compute A-parallel
   runUpdater(electromagneticACalc, tCurr, myDt, {mom0Elc, mom0Ion, mom1Elc, mom1Ion}, {aParallel1dDg})
   calcContinuous1dField(tCurr, myDt, aParallel1dDg, aParallel1d)
   -- Copy continuous field back to discontinuous field
   calcDiscontinuousField(tCurr, myDt, aParallel1d, aParallel1dDg)
   -- Compute projection of A(x)^2
   runUpdater(aSquaredCalc, tCurr, myDt, {aParallel1dDg}, {aSquared1dDg})
   calcContinuous1dField(tCurr, myDt, aSquared1dDg, aSquared1d)
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
   --numDensityElc:write( string.format("numDensityElc_%d.h5", frameNum), tCurr)
   --numDensityIon:write( string.format("numDensityIon_%d.h5", frameNum), tCurr)
   distfElc:write( string.format("distfElc_%d.h5", frameNum), tCurr)
   --distfIon:write( string.format("distfIon_%d.h5", frameNum), tCurr)
   --phi1dDg:write( string.format("phi_%d.h5", frameNum), tCurr)
   heatFluxAtEdge:write( string.format("heatFluxAtEdge_%d.h5", frameNum) ,tCurr)
   --cutoffVelocities:write( string.format("cutoffV_%d.h5", frameNum) )
   --totalEnergy:write( string.format("totalEnergy_%d.h5", frameNum) ,tCurr)
   --momentumIon:write( string.format("mom1Ion_%d.h5", frameNum), tCurr)
   --momentumElc:write( string.format("mom1Elc_%d.h5", frameNum), tCurr)
   --mom3Ion:write( string.format("mom3Ion_%d.h5", frameNum), tCurr)
   --mom3Elc:write( string.format("mom3Elc_%d.h5", frameNum), tCurr)
   --calcDiscontinuousField(0.0, 0.0, phi1d, phi1dDg)
   --phi1dDg:write(string.format("phi_%d.h5", frameNum), tCurr)
   --aParallel1dDg:write(string.format("a_%d.h5", frameNum), tCurr)
end

calcMoments(0.0, 0.0, distfElc, distfIon)

-- set initial a_parallel to zero
aParallel1dDg:clear(0.0)
aParallel1d:clear(0.0)
-- set initial a^2 to zero
aSquared1d:clear(0.0)

applyBc(0.0, 0.0, distfElc, distfIon, aParallel1dDg, cutoffVelocities)

-- calculate initial phi
calcMoments(0.0, 0.0, distElc, distIon)
runUpdater(electrostaticPhiCalc, 0.0, 0.0, {numDensityElc, numDensityIon, cutoffVelocities}, {phi1dDg})
calcContinuous1dField(0.0, 0.0, phi1dDg, phi1d)
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
