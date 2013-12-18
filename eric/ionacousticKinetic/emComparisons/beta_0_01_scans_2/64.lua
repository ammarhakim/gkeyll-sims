-- Input file for ion acoustic problem with kinetic ions and electrons
-- Electromagnetic terms added
-- Basic test case.. using real units
-- Uses constant density for ions and perturbed density for electrons

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.1

-- wave-number
knumber = 0.5

-- initial number density in each cell (1/m^3)
-- Corresponds to beta_e of 0.01
initNumDens = 9.947e19
-- temperature ratio (T_i/T_e)
Tratio = 0.25

kPerpTimesRho = 0.01
-- electron temperature (eV)
elcTemp = 250
-- ion temperature (eV)
ionTemp = elcTemp*Tratio
-- signed electron charge
elcCharge = -Lucee.ElementaryCharge
-- signed ion charge
ionCharge = Lucee.ElementaryCharge
-- electron mass
elcMass = Lucee.ElectronMass
-- ion mass
ionMass = Lucee.ProtonMass

-- permittivity of free space
epsilon0 = Lucee.Epsilon0
-- permeability of free space
mu0 = Lucee.Mu0
-- Magnetic field (Tesla)
B0 = 1.0

-- Derived parameters
-- Ion cyclotron frequency
Omega_i = ionCharge*B0/ionMass
-- Ion sound speed
c_i = math.sqrt(elcTemp*Lucee.ElementaryCharge/ionMass)
-- Sound-based gyroradius
rho_i = c_i/Omega_i
-- k_perpendicular
kPerp = kPerpTimesRho/rho_i

-- domain extents
XL, XU = -Lucee.Pi/knumber, Lucee.Pi/knumber
-- number of cells
NX, NP = 16, 64
-- compute max thermal speed to set velocity space extents
vtElc = math.sqrt(elcTemp*Lucee.ElementaryCharge/elcMass)
PL_ELC, PU_ELC = -6.0*elcMass*vtElc, 6.0*elcMass*vtElc

vtIon = math.sqrt(ionTemp*Lucee.ElementaryCharge/ionMass)
PL_ION, PU_ION = -6.0*ionMass*vtIon, 6.0*ionMass*vtIon

-- parameters to control time-stepping
tStart = 0.0
tEnd = 3e-5
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

-- updater to initialize distribution function
initDistfElc = Updater.EvalOnNodes2D {
   onGrid = gridElc,
   -- basis functions to use
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
      local elcThermal = math.sqrt(elcTemp*Lucee.ElementaryCharge/elcMass)
      local alpha = 1e-6 -- perturbation
		  local k = knumber
		  local nHat = initNumDens*(1+alpha*math.cos(k*x))
		  return maxwellian(nHat, elcMass, elcThermal, y)
	   end
}
runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})

-- updater to initialize distribution function
initDistfIon = Updater.EvalOnNodes2D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		 local ionThermal = math.sqrt(ionTemp*Lucee.ElementaryCharge/ionMass)
		 return maxwellian(initNumDens, ionMass, ionThermal, y)
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
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distf1Ion = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisElc:numNodes(),
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

-- Updater for electron Vlasov equation
vlasovSolverElc = Updater.PoissonBracket {
   onGrid = gridElc,
   basis = basisElc,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
}
vlasovSolverIon = Updater.PoissonBracket {
   onGrid = gridIon,
   basis = basisIon,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
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

-- Electron number density
numDensityElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Ion number density
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

-- field to store continous potential in 1D
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

-- updater to copy 1D (CG) field to 2D (CG) field
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
   return runUpdater(copier, curr, dt, {phi1}, {phi2})
end

-- update Vlasov equation, given appropriate updater 
function updateVlasovEqn(pbSlvr, curr, dt, distfIn, hamilIn, distfOut)
   return runUpdater(pbSlvr, curr, dt, {distfIn, hamilIn}, {distfOut})
end

-- updater to compute phi electrostatically
electrostaticPhiCalc = Updater.ElectrostaticPhiUpdater {
  -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   kPerpTimesRho = kPerpTimesRho,
   Te0 = elcTemp,
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

-- apply boundary conditions
function applyBc(curr, dt, fldElc, fldIon)
   for i,fld in ipairs({fldElc, fldIon}) do
      fld:applyPeriodicBc(0)
      fld:applyCopyBc(1, "lower")
      fld:applyCopyBc(1, "upper")
   end
end

-- dynvector for field energy
fieldEnergy = DataStruct.DynVector { numComponents = 1, }

-- to compute field energy
fieldEnergyCalc = Updater.NormGrad1D {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
}

totalElcEnergy = DataStruct.DynVector { numComponents = 1, }
totalIonEnergy = DataStruct.DynVector { numComponents = 1, }
-- to compute total energy
totalElcEnergyCalc = Updater.KineticEnergyUpdater {
   -- grid for updater
   onGrid = gridElc,
   -- basis functions to use
   basis = basisElc,
}
totalIonEnergyCalc = Updater.KineticEnergyUpdater {
   -- grid for updater
   onGrid = gridIon,
   -- basis functions to use
   basis = basisIon,
}

-- compute various diagnostics
function calcDiagnostics(curr, dt)
  runUpdater(fieldEnergyCalc, curr, dt, {phi1d}, {fieldEnergy})
  runUpdater(totalElcEnergyCalc, curr, dt, {distfElc, hamilElc}, {totalElcEnergy})
  runUpdater(totalIonEnergyCalc, curr, dt, {distfIon, hamilIon}, {totalIonEnergy})
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local statusElc, dtSuggestedElc
   local statusIon, dtSuggestedIon

   -- RK stage 1
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distfElc, hamilElc, distf1Elc)
   distf1Ion:copy(distfIon)
   if (statusElc == false) then
      return false, dtSuggestedElc
   end
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)

   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)
   -- Compute phi
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon}, {phi1dDg})
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
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distf1Elc, hamilElc, distfNewElc)
   distfNewIon:copy(distf1Ion)
   if (statusElc == false) then
      return false, dtSuggestedElc
   end
   distf1Elc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   distf1Ion:combine(3.0/4.0, distfIon, 1.0/4.0, distfNewIon)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)

   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)

   -- Compute phi
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon}, {phi1dDg})
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
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distf1Elc, hamilElc, distfNewElc)
   distfNewIon:copy(distf1Ion)
   if (statusElc == false) then
      return false, dtSuggestedElc
   end
   distf1Elc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   distf1Ion:combine(1.0/3.0, distfIon, 2.0/3.0, distfNewIon)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)

   calcMoments(tCurr, myDt, distf1Elc, distf1Ion)

   -- Compute phi
   runUpdater(electrostaticPhiCalc, tCurr, myDt, {numDensityElc, numDensityIon}, {phi1dDg})
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
   
   distfElc:copy(distf1Elc)
   distfIon:copy(distf1Ion)

   return true, dtSuggestedElc
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

-- create updater to copy a continuous field to a discontinuous field
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
  --distfElc:write(string.format("distfElc_%d.h5", frameNum), tCurr)
  --distfIon:write(string.format("distfIon_%d.h5", frameNum), tCurr)

  --numDensityElc:write(string.format("numDensityElc_%d.h5", frameNum), tCurr)
  --numDensityIon:write(string.format("numDensityIon_%d.h5", frameNum), tCurr)

  --calcDiscontinuousField(0.0, 0.0, phi1d, phi1dDg)
  --phi1dDg:write(string.format("phi_%d.h5", frameNum), tCurr)

  --aParallel1dDg:write(string.format("a_%d.h5", frameNum), tCurr)

  fieldEnergy:write( string.format("fieldEnergy_%d.h5", frameNum) )

  --totalElcEnergy:write( string.format("totalElcEnergy_%d.h5", frameNum) )
  --totalIonEnergy:write( string.format("totalIonEnergy_%d.h5", frameNum) )
end

applyBc(0.0, 0.0, distfElc, distfIon)
-- calculate initial potential
calcMoments(0.0, 0.0, distfElc, distfIon)
runUpdater(electrostaticPhiCalc, 0.0, 0.0, {numDensityElc, numDensityIon}, {phi1dDg})
calcContinuous1dField(0.0, 0.0, phi1dDg, phi1d)
-- Compute A-parallel
runUpdater(electromagneticACalc, 0.0, 0.0, {mom0Elc, mom0Ion, mom1Elc, mom1Ion}, {aParallel1dDg})
calcContinuous1dField(0.0, 0.0, aParallel1dDg, aParallel1d)
-- Copy continuous field back to discontinuous field
calcDiscontinuousField(0.0, 0.0, aParallel1d, aParallel1dDg)
-- Compute projection of A(x)^2
runUpdater(aSquaredCalc, 0.0, 0.0, {aParallel1dDg}, {aSquared1dDg})
calcContinuous1dField(0.0, 0.0, aSquared1dDg, aSquared1d)

-- calculate initial hamiltonian
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
