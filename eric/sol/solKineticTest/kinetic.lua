-- Input file for SOL problem with kinetic ions and electrons

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.01

-- physical constants
-- eletron mass (kg)
electronMass = Lucee.ElectronMass
-- electron volt to joules
eV = Lucee.ElementaryCharge
-- Deuterium ion mass (kg)
ionMass = 2.014*Lucee.ProtonMass
-- signed ion charge
ionCharge = Lucee.ElementaryCharge

-- physical input parameters

kPerpTimesRho = 0.2
-- pedestal density (1/m^3)
nPed = 5e19
-- pedestal (source) temperature (eV)
tPed = 1500
-- ELM pulse duration (seconds)
tELM = 200e-6
-- Parallel length (m)
lParallel = 40
-- Source length (m)
lSource = 25
-- Particle source proportionality factor
A = 1.2

-- Derived parameters
-- Thermal velocity of pedestal
vtPed = math.sqrt(tPed*eV/ionMass)
-- Pedestal sound speed (m/s)
cPed = math.sqrt(2*tPed*eV/ionMass)
-- Particle source
Sn   = A*nPed*cPed/lSource

-- domain extents
XL, XU = -lParallel, lParallel
-- number of cells
NX, NV = 32, 16
-- compute max thermal speed to set velocity space extents
vtElc = math.sqrt(tPed*eV/electronMass)
VL_ELC, VU_ELC = -6.0*vtElc, 6.0*vtElc

vtIon = math.sqrt(tPed*eV/ionMass)
VL_ION, VU_ION = -6.0*vtIon, 6.0*vtIon

-- L-B coefficient
lbAlpha = 0

-- parameters to control time-stepping
tStart = 0.0
tEnd = 500e-6
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
elseif (polyOrder == 2) then
   numCgNodesPerCell = 5
   numCgNodesPerCell_1d = 3
end

-- phase space grid for electrons
gridElc = Grid.RectCart2D {
   lower = {XL, VL_ELC},
   upper = {XU, VU_ELC},
   cells = {NX, NV},
}
gridIon = Grid.RectCart2D {
   lower = {XL, VL_ION},
   upper = {XU, VU_ION},
   cells = {NX, NV},
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

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdrift, vt, v)
   return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

-- Maxwellian (Right half only)
function maxwellianRight(n0, vdrift, vt, v)
  if v>=0 then
    return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
  else
    return 0
  end
end

-- Maxwellian (Left half only)
function maxwellianLeft(n0, vdrift, vt, v)
  if v <= 0 then
    return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
  else
    return 0
  end
end

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

-- updater to initialize distribution function
initDistfElc = Updater.EvalOnNodes2D {
   onGrid = gridElc,
   -- basis functions to use
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
     local backgroundTemp = 100 + 45*(1 - math.abs(x/lParallel))
     local backgroundDens = 0.7 + 0.3*(1 - math.abs(x/lParallel))
     if math.abs(x) < lSource/2 then
       backgroundTemp = backgroundTemp + 30*math.cos(math.pi*x/lSource)
       backgroundDens = backgroundDens + 0.5*math.cos(math.pi*x/lSource)
     end
     backgroundDens = backgroundDens*1e19
     
     local nHat = 2*backgroundDens
     local vTe = math.sqrt(backgroundTemp*eV/electronMass)

     return maxwellian(backgroundDens, 0.0, vTe, y)
     --if x > lSource/2 then
     --  return maxwellianRight(nHat, 0.0, vTe, y)
     --elseif x < -lSource/2 then
     --  return maxwellianLeft(nHat, 0.0, vTe, y)
     --else
     --  -- Must be between the source boundaries, use a linear combo.
     --  return ((lSource/2 + x)*maxwellianRight(nHat, 0.0, vTe, y) + (lSource/2 - x)*maxwellianLeft(nHat, 0.0, vTe, y))/lSource
     --end
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
     local backgroundTemp = 100 + 45*(1 - math.abs(x/lParallel))
     local backgroundDens = 0.7 + 0.3*(1 - math.abs(x/lParallel))
     if math.abs(x) < lSource/2 then
       backgroundTemp = backgroundTemp + 30*math.cos(math.pi*x/lSource)
       backgroundDens = backgroundDens + 0.5*math.cos(math.pi*x/lSource)
     end
     backgroundDens = backgroundDens*1e19
     
     local nHat = 2*backgroundDens
     local vTi = math.sqrt(backgroundTemp*eV/ionMass)

     return maxwellian(backgroundDens, 0.0, vTi, y)

     --if x > lSource/2 then
     --  return maxwellianRight(nHat, 0.0, vTi, y)
     --elseif x < -lSource/2 then
     --  return maxwellianLeft(nHat, 0.0, vTi, y)
     --else
       -- Must be between the source boundaries, use a linear combo.
     --  return ((lSource/2 + x)*maxwellianRight(nHat, 0.0, vTi, y) + (lSource/2 - x)*maxwellianLeft(nHat, 0.0, vTi, y))/lSource
     --end
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

-- create field to store kinetic energy term in Hamiltonian
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

-- updater to initialize electron kinetic energy term in Hamiltonian
initHamilKeElc = Updater.EvalOnNodes2D {
   onGrid = gridElc,
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = true, -- Hamiltonian is continuous
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}
runUpdater(initHamilKeElc, 0.0, 0.0, {}, {hamilKeElc})

-- updater to initialize ion kinetic energy term in Hamiltonian
initHamilKeIon = Updater.EvalOnNodes2D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = true, -- Hamiltonian is continuous
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}
runUpdater(initHamilKeIon, 0.0, 0.0, {}, {hamilKeIon})

-- particle source (ELECTRONS)
particleSourceElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
particleSourceElc:clear(0.0)
-- updater to fill out particle source
particleSourceUpdaterElc = Updater.EvalOnNodes2D {
   onGrid = gridElc,
   -- basis functions to use
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
    if math.abs(x) < lSource/2 then
      return maxwellian(Sn*math.cos(math.pi*x/lSource), 0.0, math.sqrt(tPed*eV/electronMass), y)
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
particleSourceUpdaterIon = Updater.EvalOnNodes2D {
   onGrid = gridIon,
   -- basis functions to use
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
    if math.abs(x) < lSource/2 then
      return maxwellian(Sn*math.cos(math.pi*x/lSource), 0.0, math.sqrt(tPed*eV/ionMass), y)
    else
      return 0
    end
	 end
}
particleSourceUpdaterIon:setOut( {particleSourceIon} )
-- initialize
particleSourceUpdaterIon:advance(0.0) -- time is irrelevant

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
-- Electron momentum
momentumElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Ion momentum
momentumIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Electron particle energy
ptclEnergyElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Ion particle energy
ptclEnergyIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Electron energy flux
mom3Elc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Ion energy flux
mom3Ion = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- drift velocity u(x)
driftUElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
driftUIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- thermal velocity squared
vThermSqElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
vThermSqIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- Updater to compute electron number density
numDensityCalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = gridElc,
   basis2d = basisElc,
   basis1d = basis_1d,
   moment = 0,
}
numDensityCalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = gridIon,
   basis2d = basisIon,
   basis1d = basis_1d,
   moment = 0,
}

-- Updater to compute electron momentum
momentumCalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = gridElc,
   basis2d = basisElc,
   basis1d = basis_1d,
   moment = 1,
}
momentumCalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = gridIon,
   basis2d = basisIon,
   basis1d = basis_1d,
   moment = 1,
}

-- Updater to compute electron energy
ptclEnergyCalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = gridElc,
   basis2d = basisElc,
   basis1d = basis_1d,
   moment = 2,
}
ptclEnergyCalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = gridIon,
   basis2d = basisIon,
   basis1d = basis_1d,
   moment = 2,
}

-- to compute third moment
mom3CalcElc = Updater.DistFuncMomentCalc1D {
   onGrid = gridElc,
   basis2d = basisElc,
   basis1d = basis_1d,
   moment = 3,
}
mom3CalcIon = Updater.DistFuncMomentCalc1D {
   onGrid = gridIon,
   basis2d = basisIon,
   basis1d = basis_1d,
   moment = 3,
}

vFromMomentsCalc = Updater.VelocitiesFromMomentsUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
}

-- calculate number density at current time
function calcNumDensity(calculator, curr, dt, distfIn, numDensOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {numDensOut})
end

-- calculate momentum at current time
function calcMomentum(calculator, curr, dt, distfIn, momentumOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {momentumOut})
end

-- calculate energy
function calcPtclEnergy(calculator, curr, dt, distfIn, energyOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {energyOut})
end

-- calculate energy
function calcMom3(calculator, curr, dt, distfIn, mom3Out)
   return runUpdater(calculator, curr, dt, {distfIn}, {mom3Out})
end

-- compute moments from distribution function
function calcMoments(curr, dt, distfElcIn, distfIonIn)
  -- number density
  calcNumDensity(numDensityCalcElc, curr, dt, distfElcIn, numDensityElc)
  calcNumDensity(numDensityCalcIon, curr, dt, distfIonIn, numDensityIon)

  -- momentum
  calcMomentum(momentumCalcElc, curr, dt, distfElcIn, momentumElc)
  calcMomentum(momentumCalcIon, curr, dt, distfIonIn, momentumIon)

  -- energy
  calcPtclEnergy(ptclEnergyCalcElc, curr, dt, distfElcIn, ptclEnergyElc)
  calcPtclEnergy(ptclEnergyCalcIon, curr, dt, distfIonIn, ptclEnergyIon)

  -- third moment
  calcMom3(mom3CalcElc, curr, dt, distfElcIn, mom3Elc)
  calcMom3(mom3CalcIon, curr, dt, distfIonIn, mom3Ion)

  -- Update quantities computed from moments
  runUpdater(vFromMomentsCalc, curr, dt, {numDensityElc, momentumElc, ptclEnergyElc}, {driftUElc, vThermSqElc})
  runUpdater(vFromMomentsCalc, curr, dt, {numDensityIon, momentumIon, ptclEnergyIon}, {driftUIon, vThermSqIon})
end

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

-- updater to compute phi electrostatically
electrostaticPhiCalc = Updater.ElectrostaticPhiUpdater {
  -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   kPerpTimesRho = kPerpTimesRho,
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

-- function to copy 1D field to 2D field
function copyPhi(copier, curr, dt, phi1, phi2)
   return runUpdater(copier, curr, dt, {phi1}, {phi2})
end

-- update Vlasov equation, given appropriate updater 
function updateVlasovEqn(pbSlvr, curr, dt, distfIn, hamilIn, distfOut)
   return runUpdater(pbSlvr, curr, dt, {distfIn, hamilIn}, {distfOut})
end

-- compute phi from number density
function calcPhiFromChargeDensity(curr, dt, distElcIn, distIonIn, phiOut)
   calcMoments(curr, dt, distElcIn, distIonIn)
   runUpdater(electrostaticPhiCalc, curr, dt, {numDensityElc, numDensityIon, vThermSqElc}, {phiOut})
end

-- dynvector for heat flux at edge
heatFluxAtEdge = DataStruct.DynVector { numComponents = 3, }

-- to compute total particle energy
heatFluxAtEdgeCalc = Updater.KineticHeatFluxAtEdgeUpdater {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   ionMass = ionMass,
   -- Perpendicular temperature of ions and electrons
   tPerp = tPed,
}

-- set input field
heatFluxAtEdgeCalc:setIn( {numDensityElc, numDensityIon, mom3Elc, mom3Ion, phi1dDg} )
-- set output dynvector
heatFluxAtEdgeCalc:setOut( {heatFluxAtEdge} )

-- updater to move phi from discontinuous to continuous field
phiToContCalc = Updater.ContFromDisCont1D {
   onGrid = grid_1d,
   basis = basis_1d,
}
-- compute continuous phi from discontinuous phi
function calcContinuousPhi(curr, dt, phiIn, phiOut)
   phiOut:clear(0.0)
   return runUpdater(phiToContCalc, curr, dt, {phiIn}, {phiOut})
end

-- compute hamiltonian for electrons
function calcHamiltonianElc(curr, dt, phiIn, hamilOut)
   hamilOut:clear(0.0)
   copyPhi(copyTo2DElc, curr, dt, phiIn, hamilOut)
   hamilOut:scale(Lucee.ElementaryCharge/electronMass)
   hamilOut:accumulate(1.0, hamilKeElc)
end
-- compute hamiltonian for ions
function calcHamiltonianIon(curr, dt, phiIn, hamilOut)
   hamilOut:clear(0.0)
   copyPhi(copyTo2DIon, curr, dt, phiIn, hamilOut)
   hamilOut:scale(ionCharge/ionMass)
   hamilOut:accumulate(1.0, hamilKeIon)
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

-- function to make make BC updaters
function makeBcObjElc()
   local bcLower = Updater.Bc2D {
      onGrid = gridElc,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "lower",
   }
   local bcUpper = Updater.Bc2D {
      onGrid = gridElc,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "upper",
   }
   return bcLower, bcUpper
end

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
bcLowerElc, bcUpperElc = makeBcObjElc()
bcLowerIon, bcUpperIon = makeBcObjIon()

-- apply boundary conditions
function applyBc(curr, dt, fldElc, fldIon)
   --for i,bc in ipairs({bcLowerElc, bcUpperElc}) do
   --   runUpdater(bc, curr, dt, {}, {fldElc})
   --end
   --for i,bc in ipairs({bcLowerIon, bcUpperIon}) do
   --   runUpdater(bc, curr, dt, {}, {fldIon})
   --end

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
-- set input field
fieldEnergyCalc:setIn( {phi1d} )
-- set output dynvector
fieldEnergyCalc:setOut( {fieldEnergy} )

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   heatFluxAtEdgeCalc:setCurrTime(curr)
   heatFluxAtEdgeCalc:advance(curr+dt)

   fieldEnergyCalc:setCurrTime(curr)
   fieldEnergyCalc:advance(curr+dt)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local statusElc, dtSuggestedElc
   local statusIon, dtSuggestedIon

   -- RK stage 1
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distfElc, hamilElc, distf1Elc)
   statusIon, dtSuggestedIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, distfIon, hamilIon, distf1Ion)
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)
   
   if tCurr < tELM then
    distf1Elc:accumulate(myDt, particleSourceElc)
    distf1Ion:accumulate(myDt, particleSourceIon)
   end
   
   calcPhiFromChargeDensity(tCurr, myDt, distf1Elc, distf1Ion, phi1dDg)
   calcContinuousPhi(tCurr, myDt, phi1dDg, phi1d)
   calcHamiltonianElc(tCurr, myDt, phi1d, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, hamilIon)

   -- RK stage 2
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distf1Elc, hamilElc, distfNewElc)
   statusIon, dtSuggestedIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, distf1Ion, hamilIon, distfNewIon)
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   
   if tCurr < tELM then
    distfNewElc:accumulate(myDt, particleSourceElc)
    distfNewIon:accumulate(myDt, particleSourceIon)
   end
   
   distf1Elc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   distf1Ion:combine(3.0/4.0, distfIon, 1.0/4.0, distfNewIon)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)
   
   calcPhiFromChargeDensity(tCurr, myDt, distf1Elc, distf1Ion, phi1dDg)
   calcContinuousPhi(tCurr, myDt, phi1dDg, phi1d)
   calcHamiltonianElc(tCurr, myDt, phi1d, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, hamilIon)

   -- RK stage 3
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distf1Elc, hamilElc, distfNewElc)
   statusIon, dtSuggestedIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, distf1Ion, hamilIon, distfNewIon)
   
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end

   if tCurr < tELM then
    distfNewElc:accumulate(myDt, particleSourceElc)
    distfNewIon:accumulate(myDt, particleSourceIon)
   end

   distf1Elc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   distf1Ion:combine(1.0/3.0, distfIon, 2.0/3.0, distfNewIon)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)

   distfElc:copy(distf1Elc)
   distfIon:copy(distf1Ion)

   calcPhiFromChargeDensity(tCurr, myDt, distfElc, distfIon, phi1dDg)
   calcContinuousPhi(tCurr, myDt, phi1dDg, phi1d)
   calcHamiltonianElc(tCurr, myDt, phi1d, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, hamilIon)

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

      if (tCurr+myDt > tEnd) then
	      myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	      distfElc:copy(distfDupElc)
	      distfIon:copy(distfDupIon)
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
   phi1dDg:write( string.format("phi_%d.h5", frameNum), tCurr)
   heatFluxAtEdge:write( string.format("heatFluxAtEdge_%d.h5", frameNum) )
end

applyBc(0.0, 0.0, distfElc, distfIon)
-- calculate initial Hamiltonian
calcPhiFromChargeDensity(0.0, 0.0, distfElc, distfIon, phi1dDg)
calcContinuousPhi(0.0, 0.0, phi1dDg, phi1d)
calcHamiltonianElc(0.0, 0.0, phi1d, hamilElc)
calcHamiltonianIon(0.0, 0.0, phi1d, hamilIon)
-- compute initial diagnostics
calcDiagnostics(0.0, 0.0)
-- write out initial conditions
writeFields(0, 0.0)
-- make a duplicate in case we need it
distfDupElc = distfElc:duplicate()
distfDupIon = distfIon:duplicate()

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
