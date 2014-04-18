-- Test to see effect of adiabatic elc without positivity

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.1

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
VL, VU = -6.0*vtPed, 6.0*vtPed
-- number of cells
NX, NV = 8, 32

-- L-B coefficient
lbAlpha = 0

-- parameters to control time-stepping
tStart = 0.0
tEnd = 350e-6
nFrames = 5

-- Determine number of global nodes per cell for use in creating CG
-- fields. Note that this looks a bit odd as this not the number of
-- *local* nodes but the number of nodes in each cell to give the
-- correct number of global nodes in fields.
if (polyOrder == 1) then
   numCgNodesPerCell = 1
   numCgNodesPerCell_1d = 1 -- for spatial basis
elseif (polyOrder == 2) then
   numCgNodesPerCell = 3
   numCgNodesPerCell_1d = 2 -- for spatial basis
end

-- phase space grid 
grid = Grid.RectCart2D {
   lower = {XL, VL},
   upper = {XU, VU},
   cells = {NX, NV},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.SerendipityElement {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell. One of 1, or 2. Corresponding
   -- number of nodes are 4 and 8.
   polyOrder = polyOrder,
}

-- Maxwellian
function maxwellian(vt, vdrift, x, v)
   return 1/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

-- Maxwellian (Right half only)
function maxwellianRight(vt, vdrift, x, v)
  if v >= 0 then
    return 1/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
  else
    return 0
  end
end

-- Maxwellian (Left half only)
function maxwellianLeft(vt, vdrift, x, v)
  if v <= 0 then
    return 1/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
  else
    return 0
  end
end

-- A generic function to run an updater.
--
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

-- distribution function
distf = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
distf:clear(0.0)

-- extra fields for performing RK update
distfNew = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
distf1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
diffDistf1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
dragDistf1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- Hamiltonian
hamil = DataStruct.Field2D {
   onGrid = grid,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = numCgNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- create field to store kinetic energy term in Hamiltonian
hamilKE = DataStruct.Field2D {
   onGrid = grid,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = numCgNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- DG versions of the hamiltonian
hamilKeDg = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- updater to initialize hamiltonian
initHamilKE = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = true,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}
initHamilKE:setOut( {hamilKE} )
-- initialize potential
initHamilKE:advance(0.0) -- time is irrelevant

-- create updater to copy a continuous field to a discontinuous field
copyCToD2D = Updater.CopyContToDisCont2D {
   onGrid = grid,
   basis = basis,
}
runUpdater(copyCToD2D, 0.0, 0.0, {hamilKE}, {hamilKeDg})

-- updater to initialize distribution function
initDistf = Updater.ProjectOnNodalBasis2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
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

     if x > lSource/2 then
       return nHat*maxwellianRight(vTi, 0.0, x, y)
     elseif x < -lSource/2 then
       return nHat*maxwellianLeft(vTi, 0.0, x, y)
     else
       -- Must be between the source boundaries, use a linear combo.
       return ((lSource/2 + x)*nHat*maxwellianRight(vTi, 0.0, x, y) + (lSource/2 - x)*nHat*maxwellianLeft(vTi, 0.0, x, y))/lSource
     end
   end
}
initDistf:setOut( {distf} )
-- initialize potential
initDistf:advance(0.0) -- time is irrelevant

-- particle source
particleSource = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
particleSource:clear(0.0)
-- updater to fill out particle source
particleSourceUpdater = Updater.ProjectOnNodalBasis2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
    if math.abs(x) < lSource/2 then
      if t < tELM then
        return Sn*math.cos(math.pi*x/lSource)*maxwellian(math.sqrt(tPed*eV/ionMass), 0.0, x, y)
      else
        return Sn/9*math.cos(math.pi*x/lSource)*maxwellian(math.sqrt(260*eV/ionMass), 0.0, x, y)
      end
    else
      return 0
    end
	 end
}
particleSourceUpdater:setOut( {particleSource} )
-- initialize
particleSourceUpdater:advance(0.0)
-- keeps track of whether to update particlesources again
postELM = false

-- updater for Poisson bracket
pbSlvr = Updater.PoissonBracket {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
}

-- updater for L-B drag term
lbDragSlvr = Updater.LenardBernsteinDragUpdater2D {
  onGrid = grid,
  -- basis functions to use
  basis = basis,
  -- cfl number to use
  cfl = cfl,
  -- Diffusion coefficient
  diffusionCoeff = lbAlpha,
  -- onlyIncrement is true
  onlyIncrement = true,
  -- use braginskii coeff instead of lbAlpha?
  useBraginskii = false,
  -- required if useBraginskii = true
  ionMass = ionMass,
}

-- updater for L-B drag term
lbDiffSlvr = Updater.LenardBernsteinDiffUpdater2D {
  onGrid = grid,
  -- basis functions to use
  basis = basis,
  -- cfl number to use
  cfl = cfl,
  -- Diffusion coefficient
  diffusionCoeff = lbAlpha,
  -- onlyIncrement is true
  onlyIncrement = true,
  -- use braginskii coeff instead of lbAlpha?
  useBraginskii = false,
  -- required if useBraginskii = true
  ionMass = ionMass,
}

-- spatial grid
grid_1d = Grid.RectCart1D {
   lower = {XL},
   upper = {XU},
   cells = {NX},
}

-- spatial FEM nodal basis
basis_1d = NodalFiniteElement1D.LagrangeTensor {
   -- grid on which elements should be constructured
   onGrid = grid_1d,
   -- polynomial order in each cell. One of 1, or 2. Corresponding
   -- number of nodes are 2 and 3.
   polyOrder = polyOrder,
}

-- drift velocity u(x)
driftU = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
driftU:clear(0.0)

-- thermal velocity squared
vThermSq = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
vThermSq:clear(0.0)

mom0Source = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

mom2Source = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- Ion Temperature Profile
tempProfile = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- field to store continous potential in 1D
phi1d = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numCgNodesPerCell_1d,
   -- ghost cells
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
}

-- field to store discontinous potential in 1D
phi1dDiscont = DataStruct.Field1D {
   onGrid = grid_1d,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- number density
numDensity = DataStruct.Field1D {
   onGrid = grid_1d,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- to compute number density
numDensityCalc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = grid,
   -- 2D phase-space basis functions
   basis2d = basis,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}

firstMoment = DataStruct.Field1D {
   onGrid = grid_1d,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- to compute first moment
firstMomentCalc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = grid,
   -- 2D phase-space basis functions
   basis2d = basis,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 1,
}

thirdMoment = DataStruct.Field1D {
   onGrid = grid_1d,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- to compute third moment
thirdMomentCalc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = grid,
   -- 2D phase-space basis functions
   basis2d = basis,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 3,
}

-- dynvector for total particle count
totalPtcl = DataStruct.DynVector { numComponents = 1, }

-- to compute total number of particles in domain
totalPtclCalc = Updater.IntegrateNodalField1D {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   -- are common nodes shared?
   shareCommonNodes = false, -- for DG fields common nodes not shared
}
-- set input field
totalPtclCalc:setIn( {numDensity} )
-- set output dynvector
totalPtclCalc:setOut( {totalPtcl} )

-- dynvector for number density in a cell
numDensInCell = DataStruct.DynVector { 
   numComponents = basis_1d:numNodes(),
}

-- to compute number density in a cell
numDensInCellCalc = Updater.RecordFieldInCell1D {
   -- grid for updater
   onGrid = grid_1d,
   -- index of cell to record
   cellIndex = {2},
}
-- set input field
numDensInCellCalc:setIn( {numDensity} )
-- set output dynvector
numDensInCellCalc:setOut( {numDensInCell} )

-- ptcl energy
ptclEnergy = DataStruct.Field1D {
   onGrid = grid_1d,
   location = "vertex",
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- updater to compute ptcl energy
ptclEnergyCalc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = grid,
   -- 2D phase-space basis functions
   basis2d = basis,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 2,
}

-- dynvector for total ptcl energy
totalPtclEnergy = DataStruct.DynVector { numComponents = 1, }

-- to compute total particle energy
totalPtclEnergyCalc = Updater.IntegrateNodalField1D {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   -- are common nodes shared?
   shareCommonNodes = false, -- for DG fields common nodes not shared
}
-- set input field
totalPtclEnergyCalc:setIn( {ptclEnergy} )
-- set output dynvector
totalPtclEnergyCalc:setOut( {totalPtclEnergy} )

-- Dynvector to store 0-3rd moments at left and right edges
momentsAtEdges = DataStruct.DynVector { numComponents = 8, }

momentsAtEdgesCalc = Updater.MomentsAtEdgesUpdater {
  onGrid = grid,
  basis = basis,
}

-- dynvector for heat flux at edge
heatFluxAtEdge = DataStruct.DynVector { numComponents = 3, }

-- to compute total particle energy
heatFluxAtEdgeCalc = Updater.HeatFluxAtEdgeUpdater {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   ionMass = ionMass,
   -- Perpendicular temperature of ions and electrons
   tPerp = tPed,
}

-- dynvector for total energy in system
totalEnergy = DataStruct.DynVector { numComponents = 3, }
-- updater to compute total energy
totalEnergyCalc = Updater.KineticTotalEnergyUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
  tPerp = tPed,
  ionMass = ionMass,
  electronMass = electronMass,
}

vFromMomentsCalc = Updater.VelocitiesFromMomentsUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
}

-- compute moments from distribution function
function calcMoments(curr, dt, distfIn)
   runUpdater(numDensityCalc, curr, dt, {distfIn}, {numDensity})
   runUpdater(firstMomentCalc, curr, dt, {distfIn}, {firstMoment})
   runUpdater(ptclEnergyCalc, curr, dt, {distfIn}, {ptclEnergy})
   runUpdater(thirdMomentCalc, curr, dt, {distfIn}, {thirdMoment})
end

-- updater to copy 1D field to 2D field
copyTo2D = Updater.NodalCopyFaceToInteriorUpdater {
   basis1d = basis_1d,
   basis2d = basis,
   onGrid = grid,
   -- 1D data lives along x-direction
   dir = 1,
   shareCommonNodes = true,
}

-- function to copy 1D field to 2D field
function copyPhi(curr, dt, phi1, phi2)
   copyTo2D:setCurrTime(curr)
   copyTo2D:setIn( {phi1} )
   copyTo2D:setOut( {phi2} )
   copyTo2D:advance(curr+dt)
   phi2:applyPeriodicBc(0)
end

-- updater to compute phi from number density
phiFromNumDensityCalc = Updater.BoltzmannPhiUpdater {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   -- required physical parameters
   electronMass = electronMass,
   ionMass = ionMass,
   elementaryCharge = eV,
}

-- updater to move phi to continuous field
phiToContCalc = Updater.ContFromDisCont1D {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
}

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

-- Outflow BCs
function getRepTbl(pOrder, val)
   if pOrder == 1 then
      return {val, val, val, val}
   elseif pOrder == 2 then
      return {val, val, val, val, val, val, val, val}
   end
end

function getCountTbl(pOrder, val)
   if pOrder == 1 then
      return {0, 1, 2, 3}
   elseif pOrder == 2 then
      return {0, 1, 2, 3, 4, 5, 6, 7}
   end
end

bcConst = BoundaryCondition.Const { 
   components = getCountTbl(polyOrder),
   values = getRepTbl(polyOrder, 0.0),
}

function makeBcObjIon()
   local bcLower = Updater.Bc2D {
      onGrid = grid,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "lower",
   }
   local bcUpper = Updater.Bc2D {
      onGrid = grid,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "upper",
   }
   return bcLower, bcUpper
end

-- make objects to apply BCs
bcLowerIon, bcUpperIon = makeBcObjIon()

-- apply boundary conditions to 2-D field
function applyBc(curr, dt, fldIon)
   for i,bc in ipairs({bcLowerIon, bcUpperIon}) do
      runUpdater(bc, curr, dt, {}, {fldIon})
   end

   fldIon:applyPeriodicBc(1)
   for i,fld in ipairs({fldElc, fldIon}) do
      fld:applyCopyBc(1, "lower")
      fld:applyCopyBc(1, "upper")
   end
end

-- update Poisson bracket operator
function poissonBracket(curr, dt, distfIn, hamilIn, distfOut)
   pbSlvr:setCurrTime(curr)
   pbSlvr:setIn( {distfIn, hamilIn} )
   pbSlvr:setOut( {distfOut} )
   return pbSlvr:advance(curr+dt)
end

-- update L-B drag operator
function solveDrag(curr, dt, distfIn, uIn, vThermSqIn, numDensityIn, distfOut)
   lbDragSlvr:setCurrTime(curr)
   lbDragSlvr:setIn( {distfIn, uIn, vThermSqIn, numDensityIn} )
   lbDragSlvr:setOut( {distfOut} )
   return lbDragSlvr:advance(curr+dt)
end

-- update L-B diffusion operator
function solveDiff(curr, dt, distfIn, vThermSqIn, numDensityIn, distfOut)
   lbDiffSlvr:setCurrTime(curr)
   lbDiffSlvr:setIn( {distfIn, vThermSqIn, numDensityIn} )
   lbDiffSlvr:setOut( {distfOut} )
   return lbDiffSlvr:advance(curr+dt)
end

-- compute hamiltonian
function calcHamiltonian(curr, dt, distIn, hamilOut)
   calcMoments(curr, dt, distIn)
   runUpdater(vFromMomentsCalc, curr, dt, {numDensity, firstMoment, ptclEnergy}, {driftU, vThermSq})

   phi1dDiscont:clear(0.0)
   runUpdater(phiFromNumDensityCalc, curr, dt, {numDensity, firstMoment, vThermSq}, {phi1dDiscont})

   phi1d:clear(0.0)
   runUpdater(phiToContCalc, curr, dt, {phi1dDiscont}, {phi1d})

   hamilOut:clear(0.0)
   -- Accumulate potential energy contribution to hamiltonian
   runUpdater(copyTo2D, curr, dt, {phi1d}, {hamilOut})
   hamilOut:scale(ionCharge/ionMass)
   -- Accumulate kinetic energy contribution to hamiltonian
   hamilOut:accumulate(1.0, hamilKE)
end

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   totalPtclCalc:setCurrTime(curr)
   totalPtclCalc:advance(curr+dt)

   numDensInCellCalc:setCurrTime(curr)
   numDensInCellCalc:advance(curr+dt)

   totalPtclEnergyCalc:setCurrTime(curr)
   totalPtclEnergyCalc:advance(curr+dt)

   -- compute moments at edges
   runUpdater(momentsAtEdgesCalc, curr, dt, {distf, hamilKeDg}, {momentsAtEdges})

   runUpdater(heatFluxAtEdgeCalc, curr, dt, {vThermSq, phi1dDiscont, momentsAtEdges}, {heatFluxAtEdge})

   fieldEnergyCalc:setCurrTime(curr)
   fieldEnergyCalc:advance(curr+dt)

   -- Compute temperature profile in eV
   tempProfile:copy(vThermSq)
   tempProfile:scale(ionMass/eV)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local pbStatus, pbDtSuggested
   local diffStatus, dragDtSuggested
   local dragStatus, diffDtSuggested

   if (postELM == false) then
     if tCurr + myDt > tELM then
       postELM = true
       particleSourceUpdater:advance(tCurr + myDt)
       -- Recalculate moments of particle source
       runUpdater(numDensityCalc, tCurr, myDt, {particleSource}, {mom0Source})
       runUpdater(ptclEnergyCalc, tCurr, myDt, {particleSource}, {mom2Source})
     end
   end

   -- RK stage 1
   pbStatus, pbDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {distf, hamil}, {distf1})
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf, vThermSq, numDensity, diffDistf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf, driftU, vThermSq, numDensity, dragDistf1)
   
   if (pbStatus == false or diffStatus == false or dragStatus == false) then
      return false, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
   end

   distf1:accumulate(myDt, diffDistf1)
   distf1:accumulate(myDt, dragDistf1)
   distf1:accumulate(myDt, particleSource)

   applyBc(tCurr, myDt, distf1)
   calcHamiltonian(tCurr, myDt, distf1, hamil)

   -- RK stage 2
   pbStatus, pbDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {distf1, hamil}, {distfNew})
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf1, vThermSq, numDensity, diffDistf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf1, driftU, vThermSq, numDensity, dragDistf1)

   if (pbStatus == false or diffStatus == false or dragStatus == false) then
     return false, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
   end

   distfNew:accumulate(myDt, diffDistf1)
   distfNew:accumulate(myDt, dragDistf1)
   distfNew:accumulate(myDt, particleSource)

   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)

   applyBc(tCurr, myDt, distf1)
   calcHamiltonian(tCurr, myDt, distf1, hamil)

   -- RK stage 3
   pbStatus, pbDtSuggested = runUpdater(pbSlvr, tCurr, myDt, {distf1, hamil}, {distfNew})
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf1, vThermSq, numDensity, diffDistf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf1, driftU, vThermSq, numDensity, dragDistf1)

   if (pbStatus == false or diffStatus == false or dragStatus == false) then
     return false, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
   end   
   
   distfNew:accumulate(myDt, diffDistf1)
   distfNew:accumulate(myDt, dragDistf1)
   distfNew:accumulate(myDt, particleSource)
   
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)

   applyBc(tCurr, myDt, distf1)
   distf:copy(distf1)
   calcHamiltonian(tCurr, myDt, distf, hamil)

   return true, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      distfDup:copy(distf)
      hamilDup:copy(hamil)

      if (tCurr+myDt > tEnd) then
	      myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	      distf:copy(distfDup)
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

-- write data to H5 files
function writeFields(frameNum, tCurr)
   --distf:write( string.format("distfIon_%d.h5", frame) )
   --numDensity:write( string.format("numDensity_%d.h5", frame) )
   --phi1dDiscont:write( string.format("phi_%d.h5", frame) )
   --hamil:write( string.format("hamil_%d.h5", frame) )
   --totalPtcl:write( string.format("totalPtcl_%d.h5", frame) )
   --totalPtclEnergy:write( string.format("totalPtclEnergy_%d.h5", frame) )
   --fieldEnergy:write( string.format("fieldEnergy_%d.h5", frame) )
   heatFluxAtEdge:write( string.format("heatFluxAtEdge_%d.h5", frameNum), tCurr)
   --tempProfile:write( string.format("tempProfile_%d.h5", frame) )
   --driftU:write( string.format("driftU_%d.h5", frame) )
   --firstMoment:write( string.format("mom1Ion_%d.h5", frame) , tCurr)
   --thirdMoment:write( string.format("mom3Ion_%d.h5", frame) , tCurr)
   --numDensInCell:write( string.format("numDensInCell_%d.h5", frame) )
end

-- Calculate first moment of particle source at t = 0
runUpdater(numDensityCalc, 0.0, 0.0, {particleSource}, {mom0Source})
-- Calculate second moment of particle source at t = 0
runUpdater(ptclEnergyCalc, 0.0, 0.0, {particleSource}, {mom2Source})

applyBc(0.0, 0.0, distf)
-- calculate initial Hamiltonian
calcHamiltonian(0.0, 0.0, distf, hamil)
-- compute initial diagnostics
calcDiagnostics(0.0, 0.0)
-- write out initial conditions
writeFields(0, 0.0)
-- make a duplicate in case we need it
distfDup = distf:duplicate()
hamilDup = hamil:duplicate()

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
