-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.1

-- physical constants
-- eletron mass (kg)
electronMass = 9.10938215e-31
-- electron volt to joules
eV = 1.602e-19
-- Deuterium ion mass (kg)
ionMass = 2.014*1.66e-27
-- signed ion charge
ionCharge = eV
-- permittivity of free space (F/m)
epsilon0 = 8.85418782e-12;

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
nMidplane = nPed/5;
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
NX, NV = 32, 32

-- initial elm ion source temperature
ionTemp = tPed

-- L-B coefficient
lbAlpha = 0

-- parameters to control time-stepping
tStart = 0.0
tEnd = 1000e-6
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
   numComponents = 1*numCgNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- create field to store kinetic energy term in Hamiltonian
hamilKE = DataStruct.Field2D {
   onGrid = grid,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numCgNodesPerCell,
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
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

-- updater to initialize distribution function
initDistf = Updater.EvalOnNodes2D {
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
     
     -- Figure out drift velocity
     local driftVelocity
     -- Omotani and Dudson pre-elm S_n
     local preSn = 7.4e22
     
     if math.abs(x) < lSource/2 then
      driftVelocity = lSource*preSn*math.sin(math.pi*x/lSource)/(math.pi*nMidplane)
     else
      driftVelocity = (x/math.abs(x))*lSource*preSn/(math.pi*nMidplane)
     end

		 return backgroundDens*maxwellian(math.sqrt(backgroundTemp*eV/ionMass), driftVelocity, x, y)
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
particleSourceUpdater = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
    if math.abs(x) < lSource/2 then
      return Sn*math.cos(math.pi*x/lSource)*maxwellian(math.sqrt(ionTemp*eV/ionMass), 0.0, x, y)
    else
      return 0
    end
	 end
}
particleSourceUpdater:setOut( {particleSource} )
-- initialize
particleSourceUpdater:advance(0.0) -- time is irrelevant

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
  useBraginskii = true,
  -- required if useBraginskii = true
  ionMass = ionMass,
  elementaryCharge = eV,
  epsilon0 = epsilon0,
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
  useBraginskii = true,
  -- required if useBraginskii = true
  ionMass = ionMass,
  elementaryCharge = eV,
  epsilon0 = epsilon0,
}

-- spatial grid
grid_1d = Grid.RectCart1D {
   lower = {XL},
   upper = {XU},
   cells = {NX},
}

-- spatial FEM nodal basis
basis_1d = NodalFiniteElement1D.Lobatto {
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

-- dynvector for heat flux at edge
heatFluxAtEdge = DataStruct.DynVector { numComponents = 3, }

-- to compute total particle energy
heatFluxAtEdgeCalc = Updater.HeatFluxAtEdgeUpdater {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   -- are common nodes shared?
   shareCommonNodes = false,
   ionMass = ionMass,
   elementaryCharge = eV,
}

-- set input field
heatFluxAtEdgeCalc:setIn( {firstMoment, thirdMoment, vThermSq, phi1dDiscont} )
-- set output dynvector
heatFluxAtEdgeCalc:setOut( {heatFluxAtEdge} )

vFromMomentsCalc = Updater.VelocitiesFromMomentsUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
}
vFromMomentsCalc:setIn( {numDensity, firstMoment, ptclEnergy} )
vFromMomentsCalc:setOut( {driftU, vThermSq} )

-- calculate number density
function calcNumDensity(curr, dt, distfIn, numDensOut)
   numDensityCalc:setCurrTime(curr)
   numDensityCalc:setIn( {distfIn} )
   numDensityCalc:setOut( {numDensOut} )
   numDensityCalc:advance(curr+dt)
end

-- calculate first velocity moment (moment 1)
function calcFirstMoment(curr, dt, distfIn, firstMomentOut)
   firstMomentCalc:setCurrTime(curr)
   firstMomentCalc:setIn( {distfIn} )
   firstMomentCalc:setOut( {firstMomentOut} )
   firstMomentCalc:advance(curr+dt)
end

-- calculate ptcl energy
function calcPtclEnergy(curr, dt, distfIn, energyOut)
   ptclEnergyCalc:setCurrTime(curr)
   ptclEnergyCalc:setIn( {distfIn} )
   ptclEnergyCalc:setOut( {energyOut} )
   ptclEnergyCalc:advance(curr+dt)
end

-- calculate third velocity moment (moment 3)
function calcThirdMoment(curr, dt, distfIn, thirdMomentOut)
   thirdMomentCalc:setCurrTime(curr)
   thirdMomentCalc:setIn( {distfIn} )
   thirdMomentCalc:setOut( {thirdMomentOut} )
   thirdMomentCalc:advance(curr+dt)
end

-- compute moments from distribution function
function calcMoments(curr, dt, distfIn)
   calcNumDensity(curr, dt, distfIn, numDensity)
   calcFirstMoment(curr, dt, distfIn, firstMoment)
   calcPtclEnergy(curr, dt, distfIn, ptclEnergy)
   calcThirdMoment(curr, dt, distfIn, thirdMoment)
end

-- updater to copy 1D field to 2D field
copyTo2D = Updater.Copy1DTo2DNodalField {
   -- grid for updater
   onGrid = grid,
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

-- compute phi from number density and thermal velocity
function calcPhiFromNumDensity(curr, dt, numDensIn, vtSqIn, mom1In, phiOut)
   phiOut:clear(0.0)
   phiFromNumDensityCalc:setCurrTime(curr)
   phiFromNumDensityCalc:setIn( {numDensIn, mom1In, vtSqIn} )
   phiFromNumDensityCalc:setOut( {phiOut} )
   return phiFromNumDensityCalc:advance(curr+dt)
end

-- compute continuous phi from discontinuous phi
function calcContinuousPhi(curr, dt, phiIn, phiOut)
   phiOut:clear(0.0)
   phiToContCalc:setCurrTime(curr)
   phiToContCalc:setIn( {phiIn} )
   phiToContCalc:setOut( {phiOut} )
   return phiToContCalc:advance(curr+dt)
end

function calcVelocitiesFromMoments(curr, dt)
  vFromMomentsCalc:setCurrTime(curr)
  vFromMomentsCalc:advance(curr+dt)
end

-- compute hamiltonian
function calcHamiltonian(curr, dt, distIn, hamilOut)
   calcMoments(curr, dt, distIn)
   calcVelocitiesFromMoments(curr, dt)
   calcPhiFromNumDensity(curr, dt, numDensity, vThermSq, firstMoment, phi1dDiscont)
   calcContinuousPhi(curr, dt, phi1dDiscont, phi1d)

   hamilOut:clear(0.0)
   copyPhi(curr, dt, phi1d, hamilOut) -- potential energy contribution
   hamilOut:scale(ionCharge/ionMass)
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

   heatFluxAtEdgeCalc:setCurrTime(curr)
   heatFluxAtEdgeCalc:advance(curr+dt)

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

   -- RK stage 1
   pbStatus, pbDtSuggested = poissonBracket(tCurr, myDt, distf, hamil, distf1)
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf, vThermSq, numDensity, diffDistf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf, driftU, vThermSq, numDensity, dragDistf1)

   distf1:accumulate(myDt, diffDistf1)
   distf1:accumulate(myDt, dragDistf1)
   if tCurr < tELM then
    distf1:accumulate(myDt, particleSource)
   end
   
   if (pbStatus == false or diffStatus == false or dragStatus == false) then
      return false, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
   end

   applyBc(tCurr, myDt, distf1)
   calcHamiltonian(tCurr, myDt, distf1, hamil)

   -- RK stage 2
   pbStatus, pbDtSuggested = poissonBracket(tCurr, myDt, distf1, hamil, distfNew)
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf1, vThermSq, numDensity, diffDistf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf1, driftU, vThermSq, numDensity, dragDistf1)

   distfNew:accumulate(myDt, diffDistf1)
   distfNew:accumulate(myDt, dragDistf1)
   if tCurr < tELM then
    distfNew:accumulate(myDt, particleSource)
   end

   if (pbStatus == false or diffStatus == false or dragStatus == false) then
     return false, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
   end

   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)
   applyBc(tCurr, myDt, distf1)
   calcHamiltonian(tCurr, myDt, distf1, hamil)

   -- RK stage 3
   pbStatus, pbDtSuggested = poissonBracket(tCurr, myDt, distf1, hamil, distfNew)
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf1, vThermSq, numDensity, diffDistf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf1, driftU, vThermSq, numDensity, dragDistf1)
   
   distfNew:accumulate(myDt, diffDistf1)
   distfNew:accumulate(myDt, dragDistf1)
   if tCurr < tELM then
    distfNew:accumulate(myDt, particleSource)
   end

   if (pbStatus == false or diffStatus == false or dragStatus == false) then
     return false, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
   end
   
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

      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 distf:copy(distfDup)
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
function writeFields(frame)
   --distf:write( string.format("distf_%d.h5", frame) )
   --numDensity:write( string.format("numDensity_%d.h5", frame) )
   phi1dDiscont:write( string.format("phi1d_%d.h5", frame) )
   --hamil:write( string.format("hamil_%d.h5", frame) )
   --totalPtcl:write( string.format("totalPtcl_%d.h5", frame) )
   --totalPtclEnergy:write( string.format("totalPtclEnergy_%d.h5", frame) )
   --fieldEnergy:write( string.format("fieldEnergy_%d.h5", frame) )
   heatFluxAtEdge:write( string.format("heatFluxAtEdge_%d.h5", frame) )
   --tempProfile:write( string.format("tempProfile_%d.h5", frame) )
   --vThermSq:write( string.format("vThermSq_%d.h5", frame) )
   driftU:write( string.format("driftU_%d.h5", frame) )
   firstMoment:write( string.format("mom1_%d.h5", frame) )
   --numDensInCell:write( string.format("numDensInCell_%d.h5", frame) )
end

--hamilKE:write("hamilKE.h5")
--numDensityElc:write("numDensityElc.h5")
applyBc(0.0, 0.0, distf)
-- calculate initial Hamiltonian
calcHamiltonian(0.0, 0.0, distf, hamil)
-- compute initial diagnostics
calcDiagnostics(0.0, 0.0)
-- write out initial conditions
writeFields(0)
-- make a duplicate in case we need it
distfDup = distf:duplicate()

-- parameters to control time-stepping
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
tFrame = (tEnd-tStart)/nFrames
tCurr = tStart

for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
