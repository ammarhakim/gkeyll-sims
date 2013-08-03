-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.05

-- wave-number
knumber = 0.5

-- domain extents
XL, XU = -Lucee.Pi/knumber, Lucee.Pi/knumber
VL, VU = -6.0, 6.0
-- number of cells
NX, NV = 32, 32

-- temperature ratio (T_i/T_e)
Tratio = 0.25

-- ion temperature
ionTemp = 1.0
-- electron temperature
elcTemp = ionTemp/Tratio
-- signed ion charge
ionCharge = 1.0
-- electron mass
ionMass = 1.0
-- permittivity of free space
epsilon0 = 1.0
-- L-B coefficient
lbAlpha = 0

-- parameters to control time-stepping
tStart = 0.0
tEnd = 100.0
nFrames = 1

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

-- distribution function for electrons
distfElc = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
distfElc:clear(0.0)

-- Maxwellian
function maxwellian(vt, vdrift, x, v)
   return 1/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

-- updater to initialize distribution function
initDistfElc = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		 return maxwellian(elcTemp, 0.0, x, y)
	      end
}
initDistfElc:setOut( {distfElc} )
-- initialize ion distribution function
initDistfElc:advance(0.0) -- time is irrelevant

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
hamilKE:applyPeriodicBc(0)

-- updater to initialize distribution function
initDistf = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		 local alpha = 0.01 -- perturbation
		 local k = knumber
		 return (1+alpha*math.cos(k*x))*maxwellian(ionTemp, 0.0, x, y)
	      end
}
initDistf:setOut( {distf} )
-- initialize potential
initDistf:advance(0.0) -- time is irrelevant

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

-- number density
numDensityElc = DataStruct.Field1D {
   onGrid = grid_1d,
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

   numDensOut:applyPeriodicBc(0)
end

-- compute number density for electrons
calcNumDensity(0.0, 0.0, distfElc, numDensityElc)

-- calculate first velocity moment (moment 1)
function calcFirstMoment(curr, dt, distfIn, firstMomentOut)
   firstMomentCalc:setCurrTime(curr)
   firstMomentCalc:setIn( {distfIn} )
   firstMomentCalc:setOut( {firstMomentOut} )
   firstMomentCalc:advance(curr+dt)

   firstMomentOut:applyPeriodicBc(0)
end

-- calculate ptcl energy
function calcPtclEnergy(curr, dt, distfIn, energyOut)
   ptclEnergyCalc:setCurrTime(curr)
   ptclEnergyCalc:setIn( {distfIn} )
   ptclEnergyCalc:setOut( {energyOut} )
   ptclEnergyCalc:advance(curr+dt)

   energyOut:applyPeriodicBc(0)
end

-- compute moments from distribution function
function calcMoments(curr, dt, distfIn)
   calcNumDensity(curr, dt, distfIn, numDensity)
   calcFirstMoment(curr, dt, distfIn, firstMoment)
   calcPtclEnergy(curr, dt, distfIn, ptclEnergy)
end

-- function to apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")
end

-- apply BCs to initial conditions
applyBc(distf)

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
phiFromNumDensityCalc = Updater.ContFromDisCont1D {
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

-- update Poisson bracket operator
function poissonBracket(curr, dt, distfIn, hamilIn, distfOut)
   pbSlvr:setCurrTime(curr)
   pbSlvr:setIn( {distfIn, hamilIn} )
   pbSlvr:setOut( {distfOut} )
   return pbSlvr:advance(curr+dt)
end

-- update L-B drag operator
function solveDrag(curr, dt, distfIn, uIn, distfOut)
   lbDragSlvr:setCurrTime(curr)
   lbDragSlvr:setIn( {distfIn, uIn} )
   lbDragSlvr:setOut( {distfOut} )
   return lbDragSlvr:advance(curr+dt)
end

-- update L-B diffusion operator
function solveDiff(curr, dt, distfIn, vThermSqIn, distfOut)
   lbDiffSlvr:setCurrTime(curr)
   lbDiffSlvr:setIn( {distfIn, vThermSqIn} )
   lbDiffSlvr:setOut( {distfOut} )
   return lbDiffSlvr:advance(curr+dt)
end

-- compute phi from number density
function calcPhiFromNumDensity(curr, dt, numDensIn, phiOut)
   phiOut:clear(0.0)
   phiFromNumDensityCalc:setCurrTime(curr)
   numDensIn:accumulate(-1.0, numDensityElc)
   numDensIn:scale(elcTemp)
   numDensIn:applyPeriodicBc(0)
   phiFromNumDensityCalc:setIn( {numDensIn} )
   phiFromNumDensityCalc:setOut( {phiOut} )
   return phiFromNumDensityCalc:advance(curr+dt)
end

function calcVelocitiesFromMoments(curr, dt)
  vFromMomentsCalc:setCurrTime(curr)
  vFromMomentsCalc:advance(curr+dt)

  driftU:applyPeriodicBc(0)
  vThermSq:applyPeriodicBc(0)
end

-- compute hamiltonian
function calcHamiltonian(curr, dt, distIn, hamilOut)
   calcMoments(curr, dt, distIn)
   calcVelocitiesFromMoments(curr, dt)
   calcPhiFromNumDensity(curr, dt, numDensity, phi1d)
   hamilOut:clear(0.0)
   copyPhi(curr, dt, phi1d, hamilOut) -- potential energy contribution
   hamilOut:scale(ionCharge/ionMass)
   hamilOut:accumulate(1.0, hamilKE)
   hamilOut:applyPeriodicBc(0)
end

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   totalPtclCalc:setCurrTime(curr)
   totalPtclCalc:advance(curr+dt)

   numDensInCellCalc:setCurrTime(curr)
   numDensInCellCalc:advance(curr+dt)

   totalPtclEnergyCalc:setCurrTime(curr)
   totalPtclEnergyCalc:advance(curr+dt)

   fieldEnergyCalc:setCurrTime(curr)
   fieldEnergyCalc:advance(curr+dt)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local pbStatus, pbDtSuggested
   local diffStatus, dragDtSuggested
   local dragStatus, diffDtSuggested

   -- RK stage 1
   pbStatus, pbDtSuggested = poissonBracket(tCurr, myDt, distf, hamil, distf1)
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf, vThermSq, diffDistf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf, driftU, dragDistf1)

   distf1:accumulate(myDt, diffDistf1)
   distf1:accumulate(myDt, dragDistf1)
   
   if (pbStatus == false or diffStatus == false or dragStatus == false) then
      return false, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
   end

   applyBc(distf1)
   calcHamiltonian(tCurr, myDt, distf1, hamil)

   -- RK stage 2
   pbStatus, pbDtSuggested = poissonBracket(tCurr, myDt, distf1, hamil, distfNew)
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf1, vThermSq, diffDistf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf1, driftU, dragDistf1)

   distfNew:accumulate(myDt, diffDistf1)
   distfNew:accumulate(myDt, dragDistf1)

   if (pbStatus == false or diffStatus == false or dragStatus == false) then
     return false, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
   end

   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)
   applyBc(distf1)
   calcHamiltonian(tCurr, myDt, distf1, hamil)

   -- RK stage 3
   pbStatus, pbDtSuggested = poissonBracket(tCurr, myDt, distf1, hamil, distfNew)
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf1, vThermSq, diffDistf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf1, driftU, dragDistf1)
   
   distfNew:accumulate(myDt, diffDistf1)
   distfNew:accumulate(myDt, dragDistf1)

   if (pbStatus == false or diffStatus == false or dragStatus == false) then
     return false, math.min(pbDtSuggested, diffDtSuggested, dragDtSuggested)
   end
   
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)
   applyBc(distf1)
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
   distf:write( string.format("distf_%d.h5", frame) )
   numDensity:write( string.format("numDensity_%d.h5", frame) )
   phi1d:write( string.format("phi1d_%d.h5", frame) )
   hamil:write( string.format("hamil_%d.h5", frame) )
   totalPtcl:write( string.format("totalPtcl_%d.h5", frame) )
   totalPtclEnergy:write( string.format("totalPtclEnergy_%d.h5", frame) )
   fieldEnergy:write( string.format("fieldEnergy_%d.h5", frame) )
   --numDensInCell:write( string.format("numDensInCell_%d.h5", frame) )
end

--hamilKE:write("hamilKE.h5")
-- initial dist f and num density for elec
--distfElc:write("distfElc.h5")
--numDensityElc:write("numDensityElc.h5")
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
