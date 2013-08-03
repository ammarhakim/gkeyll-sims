-- Input file for testing 1x1v L-B collision operator
-- Mainly used to look at conservation properties

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.1

-- L-B coefficient
lbAlpha = 1.0

-- wave-number
knumber = 0.5

-- domain extents
XL, XU = -0.5, 0.5
VL, VU = -6.0, 6.0
-- number of cells
NX, NV = 1, 32

-- parameters to control time-stepping
tStart = 0.0
tEnd = 80.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 1

-- phase space grid 
grid = Grid.RectCart2D {
   lower = {XL, VL},
   upper = {XU, VU},
   cells = {NX, NV},
}

-- spatial grid
grid_1d = Grid.RectCart1D {
   lower = {XL},
   upper = {XU},
   cells = {NX},
}

-- create FEM nodal basis
basis = NodalFiniteElement2D.SerendipityElement {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell
   polyOrder = polyOrder,
}

-- spatial FEM nodal basis
basis_1d = NodalFiniteElement1D.Lobatto {
   -- grid on which elements should be constructured
   onGrid = grid_1d,
   -- polynomial order in each cell. One of 1, or 2. Corresponding
   -- number of nodes are 2 and 3.
   polyOrder = polyOrder,
}

-- distribution function
distf = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
distf:clear(0.0)

-- drift velocity u(x)
driftU = DataStruct.Field1D {
   onGrid = grid_1d,
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

dragDistf1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- duplicate for use in time-stepping
distfNewDup = distfNew:duplicate()

-- updater to apply initial conditions
initField = Updater.EvalOnNodes2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
     if (math.abs(y) <= 1) then
       return (math.pi/4)*math.cos(math.pi*y/2)
     else
       return 0
     end
	 end,
}
initField:setOut( {distf} )
-- initialize
initField:advance(0.0) -- time is irrelevant

-- updater for L-B drag term
lbDragSlvr = Updater.LenardBernsteinDragUpdater2D {
  onGrid = grid,
  -- basis functions to use
  basis = basis,
  -- cfl number to use
  cfl = cfl,
  -- Diffusion coefficient
  diffusionCoeff = lbAlpha,
  onlyIncrement = true,
}

-- updater for L-B diffusion term
lbDiffSlvr = Updater.LenardBernsteinDiffUpdater2D {
  onGrid = grid,
  -- basis functions to use
  basis = basis,
  -- cfl number to use
  cfl = cfl,
  -- Diffusion coefficient
  diffusionCoeff = lbAlpha,
}

-- apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyPeriodicBc(1)
end

applyBc(distf)
distfNew:copy(distf)

-- recovery poly eval on bot edges
--recoveryBotEvals = DataStruct.Field1D {
--   onGrid = grid_1d,
--   numComponents = basis_1d:numNodes(),
--   ghost = {1, 1},
--}
-- recovery poly eval on top edges
--recoveryTopEvals = DataStruct.Field1D {
--   onGrid = grid_1d,
--   numComponents = basis_1d:numNodes(),
--   ghost = {1, 1},
--}
-- updater to compute recovery poly at edges
--recoveryEvalsCalc = Updater.RecoveryPolynomialUpdater {
--   -- 2D phase-space grid 
--   onGrid = grid,
--   -- 2D phase-space basis functions
--   basis2d = basis,
--   -- 1D spatial basis functions
--   basis1d = basis_1d,
--}

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

-- dynvector for total ptcl energy
totalPtclMomentum = DataStruct.DynVector { numComponents = 1, }

-- to compute total particle energy
totalPtclMomentumCalc = Updater.IntegrateNodalField1D {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   -- are common nodes shared?
   shareCommonNodes = false, -- for DG fields common nodes not shared
}
-- set input field
totalPtclMomentumCalc:setIn( {firstMoment} )
-- set output dynvector
totalPtclMomentumCalc:setOut( {totalPtclMomentum} )

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

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   totalPtclCalc:setCurrTime(curr)
   totalPtclCalc:advance(curr+dt)

   totalPtclMomentumCalc:setCurrTime(curr)
   totalPtclMomentumCalc:advance(curr+dt)

   totalPtclEnergyCalc:setCurrTime(curr)
   totalPtclEnergyCalc:advance(curr+dt)
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

-- calculate number density
function calcNumDensity(curr, dt, distfIn, numDensOut)
   numDensityCalc:setCurrTime(curr)
   numDensityCalc:setIn( {distfIn} )
   numDensityCalc:setOut( {numDensOut} )
   numDensityCalc:advance(curr+dt)
   
   numDensOut:applyPeriodicBc(0)
end

-- calculate first velocity moment (moment 1)
function calcFirstMoment(curr, dt, distfIn, firstMomentOut)
   firstMomentCalc:setCurrTime(curr)
   firstMomentCalc:setIn( {distfIn} )
   firstMomentCalc:setOut( {firstMomentOut} )
   firstMomentCalc:advance(curr+dt)

   firstMomentOut:applyPeriodicBc(0)
end

-- calculate ptcl energy (moment 2)
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
   -- compute recovery poly on edges
   -- recoveryEvalsCalc:setIn( {distfIn} )
   -- recoveryEvalsCalc:setOut( {recoveryBotEvals, recoveryTopEvals} )
   -- recoveryEvalsCalc:setCurrTime(curr)
   -- recoveryEvalsCalc:advance(curr+dt)

   calcVelocitiesFromMoments(curr, dt)
end

function calcVelocitiesFromMoments(curr, dt)
  vFromMomentsCalc:setCurrTime(curr)
  vFromMomentsCalc:advance(curr+dt)

  driftU:applyPeriodicBc(0)
  vThermSq:applyPeriodicBc(0)
end

function rk3(tCurr, myDt)
   local diffStatus, diffDtSuggested
   local dragStatus, dragDtSuggested

   -- RK stage 1
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf, vThermSq, distf1)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf, driftU, dragDistf1)
   distf1:accumulate(myDt, dragDistf1)
   
   if (diffStatus == false or dragStatus == false) then
      return false, math.min(diffDtSuggested, dragDtSuggested)
   end

   applyBc(distf1)
   calcMoments(tCurr, myDt, distf1)

   -- RK stage 2
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf1, vThermSq, distfNew)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf1, driftU, dragDistf1)
   distfNew:accumulate(myDt, dragDistf1)

   if (diffStatus == false or dragStatus == false) then
     return false, math.min(diffDtSuggested, dragDtSuggested)
   end

   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)
   applyBc(distf1)
   calcMoments(tCurr, myDt, distf1)

   -- RK stage 3
   diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distf1, vThermSq, distfNew)
   dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distf1, driftU, dragDistf1)
   distfNew:accumulate(myDt, dragDistf1)

   if (diffStatus == false or dragStatus == false) then
     return false, math.min(diffDtSuggested, dragDtSuggested)
   end
   
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)
   applyBc(distf1)
   distf:copy(distf1)
   calcMoments(tCurr, myDt, distf)

   return true, math.min(diffDtSuggested, dragDtSuggested)
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      distfNewDup:copy(distfNew)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	      myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	      -- time-step too large
	      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	      distfNew:copy(distfNewDup)
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

-- write data to H5 file
function writeFields(frame)
   distf:write( string.format("distf_%d.h5", frame) )
   --numDensity:write( string.format("numDensity_%d.h5", frame) )
   totalPtcl:write( string.format("totalPtcl_%d.h5", frame) )
   totalPtclEnergy:write( string.format("totalPtclEnergy_%d.h5", frame) )
   totalPtclMomentum:write( string.format("totalPtclMomentum_%d.h5", frame) )
   --firstMoment:write( string.format("mom1_%d.h5", frame) )
   --ptclEnergy:write( string.format("mom2_%d.h5", frame) )
   --vThermSq:write( string.format("vThermSq_%d.h5", frame) )
   --driftU:write( string.format("driftU_%d.h5", frame) )
end

-- compute initial driftU and vThermSq
calcMoments(0.0, 0.0, distf)
calcDiagnostics(0.0, 0.0)
-- write initial conditions
writeFields(0)

tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
