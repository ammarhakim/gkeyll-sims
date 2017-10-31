-- Input file for testing 2-D recovery DG with precomputed matrices

-- polynomial order
polyOrder = 1
-- DT factor
dtFact = 200
-- cfl number to use
cfl = 0.2*dtFact*10 -- don't ever fail
-- L-B coefficient
lbAlpha = 1.0
-- wave number
knumber = 0.5

-- domain extents
XL, XU = -Lucee.Pi/knumber, Lucee.Pi/knumber
VL, VU = -6.0, 6.0
-- number of cells
NX, NV = 4, 32

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
   -- polynomial order in each cell
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
-- Initialize vThermSq
vThermSqUpdater = Updater.ProjectOnNodalBasis1D {
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		 return 1.0
	      end
}
vThermSqUpdater:setOut( {vThermSq} )
vThermSqUpdater:advance(0.0)
-- write initial vt^2 condition
distf:write("vt_0.h5")

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
diffDistf1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
deltaDist = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
maxDist = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- to store second-order derivative of field
distLBO = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
distLBO0 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
distLBO1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
distLBO2 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- duplicate for use in time-stepping
distfNewDup = distfNew:duplicate()

function maxwellian(vt, vdrift, x, v)
   return 1/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
end

function step(x, y)
   if (math.abs(y) <= 2) then
      return 0.25
   else
      return 0
   end
end

function gaussian(x, v)
   return math.exp(-v^2)
end

-- updater to apply initial conditions
initField = Updater.ProjectOnNodalBasis2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 return step(x,y)
	      end,
}
initField:setOut( {distf} )
-- initialize
initField:advance(0.0) -- time is irrelevant

-- this set of functions determines factors which feed into RK scheme
-- (see Meyer, C. D., Balsara, D. S., & Aslam, T. D. (2014). Journal
-- of Computational Physics, 257(PA),
-- 594–626. doi:10.1016/j.jcp.2013.08.021)
function b(j)
   if (j<2) then 
      return 1.0/3.0
   else 
      return (j^2+j-2)/(2*j*(j+1))
   end
end
function a(j) return 1-b(j) end
function w1(s) return 4/(s^2+s-2) end
function mubar(s,j) 
   if (j<2) then 
      return 4/(3*(s^2+s-2)) 
   else 
      return 4*(2*j-1)/(j*(s^2+s-2))*b(j)/b(j-1)
   end
end
function mu(j) return (2*j-1)/j*b(j)/b(j-1) end
function nu(j) return -(j-1)/j*b(j)/b(j-2) end
function gbar(s,j) return -a(j-1)*mubar(s,j) end
function roundToOdd(n) return n % 2 == 0 and n+1 or n end
function calcNumStages(dhdp)
   return roundToOdd(math.ceil (math.sqrt(4*dhdp+9/4) - 1/2))+2
end

-- updater to apply initial conditions
initMax = Updater.ProjectOnNodalBasis2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 return maxwellian(1.0, 0.0, x, y)
	      end,
}
initMax:setOut( {maxDist} )
-- initialize
initMax:advance(0.0) -- time is irrelevant

-- updater for L-B drag term
lbDragSlvr = Updater.LenardBernsteinDragUpdater2D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- cfl number to use
   cfl = cfl,
   -- Diffusion coefficient
   diffusionCoeff = lbAlpha,
   -- compute increments only?
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
   -- compute increments only?
   onlyIncrement = true,
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

-- calculate number density
function calcNumDensity(curr, dt, distfIn, numDensOut)
   numDensityCalc:setCurrTime(curr)
   numDensityCalc:setIn( {distfIn} )
   numDensityCalc:setOut( {numDensOut} )
   numDensityCalc:advance(curr+dt)

   numDensOut:applyPeriodicBc(0)
end

function calcDiagnostics(curr, dt)
   totalPtclCalc:setCurrTime(curr)
   totalPtclCalc:advance(curr+dt)
end

bcConst = BoundaryCondition.Const { components = {0, 1, 2, 4},  values = {0, 0, 0, 0} }
bcFunc = BoundaryCondition.Function { 
   components = {0, 1, 2, 4},  
   bc = function (x,y,z,t)
	   local val = math.exp(-y^2/2)
	   return val,val,val,val
	end
}
bcLower = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcConst},
   -- direction to apply
   dir = 1,
   -- edge to apply on
   edge = "lower",
}
bcUpper = Updater.Bc2D {
   onGrid = grid,
   -- boundary conditions to apply
   boundaryConditions = {bcConst},
   -- direction to apply
   dir = 1,
   -- edge to apply on
   edge = "upper",
}

-- apply boundary conditions
function applyBc(fld)
   fld:applyCopyBc(0, "lower")
   fld:applyCopyBc(0, "upper")
end

applyBc(distf)
distfNew:copy(distf)

-- write initial conditions
distf:write("q_0.h5")

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

function calcLBO(tCurr, myDt, distfIn, distLBO)
   local diffStatus, diffDtSuggested = solveDiff(tCurr, myDt, distfIn, vThermSq, diffDistf1)
   local dragStatus, dragDtSuggested = solveDrag(tCurr, myDt, distfIn, driftU, dragDistf1)
   distLBO:combine(1.0, diffDistf1, 1.0, dragDistf1)
   return diffStatus and dragStatus, math.min(diffDtSuggested, dragDtSuggested)
end

print("NumStages", calcNumStages(dtFact))

-- function to take a time-step using RK-L scheme of Meyer et. al.
function rkLegendre2(tCurr, myDt)
   local dtRatio = dtFact
   local numStages = calcNumStages(dtRatio)

   -- we need this in each stage
   calcLBO(tCurr, myDt, distf, distLBO0)
   -- stage 1
   distLBO2:copy(distf)
   distLBO1:combine(1.0, distf, mubar(numStages,1)*myDt, distLBO0)
   applyBc(distLBO1)

   -- rest of stages
   for j = 2, numStages do
      calcLBO(tCurr, myDt, distLBO1, distf1)
      distLBO:combine(mu(j), distLBO1, nu(j), distLBO2, 1-mu(j)-nu(j), distf,
		 mubar(numStages,j)*myDt, distf1, gbar(numStages,j)*myDt, distLBO0)
      applyBc(distLBO)
      -- reset fields for next stage
      distLBO2:copy(distLBO1)
      distLBO1:copy(distLBO)
   end
   distf:copy(distLBO)

   return true, dtFact*0.0133333
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local status, dtSuggested

   -- RK stage 1
   deltaDist:copy(distf)
   status, dtSuggested = calcLBO(tCurr, myDt, deltaDist, distf1)
   distf1:scale(myDt)
   distf1:accumulate(1.0, distf)

   if (status == false) then
      return false, dtSuggested
   end
   applyBc(distf1)

   -- RK stage 2
   deltaDist:copy(distf1)
   status, dtSuggested = calcLBO(tCurr, myDt, deltaDist, distfNew)
   distfNew:scale(myDt)
   distfNew:combine(1.0, distf1)

   if (status == false) then
      return false, dtSuggested
   end

   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)
   applyBc(distf1)

   -- RK stage 3
   deltaDist:copy(distf1)
   status, dtSuggested = calcLBO(tCurr, myDt, deltaDist, distfNew)
   distfNew:scale(myDt)
   distfNew:combine(1.0, distf1)

   if (diffStatus == false or dragStatus == false) then
      return false, math.min(diffDtSuggested, dragDtSuggested)
   end
   
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)
   applyBc(distf1)
   distf:copy(distf1)

   return true, dtSuggested
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
      --status, dtSuggested = rk3(tCurr, myDt)
      status, dtSuggested = rkLegendre2(tCurr, myDt)
      if (status == false) then
	 -- time-step too large
	 print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 distfNew:copy(distfNewDup)
	 myDt = dtSuggested
      else
	 -- compute diagnostics
	 calcNumDensity(tCurr, myDt, distfNew, numDensity)
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

totalPtcl:write( string.format("totalPtcl_%d.h5", 0) )
-- write data to H5 file
function writeFields(frame)
   distf:write( string.format("distf_%d.h5", frame) )
   numDensity:write( string.format("numDensity_%d.h5", frame) )
   totalPtcl:write( string.format("totalPtcl_%d.h5", frame) )
end

-- parameters to control time-stepping
tStart = 0.0
tEnd = 50.0
dtSuggested = 0.0133333
nFrames = 1
tFrame = (tEnd-tStart)/nFrames -- time between frames

tCurr = tStart
for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end

