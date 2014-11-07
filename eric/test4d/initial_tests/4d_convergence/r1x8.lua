-- Input file for Poisson bracket operator

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.1
-- parameters to control time-stepping
tStart = 0.0
tEnd = 1.0
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 4
tFrame = (tEnd-tStart)/nFrames -- time between frames

-- grid on which equations are to be solved
grid = Grid.RectCart4D {
   lower = {0, 0, 0, 0},
   upper = {1, 1, 1, 1},
   cells = {8, 8, 8, 8},
}

-- create FEM nodal basis
basis = NodalFiniteElement4D.SerendipityElement {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell
   polyOrder = polyOrder,
}

-- number of nodes per cell for DG
numDgNodesPerCell = basis:numNodes()

-- solution
q = DataStruct.Field4D {
   onGrid = grid,
   numComponents = numDgNodesPerCell,
   ghost = {2, 2},
}

-- for RK time-stepping
q1 = DataStruct.Field4D {
   onGrid = grid,
   numComponents = numDgNodesPerCell,
   ghost = {2, 2},
}
-- updated solution
qNew = DataStruct.Field4D {
   onGrid = grid,
   numComponents = numDgNodesPerCell,
   ghost = {2, 2},
}

-- duplicate for use in time-stepping
qNewDup = qNew:duplicate()

function pulse(x,y,z,z1)
   local xc, yc, zc, z1c = 0.5, 0.5, 0.5, 0.5
   local r2 = (x-xc)^2 + (y-yc)^2 + (z-zc)^2 + (z1-z1c)^2
   return math.exp(-75*r2)
end

-- updater to apply initial conditions
initField = Updater.EvalOnNodes4D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,z1,t)
		 return pulse(x,y,z,z1)
	      end
}
initField:setOut( {q} )
-- initialize
initField:advance(0.0) -- time is irrelevant

-- define equation to solve
advectionEqn = HyperEquation.AdvectionTemplated4D {
   -- advection velocity
   speeds = {1.0, 1.0, 1.0, 1.0}
}

-- updater to solve hyperbolic equations
advectSlvr = Updater.NodalDgHyper4D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- equation system to solver
   equation = advectionEqn,
   -- CFL number
   cfl = cfl,
}

-- apply boundary conditions
function applyBc(fld)
   fld:applyPeriodicBc(0)
   fld:applyPeriodicBc(1)
   fld:applyPeriodicBc(2)
   fld:applyPeriodicBc(3)
end

applyBc(q)
qNew:copy(q)

-- write initial conditions
q:write("q_0.h5")

-- solve advection equation
function solveAdvection(curr, dt, qIn, qOut)
   advectSlvr:setCurrTime(curr)
   advectSlvr:setIn( {qIn} )
   advectSlvr:setOut( {qOut} )
   return advectSlvr:advance(curr+dt)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   -- RK stage 1
   local myStatus, myDtSuggested = solveAdvection(tCurr, myDt, q, q1)

   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   applyBc(q1)

   -- RK stage 2
   local myStatus, myDtSuggested = solveAdvection(tCurr, myDt, q1, qNew)

   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   q1:combine(3.0/4.0, q, 1.0/4.0, qNew)
   applyBc(q1)

   -- RK stage 3
   local myStatus, myDtSuggested = solveAdvection(tCurr, myDt, q1, qNew)

   if (myStatus == false) then
      return myStatus, myDtSuggested
   end

   q1:combine(1.0/3.0, q, 2.0/3.0, qNew)
   applyBc(q1)
   q:copy(q1)

   return myStatus, myDtSuggested
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      qNewDup:copy(qNew)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	 myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	 -- time-step too large
	 print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	 qNew:copy(qNewDup)
	 myDt = dtSuggested
      else
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
   q:write( string.format("q_%d.h5", frame) )
end

tCurr = tStart
for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
