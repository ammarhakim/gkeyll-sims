-- Swirl test for PB operator in 4d
-- Swirl motion is in x-y plane so pbEqn had to be fudged to be
-- (0,1,0,0; -1,0,0,0; 0,0,0,0; 0,0,0,0) for this test

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.1

-- period for distortion
Tperiod = 1.5

-- grid on which equations are to be solved
grid = Grid.RectCart4D {
   lower = {0, 0, 0, 0},
   upper = {1.0, 1.0, 1.0, 1.0},
   cells = {32, 32, 2, 2},
}
-- spatial grid
grid_2d = Grid.RectCart2D {
   lower = {0, 0},
   upper = {1.0, 1.0},
   cells = {32, 32},
}

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

-- create FEM nodal basis
basis = NodalFiniteElement4D.SerendipityElement {
   onGrid = grid,
   polyOrder = polyOrder,
}
basis_2d = NodalFiniteElement2D.SerendipityElement {
   onGrid = grid_2d,
   polyOrder = polyOrder,
}

-- vorticity
chi = DataStruct.Field4D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {2, 2},
}
-- clear out contents
chi:clear(0.0)

-- extra fields for performing RK update
chiNew = DataStruct.Field4D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {2, 2},
}
chi1 = DataStruct.Field4D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {2, 2},
}

-- potential
phi = DataStruct.Field4D {
   onGrid = grid,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis:numNodes(),
   ghost = {2, 2},
}
-- to store initial condition
phiStatic = DataStruct.Field4D {
   onGrid = grid,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis:numNodes(),
   ghost = {2, 2},
}

-- create updater to initialize potential
initPhi = Updater.EvalOnNodes4D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,w,t)
		 local mpi = Lucee.Pi
		 local t1 = math.sin(mpi*x)*math.sin(mpi*y)
		 return t1^2/mpi
	      end
}
initPhi:setOut( {phiStatic} )
-- initialize potential
initPhi:advance(0.0) -- time is irrelevant

-- copy over initial condition
phi:copy(phiStatic)

-- create updater to initialize vorticity
initChi = Updater.EvalOnNodes4D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,w,t)
     local x0, y0, r0 = 0.25, 0.5, 0.15
     local r = math.min(math.sqrt((x-x0)^2+(y-y0)^2), r0)/r0
     local chump = 0.25*(1+math.cos(Lucee.Pi*r))

     local xc, yc, rc = 0.75, 0.5, 0.15
     local r2 = math.min(math.sqrt((x-xc)^2+(y-yc)^2), rc)/rc
     local chump2 = 0.25*(1+math.cos(Lucee.Pi*r2))

     return chump + chump2
	 end
}
initChi:setOut( {chi} )
-- initialize potential
initChi:advance(0.0) -- time is irrelevant

-- define equation to solve
pbEqn = PoissonBracketEquation.Canonical4D {
}
-- create updater for Poisson bracket
pbSlvr = Updater.PoissonBracket4D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
   -- equation to solve
   equation = pbEqn,
}

-- write initial value
chi:write("chi_0.h5")

-- function to apply boundary conditions
function applyBc(fld)
   fld:applyCopyBc(0, "lower")
   fld:applyCopyBc(0, "upper")
   fld:applyCopyBc(1, "lower")
   fld:applyCopyBc(1, "upper")
   fld:applyCopyBc(2, "lower")
   fld:applyCopyBc(2, "upper")
   fld:applyCopyBc(3, "lower")
   fld:applyCopyBc(3, "upper")
end

-- apply BCs to initial conditions
applyBc(chi)

-- to store number density
numDensity = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {2, 2},
}
-- to compute number density
numDensityCalc = Updater.DistFuncMomentCalc2D {
   -- 4D phase-space grid 
   onGrid = grid,
   -- 4D phase-space basis functions
   basis4d = basis,
   -- 2D spatial basis functions
   basis2d = basis_2d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}

function calcDiagnostics(curr, dt)
  runUpdater(numDensityCalc, curr, dt, {chi}, {numDensity})
end

-- function to solve poisson bracket
function poissonBracket(curr, dt, chiIn, phiIn, chiOut)
   pbSlvr:setCurrTime(curr)
   pbSlvr:setIn( {chiIn, phiIn} )
   pbSlvr:setOut( {chiOut} )
   return pbSlvr:advance(curr+dt)
end

-- function to solve for potential
function poisson(curr, dt, chiIn, phiOut)
   -- this is a time-dependent potential and not a true Poisson solve
   phiOut:combine(math.cos(curr*Lucee.Pi/Tperiod), phiStatic)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local status, dtSuggested

   -- RK stage 1 (chi1 <- chi + L(chi))
   status, dtSuggested = poissonBracket(tCurr, myDt, chi, phi, chi1)

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   -- apply BCs
   applyBc(chi1)

   -- solve for potential
   poisson(tCurr+myDt, myDt, chi1, phi)

   -- RK stage 2
   status, dtSuggested = poissonBracket(tCurr, myDt, chi1, phi, chiNew)

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   chi1:combine(3.0/4.0, chi, 1.0/4.0, chiNew)

   -- apply BCs
   applyBc(chi1)

   -- solve for potential
   poisson(tCurr+myDt, myDt, chi1, phi)

   -- RK stage 3
   status, dtSuggested = poissonBracket(tCurr, myDt, chi1, phi, chiNew)

   -- check if step failed and return immediately if it did
   if (status == false) then
      return status, dtSuggested
   end

   chi1:combine(1.0/3.0, chi, 2.0/3.0, chiNew)
   -- apply BCs
   applyBc(chi1)

   -- copy over solution
   chi:copy(chi1)

   -- solve for potential
   poisson(tCurr+myDt, myDt, chi, phi)

   return status, dtSuggested
end

-- make a duplicate in case we need it
chiDup = chi:duplicate()

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   -- declare local variables
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   -- main loop
   while tCurr<=tEnd do
      -- copy chi over
      chiDup:copy(chi)

      -- if needed adjust dt to hit tEnd exactly
      if (tCurr+myDt > tEnd) then
	      myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))

      -- take a time-step
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
         -- time-step too large
         print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
         -- copy in case current solutions were messed up
         chi:copy(chiDup)
         myDt = dtSuggested
      else
	      tCurr = tCurr + myDt
	      myDt = dtSuggested
	      step = step + 1
	      -- check if done
	      if (tCurr >= tEnd) then
	        break
	      end
      end
   end

   return dtSuggested
end

-- parameters to control time-stepping
tStart = 0.0
tEnd = Tperiod
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
nFrames = 40
tFrame = (tEnd-tStart)/nFrames -- time between frames

-- write data to H5 file
function writeFields(frame)
   chi:write( string.format("chi_%d.h5", frame) )
   numDensity:write( string.format("n_%d.h5", frame) )
end

calcDiagnostics(0.0, 0.0)

tCurr = tStart
for frame = 1, nFrames do

   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   -- advance solution between frames
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   calcDiagnostics(tCurr, tFrame)
   writeFields(frame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
