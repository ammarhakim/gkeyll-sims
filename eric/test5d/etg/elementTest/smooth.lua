-- Test advection with PB updater in 5d using serendipity element

-- polynomial order
polyOrder = 1

-- spatial grid
grid = Grid.RectCart3D {
   lower = {0, 0, 0},
   upper = {1, 1, 1},
   cells = {8, 8, 8},
   periodicDirs = {0, 1, 2},
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

basis = NodalFiniteElement3D.SerendipityElement {
   -- grid on which elements should be constructured
   onGrid = grid,
   -- polynomial order in each cell
   polyOrder = polyOrder,
}

-- number of nodes per cell for DG
numDgNodesPerCell = basis:numNodes()

-- solution
q = DataStruct.Field3D {
   onGrid = grid,
   numComponents = numDgNodesPerCell,
   ghost = {1, 1},
}

qSmoothed = q:duplicate()

-- updater to apply initial conditions
initField = Updater.ProjectOnNodalBasis3D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function (x,y,z,t)
		 return math.sin(4*math.pi*y)
	 end
}
initField:setOut( {q} )
-- initialize
initField:advance(0.0) -- time is irrelevant

q:sync()

-- Smooth out initial field
smoothCalc = Updater.SimpleSmoothToC03D {
   onGrid = grid,
   basis = basis,
}

runUpdater(smoothCalc, 0.0, 0.0, {q}, {qSmoothed})

q:write( string.format("q_0.h5"), 0.0)
qSmoothed:write( string.format("qS_0.h5"), 0.0)
