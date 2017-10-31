-- Input file to solve 1D Poisson equations using FEM

-- polynomial order
polyOrder = 1

-- grid to solve equations on
grid = Grid.RectCart1D {
   lower = {0.0},
   upper = {1.0},
   cells = {36},
}

-- create FEM nodal basis
basis = NodalFiniteElement1D.Lobatto {
   onGrid = grid,
   polyOrder = polyOrder,
}

-- source term
src = DataStruct.Field1D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- number density
numDens = DataStruct.Field1D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- function to initialize source
function initSrcFunc(x,y,z)
   local n = 4*x*(1-x)
   local phiP, phiPP = 4*(1-2*x), -8
   local nP = 4*(1-2*x)
   return n*phiPP+nP*phiP
end

-- updater to apply initial conditions
initSrc = Updater.EvalOnNodes1D {
   onGrid = grid,
   basis = basis,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      return initSrcFunc(x,y,z)
   end
}
initSrc:setOut( {src} )
initSrc:advance(0.0) -- time is irrelevant
src:write("src.h5")

-- updater to apply initial conditions
initNumDens = Updater.EvalOnNodes1D {
   onGrid = grid,
   basis = basis,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local n0 = 1.0
      return 4*n0*x*(1-x)
   end
}
initNumDens:setOut( {numDens} )
initNumDens:advance(0.0) -- time is irrelevant
numDens:write("numDens.h5")

-- field to store potential
phi = DataStruct.Field1D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
phi:clear(0.0)

-- create updater to solve Poisson equation
poissonSlvr = Updater.FemGenPoisson1D {
   onGrid = grid,
   basis = basis,
   sourceNodesShared = false,
   solutionNodesShared = false,
   bcLeft = { T = "D", V = 0.0 },
   bcRight = { T = "D", V = 0.0 },
}

-- set input output fields
poissonSlvr:setIn( {src, numDens} )
poissonSlvr:setOut( {phi} )

-- solve for potential (time is irrelevant here)
status, dtSuggested = poissonSlvr:advance(0.0)
-- check if solver converged
if (status == false) then
   Lucee.logError("Poisson solver failed to converge!")
end

-- output solution
phi:write("phi.h5")
Lucee.logError(string.format("Poisson solver took %g ...", poissonSlvr:totalAdvanceTime()))

