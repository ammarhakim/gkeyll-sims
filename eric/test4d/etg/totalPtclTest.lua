-- File to test computation of number density

-- polynomial order
polyOrder = 1
-- physical parameters
elcCharge = -Lucee.ElementaryCharge
ionCharge = Lucee.ElementaryCharge
elcMass = Lucee.ElectronMass
ionMass = 2*Lucee.ProtonMass
elcTemp = 2072 -- [eV]
ionTemp = 2072 -- [ev]
B0 = 1.91 -- [T]
eV = Lucee.ElementaryCharge
R0 = 1.313 -- [m]
a = 0.4701 -- [m]
n0 = 4.992*10^(19) -- [1/m^3]
-- derived parameters
R = R0 + 0.5*a
vtElc = math.sqrt(elcTemp*eV/elcMass)
c_s = math.sqrt(elcTemp*eV/elcMass)
omega_s = math.abs(elcCharge*B0/elcMass)
rho_s = c_s/omega_s
deltaR = 32*rho_s
L_T = R/10
-- grid parameters: number of cells
N_X = 4
N_Y = 4
N_VPARA = 4
N_MU = 2
-- grid parameters: domain extent
X_LOWER = R
X_UPPER = R+deltaR
Y_LOWER = -deltaR/2
Y_UPPER = deltaR/2
VPARA_LOWER = -2*vtElc
VPARA_UPPER = 2*vtElc
MU_LOWER = 0
MU_UPPER = 4*elcMass*VPARA_UPPER*VPARA_UPPER/(2*B0)
-- full 4d phase space grid
grid_4d = Grid.RectCart4D {
   lower = {X_LOWER, Y_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {X_UPPER, Y_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_X, N_Y, N_VPARA, N_MU},
}
-- 2d spatial grid
grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
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

-- create 
basis_4d = NodalFiniteElement4D.SerendipityElement {
   onGrid = grid_4d,
   polyOrder = polyOrder,
}

basis_2d = NodalFiniteElement2D.SerendipityElement {
   onGrid = grid_2d,
   polyOrder = polyOrder,
}

-- distribution function for electrons
f = DataStruct.Field4D {
   onGrid = grid_4d,
   numComponents = basis_4d:numNodes(),
   ghost = {1, 1},
}

function bFieldProfile(x)
  return B0*R/x
end

bField2D = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}

initbField2D = Updater.EvalOnNodes2D {
   onGrid = grid_2d,
   basis = basis_2d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,t)
		 return bFieldProfile(x,y)
	 end
}
runUpdater(initbField2D, 0.0, 0.0, {}, {bField2D})

function fProfile(x,y,v,mu)
  return n0*(2*math.pi*elcTemp*eV/elcMass)^(-3/2)*math.exp(-elcMass*v^2/(2*elcTemp*eV))*math.exp(-mu*bFieldProfile(x,y)/(elcTemp*eV))
end

-- initialize electron distribution function
initElc = Updater.EvalOnNodes4D {
   onGrid = grid_4d,
   basis = basis_4d,
   shareCommonNodes = false,
   -- function to use for initialization
   evaluate = function (x,y,v,mu,t)
		 return fProfile(x,y,v,mu)
	 end
}
runUpdater(initElc, 0.0, 0.0, {}, {f})

numDensityElc = DataStruct.Field2D {
   onGrid = grid_2d,
   numComponents = basis_2d:numNodes(),
   ghost = {1, 1},
}
-- to compute number density
numDensityCalc = Updater.DistFuncMomentCalcWeighted2D {
   -- 4D phase-space grid 
   onGrid = grid_4d,
   -- 4D phase-space basis functions
   basis4d = basis_4d,
   -- 2D spatial basis functions
   basis2d = basis_2d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}

-- dynvector for total particle count
totalPtcl = DataStruct.DynVector { numComponents = 1, }

-- to compute total number of particles in domain
totalPtclCalc = Updater.IntegrateGeneralField2D {
   onGrid = grid_2d,
   basis = basis_2d,
}

function calcDiagnostics(curr, dt)
  runUpdater(numDensityCalc, curr, dt, {f, bField2D}, {numDensityElc})
  numDensityElc:scale(2*math.pi/elcMass)
  runUpdater(totalPtclCalc, curr, dt, {numDensityElc}, {totalPtcl})

end

-- write data to H5 file
function writeFields(frameNum, tCurr)
   f:write( string.format("f_%d.h5", frameNum), tCurr)
   numDensityElc:write( string.format("n_%d.h5", frameNum), tCurr)
end

runUpdater(numDensityCalc, 0, 0, {f, bField2D}, {numDensityElc})
numDensityElc:scale(2*math.pi/elcMass)
runUpdater(totalPtclCalc, 0, 0, {numDensityElc}, {totalPtcl})

print( string.format("Numerical N = %e",totalPtcl:lastInsertedData()) )
print( string.format("Desired N = %e", deltaR*deltaR*n0) )

scaleFactor = deltaR*deltaR*n0/totalPtcl:lastInsertedData()

print ( string.format("Scale factor = %f", scaleFactor) )

-- Scale fields and recalculate density
f:scale(scaleFactor)
runUpdater(numDensityCalc, 0, 0, {f, bField2D}, {numDensityElc})
numDensityElc:scale(2*math.pi/elcMass)
runUpdater(totalPtclCalc, 0, 0, {numDensityElc}, {totalPtcl})
print( string.format("Recomputed Numerical N = %e",totalPtcl:lastInsertedData()) )
