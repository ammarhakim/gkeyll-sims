-- lua file to compute something proportional to E_theta
-- mpiexec -n 4 /Users/eshi/Research/gkeyllall/par-opt/gkeyll/gkeyll -i test117.lua -pc_type lu -pc_factor_mat_solver_package superlu_dist

-- phase-space decomposition (5d)
phaseDecomp = DecompRegionCalc5D.CartProd { cuts = {4, 1, 1, 1, 1} }
-- configuration space decomposition (3d)
confDecomp = DecompRegionCalc3D.SubCartProd5D {
   decomposition = phaseDecomp,
   collectDirections = {0, 1, 2},
}
-- logical sheath decomposition for solve on 2d planes
sheathDecomp = DecompRegionCalc2D.SubCartProd3D {
   decomposition = confDecomp,
   collectDirections = {0, 1},
}

-- polynomial order
polyOrder = 1

-- cfl number to use
cfl = 0.1
-- parameters to control time-stepping
tStart = 0.0
tEnd = 5832e-6
nFrames = 5832 
tFrame = (tEnd-tStart)/nFrames -- time between frames
tCurr = tStart

-- physical parameters
eV        = Lucee.ElementaryCharge
elcCharge = -eV
ionCharge = eV
eps0      = Lucee.Epsilon0 -- permittivity of free space
ionMass   = 3.973*Lucee.ProtonMass 
elcMass   = ionMass/400 -- REDUCED ELECTRON MASS Lucee.ElectronMass
elcTemp   = 6 -- [eV]
ionTemp   = 1 -- [eV]
B0 = 0.0398  -- [T]
n0 = 2*10^18 -- [1/m^3]
-- derived parameters
vtElc   = math.sqrt(0.5*elcTemp*eV/elcMass) -- only used for grid
vtIon   = math.sqrt(ionTemp*eV/ionMass)
omega_i = math.abs(ionCharge*B0/ionMass)

c_s = math.sqrt(elcTemp*eV/ionMass)
rho_s = c_s/omega_i

deltaX  = 100*rho_s
deltaY  = 100*rho_s
R          = 40*rho_s -- [m]
L_parallel = 36*R -- [m]
-- grid parameters: number of cells
N_X = 36
N_Y = 36
N_Z = 10
N_VPARA = 10
N_MU = N_VPARA/2
-- grid parameters: domain extent
X_LOWER = -deltaX/2
X_UPPER = deltaX/2
Y_LOWER = -deltaY/2
Y_UPPER = deltaY/2
Z_LOWER = -L_parallel/2
Z_UPPER = L_parallel/2
-- source parameters
L_s = 0.5*rho_s
r_s = 20*rho_s
sourceAmplitude = 0.03*n0*c_s/R
-- total volume
totalVol = (X_UPPER-X_LOWER)*(Y_UPPER-Y_LOWER)*(Z_UPPER-Z_LOWER)

VPARA_UPPER_ELC = 4*vtElc
VPARA_LOWER_ELC = -VPARA_UPPER_ELC
MU_LOWER_ELC = 0
MU_UPPER_ELC = 0.75*elcMass*(VPARA_UPPER_ELC*VPARA_UPPER_ELC)/(2*B0)

VPARA_UPPER_ION = 4*vtIon
VPARA_LOWER_ION = -VPARA_UPPER_ION
MU_LOWER_ION = 0
MU_UPPER_ION = 0.75*ionMass*(VPARA_UPPER_ION*VPARA_UPPER_ION)/(2*B0)

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

-- full 5d phase space grids
grid_elc_5d = Grid.RectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER_ELC, MU_LOWER_ELC},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER_ELC, MU_UPPER_ELC},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
   decomposition = phaseDecomp,
}
grid_ion_5d = Grid.RectCart5D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER, VPARA_LOWER_ION, MU_LOWER_ION},
   upper = {X_UPPER, Y_UPPER, Z_UPPER, VPARA_UPPER_ION, MU_UPPER_ION},
   cells = {N_X, N_Y, N_Z, N_VPARA, N_MU},
   decomposition = phaseDecomp,
}

-- 3d spatial grid
grid_3d = Grid.RectCart3D {
   lower = {X_LOWER, Y_LOWER, Z_LOWER},
   upper = {X_UPPER, Y_UPPER, Z_UPPER},
   cells = {N_X, N_Y, N_Z},
   decomposition = confDecomp,
}

-- 2d spatial grid for sheath boundary conditions
grid_2d = Grid.RectCart2D {
   lower = {X_LOWER, Y_LOWER},
   upper = {X_UPPER, Y_UPPER},
   cells = {N_X, N_Y},
   decomposition = sheathDecomp,
}

-- create 5d basis functions
basis_elc_5d = NodalFiniteElement5D.SerendipityElement {
   onGrid = grid_elc_5d,
   polyOrder = polyOrder,
}
basis_ion_5d = NodalFiniteElement5D.SerendipityElement {
   onGrid = grid_ion_5d,
   polyOrder = polyOrder,
}
-- create 3d basis functions
basis_3d = NodalFiniteElement3D.SerendipityElement {
   onGrid = grid_3d,
   polyOrder = polyOrder,
}

-- stores result of potential solve
phi = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
phi_x = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
phi_y = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

eTheta = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

-- to compute y-derivative of a field
yDerivCalc = Updater.SOLDerivativeCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
  dir = 1,
}
-- to compute x-derivative of a field
xDerivCalc = Updater.SOLDerivativeCalc3D {
  onGrid = grid_3d,
  basis = basis_3d,
  dir = 0,
}

xField = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
initXField = Updater.EvalOnNodes3D {
   onGrid = grid_3d,
   basis = basis_3d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,t)
      return x
   end
}
runUpdater(initXField, 0, 0, {}, {xField})

yField = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
initYField = Updater.EvalOnNodes3D {
   onGrid = grid_3d,
   basis = basis_3d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,t)
      return y
   end
}
runUpdater(initYField, 0, 0, {}, {yField})

rField = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}
initRField = Updater.EvalOnNodes3D {
   onGrid = grid_3d,
   basis = basis_3d,
   shareCommonNodes = false,
   evaluate = function (x,y,z,t)
      local rVal = math.sqrt(x^2+y^2)
      if rVal > 0 then
        return 1/rVal
      else
        return 0
      end
   end
}
runUpdater(initRField, 0, 0, {}, {rField})

-- point-wise multiplication of two fields
multiply3dCalc = Updater.FieldArithmeticUpdater3D {
  onGrid = grid_3d,
  basis = basis_3d,
  evaluate = function(fld1, fld2)
    return fld1*fld2
  end
}

tempField = DataStruct.Field3D {
   onGrid = grid_3d,
   numComponents = basis_3d:numNodes(),
   ghost = {1, 1},
}

tCurr = tStart + tFrame*4000

for frame = 4000, nFrames do
  Lucee.logInfo (string.format("-- Frame %g", frame))
  Lucee.logInfo ("")
  -- read in total field at this time
  phi:read("phi_" .. frame .. ".h5")
  runUpdater(xDerivCalc, 0, 0, {phi}, {phi_x})
  runUpdater(yDerivCalc, 0, 0, {phi}, {phi_y})

  eTheta:clear(0.0)
  runUpdater(multiply3dCalc, 0, 0, {phi_y,xField}, {tempField})
  eTheta:accumulate(-1.0, tempField)
  runUpdater(multiply3dCalc, 0, 0, {phi_x,yField}, {tempField})
  eTheta:accumulate(1.0, tempField)
  -- multiply by 1/r (note that r=0 points are to be neglected)
  runUpdater(multiply3dCalc, 0, 0, {eTheta,rField}, {tempField})
  eTheta:copy(tempField)
  -- write out E_theta
  eTheta:write( string.format("eTheta_%d.h5", frame), tCurr)
  -- advance time
  tCurr = tCurr + tFrame
end
