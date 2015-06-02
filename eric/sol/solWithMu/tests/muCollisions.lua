-- testing diffusion updaters
-- 5-31-2015: Used to test collisions in mu only

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.05

-- physical constants
-- eletron mass (kg)
electronMass = Lucee.ElectronMass
-- electron volt to joules
eV = Lucee.ElementaryCharge
-- Deuterium ion mass (kg)
ionMass = 2.014*Lucee.ProtonMass
-- signed ion charge
ionCharge = Lucee.ElementaryCharge

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
-- Magnetic field
B0 = 1

-- Derived parameters
-- Thermal velocity of pedestal
vtPed = math.sqrt(tPed*eV/ionMass)
-- Pedestal sound speed (m/s)
cPed = math.sqrt(2*tPed*eV/ionMass)
-- Particle source
Sn   = A*nPed*cPed/lSource

-- number of cells
N_Z, N_VPARA, N_MU = 8, 16, 8
-- domain extents
Z_LOWER, Z_UPPER = -lParallel, lParallel
VPARA_UPPER = math.min(4, 2.5*math.sqrt(N_VPARA/4))*vtPed
--VPARA_UPPER = 6*vtPed
VPARA_LOWER = -VPARA_UPPER
MU_LOWER = 0
MU_UPPER = math.min(8, 4*math.sqrt(N_MU/2))*ionMass*vtPed*vtPed/B0

k = math.pi/VPARA_UPPER

-- parameters to control time-stepping
tStart = 0.0
tEnd = 1.0
nFrames = 5

-- phase space grid 
grid = Grid.RectCart3D {
   lower = {Z_LOWER, VPARA_LOWER, MU_LOWER},
   upper = {Z_UPPER, VPARA_UPPER, MU_UPPER},
   cells = {N_Z, N_VPARA, N_MU},
}

-- create FEM nodal basis
basis = NodalFiniteElement3D.SerendipityElement {
   onGrid = grid,
   polyOrder = polyOrder,
}

-- spatial grid
grid_1d = Grid.RectCart1D {
   lower = {Z_LOWER},
   upper = {Z_UPPER},
   cells = {N_Z},
}

-- spatial FEM nodal basis
basis_1d = NodalFiniteElement1D.LagrangeTensor {
   onGrid = grid_1d,
   polyOrder = polyOrder,
}

function maxwellian(mass, vt, vPara, mu)
  local vtPerp = math.sqrt(tPed*eV/ionMass)
  return 1/(2*Lucee.Pi*vtPerp^2*math.sqrt(2*Lucee.Pi*vt^2))*math.exp(-vPara^2/(2*vt^2))*
    math.exp(-mu*B0/(tPed*eV))
end

-- Maxwellian (Right half only)
function maxwellianRight(mass, vt, vPara, mu)
  local vtPerp = math.sqrt(tPed*eV/ionMass)
  if vPara >=0 then
    return 1/(2*Lucee.Pi*vtPerp^2*math.sqrt(2*Lucee.Pi*vt^2))*math.exp(-vPara^2/(2*vt^2))*
      math.exp(-mu*B0/(tPed*eV))
  else
    return 0
  end
end

-- Maxwellian (Left half only)
function maxwellianLeft(mass, vt, vPara, mu)
  local vtPerp = math.sqrt(tPed*eV/ionMass)
  if vPara <= 0 then
    return 1/(2*Lucee.Pi*vtPerp^2*math.sqrt(2*Lucee.Pi*vt^2))*math.exp(-vPara^2/(2*vt^2))*
      math.exp(-mu*B0/(tPed*eV))
  else
    return 0
  end
end

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

-- distribution function
distf = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
distf:clear(0.0)

-- extra fields for performing RK update
distfNew = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
distf1 = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
distfDrag = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- updater to initialize distribution function
initDistf = Updater.ProjectOnNodalBasis3D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(z,vPara,mu,t)
     local vTherm = math.sqrt(tPed*eV/ionMass)

     if mu < 4*ionMass*vTherm^2/(2*B0) then
       return 1
     else
       return 0
     end

     --return maxwellian(ionMass, vTherm, vPara, mu)
   end
}
runUpdater(initDistf, 0.0, 0.0, {}, {distf})

diffSlvr = Updater.LenardBernsteinDiff3DUpdater {
  onGrid = grid,
  basis = basis,
  cfl = cfl,
  onlyIncrement = false,
  diffusionCoeff = 1.0,
  B0 = B0,
  speciesMass = ionMass,
}

dragSlvr = Updater.LenardBernsteinDrag3DUpdater {
  onGrid = grid,
  basis = basis,
  basis1d = basis_1d,
  cfl = cfl,
  onlyIncrement = true,
  diffusionCoeff = 1.0,
}

vThermSq = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
vThermSq3d = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- updater to initialize distribution function
initVThermSq = Updater.EvalOnNodes1D {
   onGrid = grid_1d,
   basis = basis_1d,
   shareCommonNodes = false,
   evaluate = function(z,t)
     return vtPed*vtPed
   end
}
runUpdater(initVThermSq, 0.0, 0.0, {}, {vThermSq})

copy1dTo3d = Updater.NodalCopy1DTo3DFieldUpdater {
  onGrid = grid,
  polyOrder = polyOrder,
  basis1d = basis_1d,
  basis3d = basis
}
runUpdater(copy1dTo3d, 0.0, 0.0, {vThermSq}, {vThermSq3d})

uParallel = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- updater to initialize distribution function
initUParallel = Updater.EvalOnNodes1D {
   onGrid = grid_1d,
   basis = basis_1d,
   shareCommonNodes = false,
   evaluate = function(z,t)
     return 0
   end
}
runUpdater(initUParallel, 0.0, 0.0, {}, {uParallel})

-- apply boundary conditions to 2-D field
function applyBc(curr, dt, fldIon)
   --for i,bc in ipairs({bcLowerIon, bcUpperIon}) do
   --   runUpdater(bc, curr, dt, {}, {fldIon})
   --end
  --fldIon:applyPeriodicBc(0)
  --fldIon:applyPeriodicBc(1)
  --fldIon:applyPeriodicBc(2)
end

-- compute various diagnostics
function calcDiagnostics(curr, dt)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local diffStatus, diffDtSuggested
   local dragStatus, dragDtSuggested

   -- RK stage 1
   diffStatus, diffDtSuggested = runUpdater(diffSlvr, tCurr, myDt, {distf, vThermSq3d}, {distf1})
   dragStatus, dragDtSuggested = runUpdater(dragSlvr, tCurr, myDt, {distf, uParallel}, {distfDrag})
   distf1:accumulate(myDt, distfDrag)

   if (diffStatus == false or dragStatus == false) then
     return false, math.min(diffDtSuggested, dragDtSuggested)
   end  

   applyBc(tCurr, myDt, distf1)

   -- RK stage 2
   diffStatus, diffDtSuggested = runUpdater(diffSlvr, tCurr, myDt, {distf1, vThermSq3d}, {distfNew})
   dragStatus, dragDtSuggested = runUpdater(dragSlvr, tCurr, myDt, {distf1, uParallel}, {distfDrag})
   distfNew:accumulate(myDt, distfDrag)
   
   if (diffStatus == false or dragStatus == false) then
     return false, math.min(diffDtSuggested, dragDtSuggested)
   end  

   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)
   applyBc(tCurr, myDt, distf1)


   -- RK stage 3
   diffStatus, diffDtSuggested = runUpdater(diffSlvr, tCurr, myDt, {distf1, vThermSq3d}, {distfNew})
   dragStatus, dragDtSuggested = runUpdater(dragSlvr, tCurr, myDt, {distf1, uParallel}, {distfDrag})
   distfNew:accumulate(myDt, distfDrag)
   
   if (diffStatus == false or dragStatus == false) then
     return false, math.min(diffDtSuggested, dragDtSuggested)
   end  
   
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)
   applyBc(tCurr, myDt, distf1)
   
   distf:copy(distf1)

   return true, math.min(diffDtSuggested, dragDtSuggested)
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
	      --calcDiagnostics(tCurr, myDt)

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
function writeFields(frameNum, tCurr)
   distf:write( string.format("distfIon_%d.h5", frameNum), tCurr)
end

applyBc(0.0, 0.0, distf)
-- compute initial diagnostics
--calcDiagnostics(0.0, 0.0)
-- write out initial conditions
writeFields(0, 0.0)
-- make a duplicate in case we need it
distfDup = distf:duplicate()

-- parameters to control time-stepping
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
tFrame = (tEnd-tStart)/nFrames
tCurr = tStart

for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame, tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
