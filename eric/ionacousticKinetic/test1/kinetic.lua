-- Input file for ion acoustic problem with kinetic ions and electrons
-- Using a different initial condition for electons to try to minimze high
-- frequency oscillations

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.1

-- wave-number
knumber = 0.5

-- initial number density in each cell
initNumDens = 1.0
-- temperature ratio (T_i/T_e)
Tratio = 0.5

kPerpTimesRho = 0.2
-- ion temperature
ionTemp = 1.0
-- electron temperature
elcTemp = ionTemp/Tratio
-- signed electron charge
elcCharge = -1.0
-- signed ion charge
ionCharge = 1.0
-- electron mass
elcMass = Lucee.ElectronMass/Lucee.ProtonMass
-- ion mass
ionMass = 1
-- electron drift speed
elcDrift = 0.0
-- ion drift speed
ionDrift = 0.0

-- permittivity of free space
epsilon0 = 1.0

-- domain extents
XL, XU = -Lucee.Pi/knumber, Lucee.Pi/knumber
VL, VU = -6.0, 6.0
-- number of cells
NX, NV = 32, 32
-- compute max thermal speed to set velocity space extents
vtElc = math.sqrt(elcTemp/elcMass)
VL_ELC, VU_ELC = -6.0*vtElc, 6.0*vtElc

vtIon = math.sqrt(ionTemp/ionMass)
VL_ION, VU_ION = -6.0*vtIon, 6.0*vtIon

-- parameters to control time-stepping
tStart = 0.0
tEnd = 20.0
nFrames = 10

-- Determine number of global nodes per cell for use in creating CG
-- fields. Note that this looks a bit odd as this not the number of
-- *local* nodes but the number of nodes in each cell to give the
-- correct number of global nodes in fields.
if (polyOrder == 1) then
   numCgNodesPerCell = 1
   numCgNodesPerCell_1d = 1
elseif (polyOrder == 2) then
   numCgNodesPerCell = 3
   numCgNodesPerCell_1d = 2
elseif (polyOrder == 2) then
   numCgNodesPerCell = 5
   numCgNodesPerCell_1d = 3
end

-- phase space grid for electrons
gridElc = Grid.RectCart2D {
   lower = {XL, VL_ELC},
   upper = {XU, VU_ELC},
   cells = {NX, NV},
}
gridIon = Grid.RectCart2D {
   lower = {XL, VL_ION},
   upper = {XU, VU_ION},
   cells = {NX, NV},
}

-- create FEM nodal basis
basisElc = NodalFiniteElement2D.Serendipity {
   onGrid = gridElc,
   polyOrder = polyOrder,
}

basisIon = NodalFiniteElement2D.Serendipity {
   onGrid = gridIon,
   polyOrder = polyOrder,
}

-- distribution function for electrons
distfElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
-- distribution function for ions
distfIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numNodes(),
   ghost = {1, 1},
}

-- Maxwellian with number density 'n0', drift-speed 'vdrift' and
-- thermal speed 'vt' = \sqrt{T/m}, where T and m are species
-- temperature and mass respectively.
function maxwellian(n0, vdrift, vt, v)
   return n0/math.sqrt(2*Lucee.Pi*vt^2)*math.exp(-(v-vdrift)^2/(2*vt^2))
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

-- updater to initialize distribution function
initDistfElc = Updater.EvalOnNodes2D {
   onGrid = gridElc,
   -- basis functions to use
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		  local ionThermal = math.sqrt(ionTemp/ionMass)
      local alpha = 0.01 -- perturbation
		  local k = knumber
		  return ((1+alpha*math.cos(k*x))*maxwellian(initNumDens, ionDrift, ionThermal, y)
		    + kPerpTimesRho^2*initNumDens)/(1+kPerpTimesRho^2)
	   end
}
runUpdater(initDistfElc, 0.0, 0.0, {}, {distfElc})

-- updater to initialize distribution function
initDistfIon = Updater.EvalOnNodes2D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(x,y,z,t)
		 local ionThermal = math.sqrt(ionTemp/ionMass)
		 local alpha = 0.01 -- perturbation
		 local k = knumber
		 return (1+alpha*math.cos(k*x))*maxwellian(initNumDens, ionDrift, ionThermal, y)
	  end
}
runUpdater(initDistfIon, 0.0, 0.0, {}, {distfIon})

-- extra fields for performing RK update
distfNewElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distf1Elc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distfNewIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}
distf1Ion = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisElc:numNodes(),
   ghost = {1, 1},
}

-- Electron Hamiltonian
hamilElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
-- Ion Hamiltonian
hamilIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- create field to store kinetic energy term in Hamiltonian
hamilKeElc = DataStruct.Field2D {
   onGrid = gridElc,
   numComponents = basisElc:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}
hamilKeIon = DataStruct.Field2D {
   onGrid = gridIon,
   numComponents = basisIon:numExclusiveNodes(),
   ghost = {1, 1},
   -- ghost cells to write
   writeGhost = {0, 1} -- write extra layer on right to get nodes
}

-- updater to initialize electron kinetic energy term in Hamiltonian
initHamilKeElc = Updater.EvalOnNodes2D {
   onGrid = gridElc,
   basis = basisElc,
   -- are common nodes shared?
   shareCommonNodes = true, -- Hamiltonian is continuous
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}
runUpdater(initHamilKeElc, 0.0, 0.0, {}, {hamilKeElc})

-- updater to initialize ion kinetic energy term in Hamiltonian
initHamilKeIon = Updater.EvalOnNodes2D {
   onGrid = gridIon,
   basis = basisIon,
   -- are common nodes shared?
   shareCommonNodes = true, -- Hamiltonian is continuous
   -- function to use for initialization
   evaluate = function (x,y,z,t)
      local v = y
      return v^2/2
   end
}
runUpdater(initHamilKeIon, 0.0, 0.0, {}, {hamilKeIon})

-- Updater for electron Vlasov equation
vlasovSolverElc = Updater.PoissonBracket {
   onGrid = gridElc,
   basis = basisElc,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
}
vlasovSolverIon = Updater.PoissonBracket {
   onGrid = gridIon,
   basis = basisIon,
   -- cfl number to use
   cfl = cfl,
   -- flux type: one of "upwind" (default) or "central"
   fluxType = "upwind",
}

-- spatial grid
grid_1d = Grid.RectCart1D {
   lower = {XL},
   upper = {XU},
   cells = {NX},
}

-- spatial FEM nodal basis
basis_1d = NodalFiniteElement1D.LagrangeTensor {
   onGrid = grid_1d,
   polyOrder = polyOrder,
}

-- Electron number density
numDensityElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Ion number density
numDensityIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Electron momentum
momentumElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Ion momentum
momentumIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Electron particle energy
ptclEnergyElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- Ion particle energy
ptclEnergyIon = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- charge density
chargeDensity = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- Updater to compute electron number density
numDensityCalcElc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridElc,
   -- 2D phase-space basis functions
   basis2d = basisElc,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}
numDensityCalcIon = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridIon,
   -- 2D phase-space basis functions
   basis2d = basisIon,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 0,
}

-- Updater to compute electron momentum
momentumCalcElc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridElc,
   -- 2D phase-space basis functions
   basis2d = basisElc,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 1,
}
momentumCalcIon = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridIon,
   -- 2D phase-space basis functions
   basis2d = basisIon,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 1,
}

-- Updater to compute electron momentum
ptclEnergyCalcElc = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridElc,
   -- 2D phase-space basis functions
   basis2d = basisElc,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 2,
}
ptclEnergyCalcIon = Updater.DistFuncMomentCalc1D {
   -- 2D phase-space grid 
   onGrid = gridIon,
   -- 2D phase-space basis functions
   basis2d = basisIon,
   -- 1D spatial basis functions
   basis1d = basis_1d,
   -- desired moment (0, 1 or 2)
   moment = 2,
}

-- calculate number density at current time
function calcNumDensity(calculator, curr, dt, distfIn, numDensOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {numDensOut})
end

-- calculate momentum at current time
function calcMomentum(calculator, curr, dt, distfIn, momentumOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {momentumOut})
end

-- calculate energy
function calcPtclEnergy(calculator, curr, dt, distfIn, energyOut)
   return runUpdater(calculator, curr, dt, {distfIn}, {energyOut})
end

-- compute moments from distribution function
function calcMoments(curr, dt, distfElcIn, distfIonIn)
  -- number density
  calcNumDensity(numDensityCalcElc, curr, dt, distfElcIn, numDensityElc)
  calcNumDensity(numDensityCalcIon, curr, dt, distfIonIn, numDensityIon)

  -- momentum
  calcMomentum(momentumCalcElc, curr, dt, distfElcIn, momentumElc)
  calcMomentum(momentumCalcIon, curr, dt, distfIonIn, momentumIon)

  -- energy
  calcPtclEnergy(ptclEnergyCalcElc, curr, dt, distfElcIn, ptclEnergyElc)
  calcPtclEnergy(ptclEnergyCalcIon, curr, dt, distfIonIn, ptclEnergyIon)
end

-- field to store continous potential in 1D
phi1d = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numExclusiveNodes(),
   ghost = {1, 1},
   -- write ghosts
   writeGhost = {0, 1},
}
phi1dDg = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- updater to copy 1D field to 2D field
copyTo2DElc = Updater.Copy1DTo2DNodalField {
   -- grid for updater
   onGrid = gridElc,
}
copyTo2DIon = Updater.Copy1DTo2DNodalField {
   -- grid for updater
   onGrid = gridIon,
}

-- function to copy 1D field to 2D field
function copyPhi(copier, curr, dt, phi1, phi2)
   return runUpdater(copier, curr, dt, {phi1}, {phi2})
end

-- update Vlasov equation, given appropriate updater 
function updateVlasovEqn(pbSlvr, curr, dt, distfIn, hamilIn, distfOut)
   return runUpdater(pbSlvr, curr, dt, {distfIn, hamilIn}, {distfOut})
end

-- compute phi from number density
function calcPhiFromChargeDensity(curr, dt, distElcIn, distIonIn, phiOut)
   calcMoments(curr, dt, distElcIn, distIonIn)
   -- charge density: rho_c/epsilon0
   chargeDensity:combine(1.0/epsilon0, numDensityIon, -1.0/epsilon0, numDensityElc)
   chargeDensity:scale(elcTemp/(ionCharge*kPerpTimesRho^2*initNumDens))
   phiOut:copy(chargeDensity)
end

-- create updater to initialize chi
copyCToD = Updater.CopyContToDisCont1D {
   onGrid = grid_1d,
   basis = basis_1d,
}

-- updater to move phi to continuous field
phiToContCalc = Updater.ContFromDisCont1D {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
}

function copyPotential(tCurr, dt, cgIn, dgOut)
   runUpdater(copyCToD, tCurr, dt, {cgIn}, {dgOut})
end

-- compute continuous phi from discontinuous phi
function calcContinuousPhi(curr, dt, phiIn, phiOut)
   phiOut:clear(0.0)
   phiToContCalc:setCurrTime(curr)
   phiToContCalc:setIn( {phiIn} )
   phiToContCalc:setOut( {phiOut} )
   return phiToContCalc:advance(curr+dt)
end

-- compute hamiltonian for electrons
function calcHamiltonianElc(curr, dt, phiIn, hamilOut)
   hamilOut:clear(0.0)
   copyPhi(copyTo2DElc, curr, dt, phiIn, hamilOut)
   hamilOut:scale(elcCharge/elcMass)
   hamilOut:accumulate(1.0, hamilKeElc)
end
-- compute hamiltonian for ions
function calcHamiltonianIon(curr, dt, phiIn, hamilOut)
   hamilOut:clear(0.0)
   copyPhi(copyTo2DIon, curr, dt, phiIn, hamilOut)
   hamilOut:scale(ionCharge/ionMass)
   hamilOut:accumulate(1.0, hamilKeIon)
end

-- apply boundary conditions
function applyBc(curr, dt, fldElc, fldIon)
   for i,fld in ipairs({fldElc, fldIon}) do
      fld:applyPeriodicBc(0)
      fld:applyCopyBc(1, "lower")
      fld:applyCopyBc(1, "upper")
   end
end

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

-- compute various diagnostics
function calcDiagnostics(curr, dt)
   fieldEnergyCalc:setCurrTime(curr)
   fieldEnergyCalc:advance(curr+dt)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local statusElc, dtSuggestedElc
   local statusIon, dtSuggestedIon

   -- RK stage 1
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distfElc, hamilElc, distf1Elc)
   statusIon, dtSuggestedIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, distfIon, hamilIon, distf1Ion)
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)
   calcPhiFromChargeDensity(tCurr, myDt, distf1Elc, distf1Ion, phi1dDg)
   calcContinuousPhi(tCurr, myDt, phi1dDg, phi1d)
   calcHamiltonianElc(tCurr, myDt, phi1d, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, hamilIon)

   -- RK stage 2
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distf1Elc, hamilElc, distfNewElc)
   statusIon, dtSuggestedIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, distf1Ion, hamilIon, distfNewIon)
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   distf1Elc:combine(3.0/4.0, distfElc, 1.0/4.0, distfNewElc)
   distf1Ion:combine(3.0/4.0, distfIon, 1.0/4.0, distfNewIon)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)
   calcPhiFromChargeDensity(tCurr, myDt, distf1Elc, distf1Ion, phi1dDg)
   calcContinuousPhi(tCurr, myDt, phi1dDg, phi1d)
   calcHamiltonianElc(tCurr, myDt, phi1d, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, hamilIon)

   -- RK stage 3
   statusElc, dtSuggestedElc = updateVlasovEqn(vlasovSolverElc, tCurr, myDt, distf1Elc, hamilElc, distfNewElc)
   statusIon, dtSuggestedIon = updateVlasovEqn(vlasovSolverIon, tCurr, myDt, distf1Ion, hamilIon, distfNewIon)
   if (statusElc == false) or (statusIon == false) then
      return false, math.min(dtSuggestedElc, dtSuggestedIon)
   end
   distf1Elc:combine(1.0/3.0, distfElc, 2.0/3.0, distfNewElc)
   distf1Ion:combine(1.0/3.0, distfIon, 2.0/3.0, distfNewIon)
   applyBc(tCurr, myDt, distf1Elc, distf1Ion)

   distfElc:copy(distf1Elc)
   distfIon:copy(distf1Ion)

   calcPhiFromChargeDensity(tCurr, myDt, distfElc, distfIon, phi1dDg)
   calcContinuousPhi(tCurr, myDt, phi1dDg, phi1d)
   calcHamiltonianElc(tCurr, myDt, phi1d, hamilElc)
   calcHamiltonianIon(tCurr, myDt, phi1d, hamilIon)

   return true, math.min(dtSuggestedElc, dtSuggestedIon)
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      distfDupElc:copy(distfElc)
      distfDupIon:copy(distfIon)

      if (tCurr+myDt > tEnd) then
	      myDt = tEnd-tCurr
      end

      Lucee.logInfo (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	      distfElc:copy(distfDupElc)
	      distfIon:copy(distfDupIon)
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
function writeFields(frameNum, tCurr)
  distfElc:write(string.format("distfElc_%d.h5", frameNum), tCurr)
  distfIon:write(string.format("distfIon_%d.h5", frameNum), tCurr)

  numDensityElc:write(string.format("numDensityElc_%d.h5", frameNum), tCurr)
  numDensityIon:write(string.format("numDensityIon_%d.h5", frameNum), tCurr)

  momentumElc:write(string.format("momentumElc_%d.h5", frameNum), tCurr)
  momentumIon:write(string.format("momentumIon_%d.h5", frameNum), tCurr)

  ptclEnergyElc:write(string.format("ptclEnergyElc_%d.h5", frameNum), tCurr)
  ptclEnergyIon:write(string.format("ptclEnergyIon_%d.h5", frameNum), tCurr)

  copyPotential(0.0, 0.0, phi1d, phi1dDg)
  phi1dDg:write(string.format("phi_%d.h5", frameNum), tCurr)

  fieldEnergy:write( string.format("fieldEnergy_%d.h5", frameNum) )
end

-- calculate initial potential and Hamiltonians
applyBc(0.0, 0.0, distfElc, distfIon)
calcPhiFromChargeDensity(0.0, 0.0, distfElc, distfIon, phi1dDg)
calcHamiltonianElc(0.0, 0.0, phi1d, hamilElc)
calcHamiltonianIon(0.0, 0.0, phi1d, hamilIon)
-- compute initial diagnostics
calcDiagnostics(0.0, 0.0)
-- write out initial conditions
writeFields(0, 0.0)
-- make a duplicate in case we need it
distfDupElc = distfElc:duplicate()
distfDupIon = distfIon:duplicate()

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
