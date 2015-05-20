-- Used for figures in SOL paper

-- polynomial order
polyOrder = 2

-- cfl number to use
cfl = 0.1

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
N_Z, N_VPARA, N_MU = 8, 32, 16
-- domain extents
Z_LOWER, Z_UPPER = -lParallel, lParallel
VPARA_UPPER = math.min(4, 2.5*math.sqrt(N_VPARA/4))*vtPed
VPARA_LOWER = -VPARA_UPPER
MU_LOWER = 0
MU_UPPER = math.min(8, 4*math.sqrt(N_MU/2))*ionMass*vtPed*vtPed/B0

-- parameters to control time-stepping
tStart = 0.0
tEnd = 350e-6
nFrames = 5

-- Determine number of global nodes per cell for use in creating CG
-- fields. Note that this looks a bit odd as this not the number of
-- *local* nodes but the number of nodes in each cell to give the
-- correct number of global nodes in fields.
if (polyOrder == 1) then
   numCgNodesPerCell = 1
   numCgNodesPerCell_1d = 1 -- for spatial basis
elseif (polyOrder == 2) then
   numCgNodesPerCell = 3
   numCgNodesPerCell_1d = 2 -- for spatial basis
end

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
-- unit density ion distribution function
distfUnit = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- Hamiltonian (DG, but should be continuous)
hamil = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- DG version of the kinetic energy part of the hamiltonian
hamilKeDg = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

-- updater to initialize hamiltonian
initHamilKE = Updater.EvalOnNodes3D {
   onGrid = grid,
   basis = basis,
   shareCommonNodes = false,
   evaluate = function (z,vPara,mu,t)
      return ionMass*vPara^2/2 + mu*B0
   end
}
initHamilKE:setOut( {hamilKeDg} )
-- initialize potential
initHamilKE:advance(0.0) -- time is irrelevant

-- Return initial ne(x) in 1/m^3
numDensityElc = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- number density for ions
numDensity = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
numDensityIon3dDg = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}

function initialElcDens(x)
  local backgroundDens = 0.7 + 0.3*(1 - math.abs(x/lParallel))
  if math.abs(x) < lSource/2 then
       backgroundDens = backgroundDens + 0.5*math.cos(math.pi*x/lSource)
     end
  backgroundDens = backgroundDens*1e19
  return backgroundDens
end

-- Return initial Ti(x) in eV
function initialIonTemp(x)
  local backgroundTemp = 100 + 45*(1 - math.abs(x/lParallel))
  if math.abs(x) < lSource/2 then
       backgroundTemp = backgroundTemp + 30*math.cos(math.pi*x/lSource)
     end
  return backgroundTemp
end

initialElcDensUpdater = Updater.ProjectOnNodalBasis1D {
   onGrid = grid_1d,
   basis = basis_1d,
   shareCommonNodes = false,
   evaluate = function(x,y,z,t)
     return initialElcDens(x)
	 end
}
runUpdater(initialElcDensUpdater, 0.0, 0.0, {}, {numDensityElc})

-- updater to initialize distribution function
initDistf = Updater.ProjectOnNodalBasis3D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(z,vPara,mu,t)
     local backgroundTemp = initialIonTemp(z)
     local nHat = 2
     -- takes into account average v_ti over a cosine
     local vTi = math.sqrt(2*math.pi/(math.pi-1)*backgroundTemp*eV/ionMass)

     if z > lSource/2 then
       return nHat*maxwellianRight(ionMass, vTi, vPara, mu)
     elseif z < -lSource/2 then
       return nHat*maxwellianLeft(ionMass, vTi, vPara, mu)
     else
       -- Must be between the source boundaries, use a linear combo.
       return ((lSource/2 + z)*nHat*maxwellianRight(ionMass, vTi, vPara, mu) + 
          (lSource/2 - z)*nHat*maxwellianLeft(ionMass, vTi, vPara, mu))/lSource
     end
   end
}
runUpdater(initDistf, 0.0, 0.0, {}, {distfUnit})

-- Calculate initial ion density field
initIonDensityCalc = Updater.SOLIonDensityInitialization {
  onGrid = grid_1d,
  basis = basis_1d,
  kPerpTimesRho = 0.2,
}
runUpdater(initIonDensityCalc, 0.0, 0.0, {numDensityElc}, {numDensity})

-- Copy numDensityIon from a 1d to a 3d field
copyTo3DIonDg = Updater.NodalCopy1DTo3DFieldUpdater {
  onGrid = grid,
  basis1d = basis_1d,
  basis3d = basis,
  polyOrder = polyOrder,
}
runUpdater(copyTo3DIonDg, 0.0, 0.0, {numDensity}, {numDensityIon3dDg})

-- Multiply with gaussian to get total distfIon
multiply3dFields = Updater.FieldArithmeticUpdater3D {
  onGrid = grid,
  basis = basis,
  evaluate = function(An,Bn,t)
    return An*Bn
  end,
}
runUpdater(multiply3dFields, 0.0, 0.0, {numDensityIon3dDg, distfUnit}, {distf})

-- particle source
particleSource = DataStruct.Field3D {
   onGrid = grid,
   numComponents = basis:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
particleSource:clear(0.0)
-- updater to fill out particle source
particleSourceUpdater = Updater.ProjectOnNodalBasis3D {
   onGrid = grid,
   -- basis functions to use
   basis = basis,
   -- are common nodes shared?
   shareCommonNodes = false, -- In DG, common nodes are not shared
   -- function to use for initialization
   evaluate = function(z,vPara,mu,t)
    if math.abs(z) < lSource/2 then
      if t < tELM then
        return Sn*math.cos(math.pi*z/lSource)*maxwellian(ionMass, math.sqrt(tPed*eV/ionMass), vPara, mu)
      else
        return Sn/9*math.cos(math.pi*z/lSource)*maxwellian(ionMass, math.sqrt(260*eV/ionMass), vPara, mu)
      end
    else
      return 0
    end
	 end
}
particleSourceUpdater:setOut( {particleSource} )
-- initialize
particleSourceUpdater:advance(0.0)
-- keeps track of whether to update particlesources again
postELM = false

-- define equation to solve
pbEqn = PoissonBracketEquation.SOL3D {
  speciesMass = ionMass,
}
pbSlvrNew = Updater.PoissonBracketSimple3D {
  onGrid = grid,
  basis = basis,
  cfl = cfl,
  equation = pbEqn,
  updateDirections = {0,1},
  zeroFluxDirections = {1,2},
  onlyIncrement = false,
  fluxType = "upwind",
}

-- drift velocity u_parallel(x)
driftUPara = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- drift velocity u_perp(x)
driftUPerp = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- thermal velocity squared
vThermSq = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

vThermSqPara = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

vThermSqPerp = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- clear out contents
vThermSq:clear(0.0)

-- stuff to integrate vThermSq
integrateVThermSq = Updater.IntegrateGeneralField1D {
  onGrid = grid_1d,
  basis = basis_1d,
}
integratedVThermSq = DataStruct.DynVector { numComponents = 1, }

mom0Source = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

mom2Source = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- Ion Temperature Profile
tempProfile = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- field to store continous potential in 1D
phi1d = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = 1*numCgNodesPerCell_1d,
   -- ghost cells
   ghost = {1, 1},
}

-- field to store discontinous potential in 1D
phi1dDiscont = DataStruct.Field1D {
   onGrid = grid_1d,
   -- numNodesPerCell is number of global nodes stored in each cell
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}

-- to compute number density
mom0Calc = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = grid,
   basis3d = basis,
   basis1d = basis_1d,
   moment = 0,
}

mom1Para = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- to compute first moment
mom1Calc = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = grid,
   basis3d = basis,
   basis1d = basis_1d,
   moment = 1,
   momentDirection = 1,
}

mom1Mu = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
tPerpProfile = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- to compute first moment in mu
mom1MuCalc = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = grid,
   basis3d = basis,
   basis1d = basis_1d,
   moment = 1,
   momentDirection = 2,
}

-- ptcl energy
mom2Para = DataStruct.Field1D {
   onGrid = grid_1d,
   numComponents = basis_1d:numNodes(),
   ghost = {1, 1},
}
-- updater to compute ptcl energy
mom2CalcPara = Updater.DistFuncMomentCalc1DFrom3D {
   onGrid = grid,
   basis3d = basis,
   basis1d = basis_1d,
   moment = 2,
   momentDirection = 1,
}

-- Dynvector to store 0-3rd moments at left and right edges
momentsAtEdges = DataStruct.DynVector { numComponents = 8, }

momentsAtEdgesCalc = Updater.MomentsAtEdges3DUpdater {
  onGrid = grid,
  basis = basis,
  scaleFactor = 2*math.pi*B0/ionMass,
}

-- dynvector for heat flux at edge
heatFluxAtEdge = DataStruct.DynVector { numComponents = 3, }

-- to compute total particle energy
heatFluxAtEdgeCalc = Updater.HeatFluxAtEdge3DUpdater {
   onGrid = grid_1d,
   basis = basis_1d,
   ionMass = ionMass,
}

-- to calculate total temperature
vFromMomentsCalc = Updater.VelocitiesFromMoments3DUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
  scaleFactor = 2*math.pi*B0/ionMass,
}

-- to calculate parallel temperature
vParaFromMomentsCalc = Updater.VelocitiesFromMomentsUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
}

-- to calculate parallel temperature
vPerpFromMomentsCalc = Updater.VelocitiesFromMomentsUpdater {
  onGrid = grid_1d,
  basis = basis_1d,
}

-- compute moments from distribution function
function calcMoments(curr, dt, distfIn)
   runUpdater(mom0Calc, curr, dt, {distfIn}, {numDensity})
   numDensity:scale(2*math.pi*B0/ionMass)
   runUpdater(mom1Calc, curr, dt, {distfIn}, {mom1Para})
   mom1Para:scale(2*math.pi*B0/ionMass)
   runUpdater(mom2CalcPara, curr, dt, {distfIn}, {mom2Para})
   mom2Para:scale(2*math.pi*B0/ionMass)

   -- should be <v_\perp^2>
   runUpdater(mom1MuCalc, curr, dt, {distfIn}, {mom1Mu})
   mom1Mu:scale(2*math.pi*B0/ionMass*2*B0/ionMass)
   -- need to divide by density
   tPerpProfile:copy(mom1Mu)
   tPerpProfile:scale(0.5*ionMass/eV)
end

-- updater to copy 1d CG field to 1d DG field
copyCToD1d = Updater.CopyContToDisCont1D {
   onGrid = grid_1d,
   basis = basis_1d,
}

-- updater to compute phi from number density
adiabaticPhiCalc = Updater.BoltzmannPhiUpdater {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
   -- required physical parameters
   electronMass = electronMass,
   ionMass = ionMass,
   elementaryCharge = eV,
}

-- updater to move phi to continuous field
phiToContCalc = Updater.ContFromDisCont1D {
   -- grid for updater
   onGrid = grid_1d,
   -- basis functions to use
   basis = basis_1d,
}

-- Outflow BCs
function getRepTbl(pOrder, val)
   if pOrder == 1 then
      return {val, val, val, val, val, val, val, val}
   elseif pOrder == 2 then
      return {val, val, val, val, val, val, val, val,
              val, val, val, val,
              val, val, val, val, val, val, val, val}
   end
end

function getCountTbl(pOrder, val)
   if pOrder == 1 then
      return {0, 1, 2, 3, 4, 5, 6, 7}
   elseif pOrder == 2 then
      return {0, 1, 2, 3, 4, 5, 6, 7,
              8, 9, 10, 11,
              12, 13, 14, 15, 16, 17, 18, 19}
   end
end

bcConst = BoundaryCondition.Const { 
   components = getCountTbl(polyOrder),
   values = getRepTbl(polyOrder, 0.0),
}

-- No inflow boundary conditions?!
function makeBcObjIon()
   local bcLower = Updater.Bc3D {
      onGrid = grid,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "lower",
   }
   local bcUpper = Updater.Bc3D {
      onGrid = grid,
      boundaryConditions = {bcConst},
      dir = 0,
      edge = "upper",
   }
   return bcLower, bcUpper
end

-- make objects to apply BCs
bcLowerIon, bcUpperIon = makeBcObjIon()

-- apply boundary conditions to 2-D field
function applyBc(curr, dt, fldIon)
   --for i,bc in ipairs({bcLowerIon, bcUpperIon}) do
   --   runUpdater(bc, curr, dt, {}, {fldIon})
   --end

   runUpdater(bcLowerIon, curr, dt, {}, {fldIon})
   runUpdater(bcUpperIon, curr, dt, {}, {fldIon})
   -- (zero flux bc's applied via poisson updater)
end

-- compute hamiltonian
function calcHamiltonian(curr, dt, distIn, hamilOut)
   calcMoments(curr, dt, distIn)
   runUpdater(vFromMomentsCalc, curr, dt, {numDensity, mom1Para, mom1Mu, mom2Para}, {driftUPerp, vThermSqPerp})
   -- to calculate parallel temperature (used for heat flux)
   runUpdater(vParaFromMomentsCalc, curr, dt, {numDensity, mom1Para, mom2Para}, {driftUPara, vThermSqPara})

   phi1dDiscont:clear(0.0)
   runUpdater(adiabaticPhiCalc, curr, dt, {numDensity, mom1Para, vThermSqPara}, {phi1dDiscont})

   -- Project discontinuous potential on a continuous basis
   phi1d:clear(0.0)
   runUpdater(phiToContCalc, curr, dt, {phi1dDiscont}, {phi1d})

   -- Copy the continuous phi back onto a dg field
   runUpdater(copyCToD1d, curr, dt, {phi1d}, {phi1dDiscont})

   hamilOut:clear(0.0)
   -- Accumulate potential energy contribution to hamiltonian
   runUpdater(copyTo3DIonDg, curr, dt, {phi1dDiscont}, {hamilOut})
   hamilOut:scale(ionCharge)
   -- Accumulate kinetic energy contribution to hamiltonian
   hamilOut:accumulate(1.0, hamilKeDg)
end

-- compute various diagnostics
function calcDiagnostics(curr, dt)
  -- compute average temperature
  runUpdater(integrateVThermSq, curr, dt, {vThermSqPara}, {integratedVThermSq})
  -- compute moments at edges
  runUpdater(momentsAtEdgesCalc, curr, dt, {distf, hamilKeDg}, {momentsAtEdges})
  runUpdater(heatFluxAtEdgeCalc, curr, dt, {phi1dDiscont, momentsAtEdges, integratedVThermSq,
    numDensity, tPerpProfile}, {heatFluxAtEdge})
   -- Compute temperature profile in eV
  tempProfile:copy(vThermSq)
  tempProfile:scale(ionMass/eV)
end

-- function to take a time-step using SSP-RK3 time-stepping scheme
function rk3(tCurr, myDt)
   local pbStatus, pbDtSuggested

   if (postELM == false) then
     if tCurr + myDt > tELM then
       postELM = true
       particleSourceUpdater:advance(tCurr + myDt)
       -- Recalculate moments of particle source
       runUpdater(mom0Calc, tCurr, myDt, {particleSource}, {mom0Source})
       mom0Source:scale(2*math.pi*B0/ionMass)
       runUpdater(mom2CalcPara, tCurr, myDt, {particleSource}, {mom2Source})
       mom2Source:scale(2*math.pi*B0/ionMass)
     end
   end

   -- RK stage 1
   pbStatus, pbDtSuggested = runUpdater(pbSlvrNew, tCurr, myDt, {distf, hamil}, {distf1})
   if (pbStatus == false) then
     return false, pbDtSuggested
   end  

   distf1:accumulate(myDt, particleSource)

   applyBc(tCurr, myDt, distf1)
   calcHamiltonian(tCurr, myDt, distf1, hamil)

   -- RK stage 2
   pbStatus, pbDtSuggested = runUpdater(pbSlvrNew, tCurr, myDt, {distf1, hamil}, {distfNew})
   if (pbStatus == false) then
     return false, pbDtSuggested
   end  

   distfNew:accumulate(myDt, particleSource)
   distf1:combine(3.0/4.0, distf, 1.0/4.0, distfNew)

   applyBc(tCurr, myDt, distf1)
   calcHamiltonian(tCurr, myDt, distf1, hamil)

   -- RK stage 3
   pbStatus, pbDtSuggested = runUpdater(pbSlvrNew, tCurr, myDt, {distf1, hamil}, {distfNew})
   if (pbStatus == false) then
     return false, pbDtSuggested
   end   
   
   distfNew:accumulate(myDt, particleSource)
   distf1:combine(1.0/3.0, distf, 2.0/3.0, distfNew)

   applyBc(tCurr, myDt, distf1)
   distf:copy(distf1)
   calcHamiltonian(tCurr, myDt, distf, hamil)

   return true, pbDtSuggested
end

-- function to advance solution from tStart to tEnd
function advanceFrame(tStart, tEnd, initDt)
   local step = 1
   local tCurr = tStart
   local myDt = initDt
   local status, dtSuggested

   while tCurr<=tEnd do
      distfDup:copy(distf)
      hamilDup:copy(hamil)

      if (tCurr+myDt > tEnd) then
	      myDt = tEnd-tCurr
      end

      print (string.format("Taking step %d at time %g with dt %g", step, tCurr, myDt))
      status, dtSuggested = rk3(tCurr, myDt)

      if (status == false) then
	      print (string.format("** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
	      distf:copy(distfDup)
	      hamil:copy(hamilDup)
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
   --distf:write( string.format("distfIon_%d.h5", frame) )
   numDensity:write( string.format("numDensity_%d.h5", frameNum), tCurr)
   phi1dDiscont:write( string.format("phi_%d.h5", frameNum), tCurr)
   --hamil:write( string.format("hamil_%d.h5", frame) )
   --totalPtcl:write( string.format("totalPtcl_%d.h5", frame) )
   --fieldEnergy:write( string.format("fieldEnergy_%d.h5", frame) )
   heatFluxAtEdge:write( string.format("heatFluxAtEdge_%d.h5", frameNum), tCurr)
   momentsAtEdges:write( string.format("momentsAtEdges_%d.h5", frameNum), tCurr)
   tempProfile:write( string.format("tempProfile_%d.h5", frameNum), tCurr)
   tPerpProfile:write( string.format("tPerpProfile_%d.h5", frameNum), tCurr)
   --tempProfile:write( string.format("tempProfile_%d.h5", frame) )
   --driftU:write( string.format("driftU_%d.h5", frame) )
   --mom1Para:write( string.format("mom1Ion_%d.h5", frame) , tCurr)
   --numDensInCell:write( string.format("numDensInCell_%d.h5", frame) )
end

-- Calculate first moment of particle source at t = 0
runUpdater(mom0Calc, 0.0, 0.0, {particleSource}, {mom0Source})
mom0Source:scale(2*math.pi*B0/ionMass)
-- Calculate second moment of particle source at t = 0
runUpdater(mom2CalcPara, 0.0, 0.0, {particleSource}, {mom2Source})
mom2Source:scale(2*math.pi*B0/ionMass)

applyBc(0.0, 0.0, distf)
-- calculate initial Hamiltonian
calcHamiltonian(0.0, 0.0, distf, hamil)
-- compute initial diagnostics
calcDiagnostics(0.0, 0.0)
-- write out initial conditions
writeFields(0, 0.0)
-- make a duplicate in case we need it
distfDup = distf:duplicate()
hamilDup = hamil:duplicate()

-- parameters to control time-stepping
dtSuggested = 0.1*tEnd -- initial time-step to use (will be adjusted)
tFrame = (tEnd-tStart)/nFrames
tCurr = tStart

-- try to calculate temperature
--calcMoments(0.0, 0.0, distfUnit)
--runUpdater(vFromMomentsCalc, 0.0, 0.0, {numDensity, mom1Para, mom1Mu, mom2Para}, {driftU, vThermSq})
--runUpdater(vParaFromMomentsCalc, 0.0, 0.0, {numDensity, mom1Para, mom2Para}, {driftU, vThermSqPara})
--vThermSq:scale(ionMass/eV)
--vThermSq:write( string.format("temp_%d.h5", 0), 0)
--numDensity:write( string.format("n_%d.h5", 0), 0)
--driftU:write( string.format("u_%d.h5", 0), 0)
--tPerpProfile:write( string.format("tPerp_%d.h5", 0), 0)

--calcDiagnostics(0.0, 0.0)

for frame = 1, nFrames do
   Lucee.logInfo (string.format("-- Advancing solution from %g to %g", tCurr, tCurr+tFrame))
   dtSuggested = advanceFrame(tCurr, tCurr+tFrame, dtSuggested)
   writeFields(frame, tCurr+tFrame)
   tCurr = tCurr+tFrame
   Lucee.logInfo ("")
end
