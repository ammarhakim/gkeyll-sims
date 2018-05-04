-- 4D GK ETG nonlinear calculation
--
-- Plasma ------------------------------------------------------------------------
local Plasma = require "App.PlasmaOnCartGrid"
local Constants = require "Lib.Constants"

-- physical parameters
eV = Constants.ELEMENTARY_CHARGE
qe = -eV
qi = eV
me = Constants.ELECTRON_MASS
mi = 2.014*Constants.PROTON_MASS -- (deuterium ions)
Te0 = 2072*eV 
Ti0 = 2072*eV 
B0 = 1.91   -- [T]
R0 = 1.313  -- [m]
a  = 0.4701 -- [m]
n0 = 4.992*10^(19) -- [1/m^3]
-- derived parameters
R       = R0 + 0.5*a
vte  	= math.sqrt(Te0/me)
c_s     = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
omega_ce = math.abs(qe*B0/me)
rho_s   = c_s/omega_ci
rho_e   = vte/omega_ce
deltaR  = 32*rho_e
L_T     = R/4
ky_min  = 2*math.pi/deltaR
-- velocity grid parameters
N_VPAR, N_MU = 4, 2
VPAR_UPPER = math.min(4, 2.5*math.sqrt(N_VPAR/4))*vte
VPAR_LOWER = -VPAR_UPPER
MU_LOWER = 0
MU_UPPER = math.min(16, 4*math.sqrt(N_MU/2))*me*vte*vte/B0

-- background magnetic field profile
function Bmag(x) 
   return B0*R/x
end
-- background electron temperature profile
function Te(x)
   return Te0*(1-(x-R)/L_T)
end

plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 2.5e-5, -- end time
   nFrame = 9, -- number of output frames
   lower = {R, -deltaR/2}, -- configuration space lower left
   upper = {R+deltaR, deltaR/2}, -- configuration space upper right
   cells = {8, 8}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 1.0,

   -- decomposition for configuration space
   decompCuts = {1, 2}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1,2}, -- periodic directions

   -- gyrokinetic electrons
   electron = Plasma.GkSpecies {
      charge = qe,
      mass = me,
      -- velocity space grid
      lower = {VPAR_LOWER, MU_LOWER},
      upper = {VPAR_UPPER, MU_UPPER},
      cells = {N_VPAR, N_MU},
      decompCuts = {1, 1},
      -- initial conditions
      initBackground = function (t, xn, self)
         local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
         return self:Maxwellian(xn, n0, Te(x))
      end,
      init = function (t, xn, self)
         local x, y, vpar, mu = xn[1], xn[2], xn[3], xn[4]
         local x0 = R+deltaR/2
         local sigma = deltaR/4
         local perturb = 1e-2*rho_e/L_T*math.cos(ky_min*y)*math.exp(-(x-x0)^2/(2*sigma^2))
         return self:Maxwellian(xn, n0*(1+perturb), Te(x))
      end,
      fluctuationBCs = true, -- only apply BCs to fluctuations
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkDens"}, 
   },

   -- adiabatic ions
   adiabaticIon = Plasma.GkSpecies {
      charge = qi,
      mass = mi,
      -- velocity space grid
      lower = {VPAR_LOWER*math.sqrt(me/mi), MU_LOWER},
      upper = {VPAR_UPPER*math.sqrt(me/mi), MU_UPPER},
      cells = {N_VPAR, N_MU},
      decompCuts = {1, 1},
      -- initial conditions
      init = function (t, xn, self)
         return self:Maxwellian(xn, n0, Ti0)
      end,
      evolve = false, -- evolve species?
      diagnosticMoments = {"GkDens"}, 
   },

   -- field solver
   field = Plasma.GkField {
      evolve = true, -- evolve fields?
      adiabatic = {response = "ion", charge = qi, dens = n0, temp = Ti0},
      polarizationWeight = mi*n0/B0^2,
      discontinuous = false
   },

   -- magnetic geometry 
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x = xn[1]
         return Bmag(x)
      end,

      -- bcurvY = 1/B*curl(bhat).grad(y)
      bcurvY = function (t, xn)
         local x = xn[1]
         return -1/(B0*R)
      end,

      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
