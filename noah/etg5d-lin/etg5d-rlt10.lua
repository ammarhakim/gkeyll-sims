-- 5D GK ETG linear instability calculation with kz=0. results match corresponding 4D ETG case.
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
q = 8
-- derived parameters
R       = R0 + 0.5*a
L_parallel = 2*math.pi*R*q
vte  	= math.sqrt(Te0/me)
c_s     = math.sqrt(Te0/mi)
omega_ci = math.abs(qi*B0/mi)
omega_ce = math.abs(qe*B0/me)
rho_s   = c_s/omega_ci
rho_e   = vte/omega_ce
deltaR  = 32*rho_e
L_T     = R/10
ky_min  = 2*math.pi/deltaR
kz_min  = 2*math.pi/L_parallel
-- velocity grid parameters
N_VPAR, N_MU = 16, 8
VPAR_UPPER = math.min(4, 2.5*math.sqrt(N_VPAR/4))*vte
VPAR_LOWER = -VPAR_UPPER
MU_LOWER = 0
MU_UPPER = math.min(16, 4*math.sqrt(N_MU/2))*me*vte*vte/B0


plasmaApp = Plasma.App {
   logToFile = true,

   tEnd = 1e-6, -- end time
   nFrame = 9, -- number of output frames
   lower = {R, -deltaR/2, 0}, -- configuration space lower left
   upper = {R+deltaR, deltaR/2, L_parallel}, -- configuration space upper right
   cells = {1, 8, 8}, -- configuration space cells
   basis = "serendipity", -- one of "serendipity" or "maximal-order"
   polyOrder = 1, -- polynomial order
   timeStepper = "rk3", -- one of "rk2" or "rk3"
   cflFrac = 1.0,

   -- decomposition for configuration space
   decompCuts = {1, 2, 1}, -- cuts in each configuration direction
   useShared = false, -- if to use shared memory

   -- boundary conditions for configuration space
   periodicDirs = {1,2,3}, -- periodic directions

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
      initBackground = {"maxwellian", 
              density = function (t, xn)
                 return n0
              end,
              temperature = function (t, xn)
                 local x = xn[1]
                 return Te0*(1-(x-R)/L_T)
              end,
             },
      init = {"maxwellian", 
              density = function (t, xn)
                 local x, y, z = xn[1], xn[2], xn[3]
                 local perturb = 1e-3*rho_e/L_T*math.cos(ky_min*y + kz_min*z)
                 return n0*(1+perturb)
              end,
              temperature = function (t, xn)
                 local x = xn[1]
                 return Te0*(1-(x-R)/L_T)
              end,
             },
      fluctuationBCs = true, -- only apply BCs to fluctuations
      evolve = true, -- evolve species?
      diagnosticMoments = {"GkDens"}, 
   },

   -- adiabatic ions
   adiabaticIon = Plasma.AdiabaticSpecies {
      charge = qi,
      mass = mi,
      temp = Ti0,
      -- initial conditions
      init = function (t, xn)
         return n0
      end,
      evolve = false, -- evolve species?
   },

   -- field solver
   field = Plasma.GkField {
      evolve = true, -- evolve fields?
   },

   -- magnetic geometry 
   funcField = Plasma.GkGeometry {
      -- background magnetic field
      bmag = function (t, xn)
         local x = xn[1]
         return B0*R/x
      end,
      -- geometry is not time-dependent
      evolve = false,
   },
}
-- run application
plasmaApp:run()
