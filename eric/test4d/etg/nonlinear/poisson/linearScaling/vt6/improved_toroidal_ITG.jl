#= Plot ITG/ETG dispersion function, using Julia.

execute in Julia by typing 

    include("toroidal_ITG.jl")

This finds approximate roots of the toroidal dispersional relation
D(omega) = 0, from Eq. 2.63, 2.67-2.69 of Beer's thesis (and Biglari et
al. 1989).  omega here is in units of

     omega_d = (v_t/R) k_y rho

This doesn't use any rigorous root finders, it just plots
|D(\omega)| over a grid in the complex plane, and prints out local
minima if fintds.

Note, because of the unusual branch cut issues for the toroidal
dispersion relation, the only physically meaningful roots are for
unstable modes.  (Damped dynamics are dominated by the branch cut, so
one gets algebraic decay instead of exponential decay.)

GWH 2014.10.24

Modified to include polarization charge 8.26.2015
Modified to do an iterative search for a global minimum assumed
to exist in starting interval 9.16.2015

=#

# -----------------------------------------------------------------
# Main inputs:

R_Ln = 0; R_LT = 20.235226; Ts_ov_Ta = 2047.705670/2072; kPerpRhoS = 2*pi/32

# Search for roots in the range of frequencies (real and imaginary parts
#  of omega/omega_d):
omega_min=4; omega_max=6
gamma_min=3; gamma_max=6

# Number of grid points to plot in each direction:
N_plot = 100+1

# end of inputs
# -----------------------------------------------------------------

# Plasma dispersion function Z(zeta), expressed in terms of the
# complementary error function with a complex argument:

Z(zeta) = im*sqrt(pi)*w(zeta)

# w(zeta) is sometimes called the Voight function or Faddeeva's function, and
# is related to the complementary error function by:

w(zeta) = exp(-zeta^2) * erfc(-im*zeta)

# Need the complex() call in order for sqrt to work with a negative real 
# argument:
R_0(zeta) = 1 - (zeta/2) *( Z(sqrt2(complex(zeta)/2)) )^2

R_1(zeta) = (1/2.) *( Z(sqrt2(complex(zeta)/2)) )^2

R_2(zeta) = (zeta/2 - 1/2)*( Z(sqrt2(complex(zeta)/2)) )^2 +
            sqrt2(complex(zeta)/2)*Z(sqrt2(complex(zeta)/2))

# try redefining the branch cuts:
# (theta_branch is angle by which to rotate the branch cut)
theta_branch=0 # traditional branch cut (on -real axis)
# theta_branch=pi/2 # along the -imaginary axis
# theta_branch=pi*(1-2*eps()) # just below +real axis

sqrt2(z) = exp(im*theta_branch/2)*sqrt(exp(-im*theta_branch)*z)

#=

Physical results can't depend on the choice of a branch cut.
Instabilities are uniquely determined, independent of the choice of
the branch cut (in the lower half plane).  But for damped modes, one
has to do a Laplace transform, and one will pick up a dominant
contribution from deforming around the branch cut.  Might be simplest
to move the branch cut to just below the positive real axis in that
case, because then the poles all disappear?  But don't bother worrying
about damped responses right now.

=#

x=linspace(omega_min, omega_max, N_plot)
y=linspace(gamma_min, gamma_max, N_plot)

z = x .+ y'*im

# Define dispersion relation D(omega) = 0,
# from Eq. 2.63, 2.67-2.69 of Beer's thesis:
disp(z) = R_0(z) + R_Ln*R_1(z) + R_LT*R_2(z) + Ts_ov_Ta*(1 + 1/Ts_ov_Ta*kPerpRhoS^2)

# disp takes 1 complex number as an argument.  Could define a function
# that can apply this to a whole array of complex numbers elementwise,
# but debug that another time:
#
# disp = broadcast_function(disp1)

derr = Array(Float64, N_plot, N_plot)

# broadcast applies the function disp elementwise to z:
# derr = abs(broadcast(disp,z))

targetTol = 1e-8
err = 1
iterations = 0

while err > targetTol
  # Evaluate dispersion relation at every point on grid
  for j = 1:N_plot
      for i = 1:N_plot
         derr[i,j] = abs(disp(z[i,j]))
      end
  end

  val,index = findmin(derr)
  i,j = ind2sub(size(derr), index)

  err = sqrt((x[i+1]-x[i])^2 + (y[j+1]-y[j])^2)/2
  
  @printf("iteration %d: local min = %10.8f at omega_r, gamma = %10.8f % 10.8f +/-%10.8f\n",
            iterations,derr[i,j], x[i], y[j], err)
  iterations = iterations + 1

  # refine grid
  x=linspace(x[i-1], x[i+1], N_plot)
  y=linspace(y[j-1], y[j+1], N_plot)
  z = x .+ y'*im
end


