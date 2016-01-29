using tssm
include("time_stepper.jl")
nx = 1024
xmin = -100.0
xmax =  100.0

function harmonic_trap_1D (x)
0.5*x^2
end

function weak_periodic_1D(x)
1.4 * cos(11.46*x) 
end


  include("groundstate.jl")
  m=SchroedingerReal1D(nx,xmin,xmax,hbar=1.0,mass=1.0,potential=harmonic_trap_1D,
                     cubic_coupling=390.0,boundary_conditions=periodic)
  psi=wave_function(m)
  groundstate!(psi, extrapolation_order=2)
  x = get_nodes(psi)
  to_real_space!(psi)
  u = get_data(psi, true)
  using PyPlot
  figure(2)
  hold(false)
  plot(x, abs(u[1:end-2]).^2)
  title("groundstate")
  hold(true)
  save(psi,"groundstate.h5")

  m=Schroedinger1D(4*nx,4*xmin,4*xmax,hbar=1.0,mass=1.0,potential=weak_periodic_1D,
                     cubic_coupling=390.0,boundary_conditions=periodic)
  psi=wave_function(m)
  load!(psi,"groundstate.h5")
  
  palindromic_scheme_34 = PalindromicScheme( 
          ( 0.268330095781759925,  0.919661523017399857, 
           -0.187991618799159782, -0.187991618799159782, 
            0.919661523017399857,  0.268330095781759925 ),
            3 )

  x = get_nodes(psi)
  to_real_space!(psi)
  u = get_data(psi, true)
  
  using PyPlot
  
  plotdata=abs(u).^2
  steps=[0.0]
  mytime=[0.0]
  told=0.0
  nsteps=0
  tend=10.0
  t0=0.0
  out=1
  tol=1e-4
  
  tic()
  for t in adaptive_time_stepper(psi, t0, tend, 0.01, tol, palindromic_scheme_34, "AB")
  push!(steps,t-told)
  push!(mytime,t)
  told=t
  nsteps=nsteps+1
  
  if out==1
  to_real_space!(psi)
  plotdata=[plotdata abs(u).^2]
     @printf("t=%5.3f\n", t)
     #figure(2)
     #hold(false)
     #plot(x, abs(u).^2)
     #hold(true)
  end
  end
  toc()

  figure(1)
  hold(false)
  if out==1
  subplot(2, 1, 1)
  pcolormesh(mytime,x,plotdata)
  xlabel("t")
  ylabel("x")
  #colorbar()
  axis([t0,tend,4*xmin,4*xmax])
  title("abs(u)^2")

  subplot(2, 1, 2)
  end
  plot(mytime[1:end-2], steps[2:end-1])
  xlabel("t")
  ylabel("stepsize")
  savefig("step.png", bbox_inches="tight")
  #colorbar( )
  stepmin=minimum(steps[2:end])
  @printf("dtmin=%1.20f  nsteps=%d\n", stepmin,nsteps)
  #compute again with constant stepsize.
  load!(psi,"groundstate.h5")
  nsteps=0
  tic()
  for t in equidistant_time_stepper(psi, t0, tend, stepmin, palindromic_scheme_34.scheme, "AB")
  nsteps=nsteps+1
  if out==1
    @printf("t=%5.3f\n", t)
  end
  end
  toc()
  @printf("dtmin=%1.20f  nsteps=%d\n", stepmin,nsteps)

