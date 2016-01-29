using tssm

#include("tssm_schroedinger.jl")
include("time_stepper.jl")

const k = -1.0
xmin=-16.0
xmax=16.0
Nx=1024
m=Schroedinger1D(Nx,xmin, xmax,cubic_coupling=k)


psi=wave_function(m)
const a=[2.0, 2.0,  0.0,  0.0] #aplitude
const b=[1.0,-3.0,  1.0, -9.0] #speed
const c=[5.0,-5.0,  0.0, -10.0] #offset x

function soliton(x)
    (a[1]/cosh(a[1]*(x-c[1]))) * exp(-b[1]*x*1im)+
    (a[2]/cosh(a[2]*(x-c[2]))) * exp(-b[2]*x*1im)+
    (a[3]/cosh(a[3]*(x-c[3]))) * exp(-b[3]*x*1im)+
    (a[4]/cosh(a[4]*(x-c[4]))) * exp(-b[4]*x*1im)
end

set!(psi, soliton)

palindromic_scheme_34 = PalindromicScheme( 
          ( 0.268330095781759925,  0.919661523017399857, 
           -0.187991618799159782, -0.187991618799159782, 
            0.919661523017399857,  0.268330095781759925 ),
            3 )

palindromic_scheme_56 = PalindromicScheme(
          ( 0.201651044312324230,   0.578800656272664932, 
            0.562615975356569200,   0.273128836056524479, 
            0.253874038247554845,  -0.102733803148432142, 
           -0.835351693190370636,   0.068014946093165092, 
            0.068014946093165092,  -0.835351693190370636,
           -0.102733803148432142,   0.253874038247554845, 
            0.273128836056524479,   0.562615975356569200, 
            0.578800656272664932,   0.201651044312324230 ),
            
            5 )

  x = get_nodes(psi)
  to_real_space!(psi)
  u = get_data(psi, true)
  

  
  plotdata=abs(u).^2
  steps=[0.0]
  mytime=[0.0]
  told=0.0
  nsteps=0
  tend=10.0
  t0=0.0
  out=1
  tol=1e-8
  
  tic()
  for t in adaptive_time_stepper2(psi, t0, tend, 0.01, tol, palindromic_scheme_56, "AB")
       #readline(STDIN)
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
  using PyPlot
  figure(1)
  hold(false)
  if out==1
  subplot(2, 1, 1)
  pcolormesh(mytime,x,plotdata)
  xlabel("t")
  ylabel("x")
  #colorbar()
  axis([t0,tend,xmin,xmax])
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
  set!(psi, soliton)
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
  
