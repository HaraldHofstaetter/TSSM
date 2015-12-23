using tssm

include("coupled_schroedinger.jl")
include("time_stepper.jl")


const alpha=1.0
const e=0.8
const v=1.1
const delta=0.5

m=CoupledSchroedinger1D(512,-20, 60, delta, e)
psi=wave_function(m)

function ex_sol1(x, t)
    sqrt(2*alpha/(1+e)) * sech(sqrt(2*alpha)*(x-v*t)) * exp(1im*( (v-delta)*x - ( 0.5*(v^2-delta^2)-alpha )*t ))
end

function ex_sol2(x, t)
    sqrt(2*alpha/(1+e)) * sech(sqrt(2*alpha)*(x-v*t)) * exp(1im*( (v+delta)*x - ( 0.5*(v^2-delta^2)-alpha )*t ))
end

ex_sol = Function2(ex_sol1, ex_sol2)

set!(psi, ex_sol, 0.0)

embedded_scheme = EmbeddedScheme(
          ( 0.475018345144539497,   -0.402020995028838599,
            0.021856594741098449,    0.345821780864741783,
           -0.334948298035883491,    0.400962967485371350,
            0.512638174652696736,    0.980926531879316517,
           -0.00596874002994121298, -1.28006227193004091,
           -0.0418443254169191836,   0.709119365029440952,
            0.0616955797702736790,   0.132683037231661766,
            0.31155266917413552658,  0.112569584468347105 ),

          ( 0.475018345144539497, -0.402020995028838599,
            0.021856594741098449,  0.345821780864741783,
           -0.334948298035883491,  0.400962967485371350,
            0.512638174652696736,  0.980926531879316517,
           -0.011978701020553904, -1.362064898669775624,
           -0.032120004263046859,  0.923805029000837468,
            0.369533888781149572,  0.112569584468347105 ),

            4 )

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

palindromic_scheme_34 = PalindromicScheme( 
          ( 0.268330095781759925,  0.919661523017399857, 
           -0.187991618799159782, -0.187991618799159782, 
            0.919661523017399857,  0.268330095781759925 ),
            3 )


scheme = embedded_scheme.scheme1

#local_orders(psi, ex_sol, 0.0, 1.0, embedded_scheme.scheme1, "AB", 12 )
#local_orders(psi, ex_sol, 0.0, 1.0, embedded_scheme.scheme2, "AB", 12 )

#local_orders(psi, ex_sol, 0.0, 0.1, palindromic_scheme_34, "AB", 10 );
#local_orders(psi, ex_sol, 0.0, 0.1, palindromic_scheme_56, "AB", 10 );

#psi_ex=wave_function(m)
#set!(psi_ex, ex_sol, 5.0)
#global_orders(psi, psi_ex, 0.0, 5.0, 0.1, palindromic_scheme_34.scheme, "AB", 8 )
#global_orders(psi, psi_ex, 0.0, 5.0, 0.1, palindromic_scheme_56.scheme, "AB", 8 )


 
  x = get_nodes(psi.psi1)
  u = get_data(psi.psi1, true)
  
  using PyPlot
  
  hold(false)
  
  psi_ex=wave_function(m)
  u_ex = get_data(psi_ex.psi1, true)
  
  
  for t in adaptive_time_stepper(psi, 0.0, 100.0, 0.01, 1e-4, palindromic_scheme_34, "AB")
  #for t in equidistant_time_stepper(psi, 0.0, 100.0, 0.1, palindromic_scheme_34.scheme, "AB")
     @printf("t=%5.3f\n", t)
     to_real_space!(psi.psi1)
     plot(x, abs(u).^2)
     hold(true)
     set!(psi_ex, ex_sol, t)
     plot(x, abs(u_ex).^2)
     hold(false)
     axis([-20,60,0,1.2])
  #   #readline(STDIN)
  end
  
