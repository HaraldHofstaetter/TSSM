#ifdef _QUADPRECISION_
module tssmq_splitting_schemes
    use tssmq_base
#else
module tssm_splitting_schemes
    use tssm_base
#endif    

!Source: Coefficients of various splitting methods
!        http://www.asc.tuwien.ac.at/~winfried/splitting/index.php

   
!   AB scheme, real: Lie-Trotter
!   s=1, p=1
!   LEM: 1.0
    real(prec), target :: coeffs_lie_trotter(2) =  (/ 1.0_prec, 1.0_prec /) 
    
!   AB scheme, real: Strang
!   s=2, p=2, symmetric
!   LEM: 0.6    
    real(prec), target :: coeffs_strang(3) =  (/ 0.5_prec, 1.0_prec, 0.5_prec /) 

!   AB scheme, real: best 2-stage 2nd order
!   s=2, p=2, optimized; palindromic
!   LEM: 0.1 
!   optimum is global (over the domain of real or complex coefficients)
    real(prec), target :: coeffs_4(4) = &
                (/ 1.0_prec-0.5_prec*sqrt(2.0/prec), 0.5_prec*sqrt(2.0/prec), &
                   0.5_prec*sqrt(2.0/prec),          1.0_prec-0.5_prec*sqrt(2.0/prec) /) 

!   AB scheme, real: AM 3-2
!   s=3, p=2, symmetric (Auzinger 2015)
!   best symmetric 3-stage 2nd order scheme
!   see also McLachlan (1995): similar coefficients
!   LEM: 0.05
    real(prec), target :: coeffs_5(5) = &
                (/ 0.193183327503783685_prec,  0.5_prec, &
                   0.613633344992432630_prec,  0.5_prec, &
                   0.193183327503783685_prec /)

!   AB scheme, real: Milne 2/2 rational
!   (I) s=2, p=2, 'compromise' (Auzinger 2015)
!   LEM: 0.2                  
    real(prec), target :: coeffs_6(4) = &
                (/ 3.0_prec/8.0_prec, 4.0_prec/5.0_prec, &
                   5.0_prec/8.0_prec, 1.0_prec/5.0_prec /)

!   Associated scheme (Milne 2/2 rational)
!   (II) s=2, p=2
!   LEM: 2.1
!   'Milne device': With κ=1/9, κ*((II)-(I)) is local error estimator for (I)
    real(prec), target :: coeffs_6A(4) = &
                (/ 3.0_prec/4.0_prec,  2.0_prec, &
                   1.0_prec/4.0_prec, -1.0_prec /)

!   AB scheme, real: Emb 3/2 AKS 
!   s=3, p=3, optimized; palindromic (Auzinger/Ketcheson 2014)
!   see also Suzuki (1992)
!   LEM: 0.25
    real(prec), target :: coeffs_7(6) = &
                (/ 0.268330095781759925_prec,  0.919661523017399857_prec, &
                  -0.187991618799159782_prec, -0.187991618799159782_prec, &
                   0.919661523017399857_prec,  0.268330095781759925_prec /)

!   Associated scheme (Emb 3/2 AKS)
!   s=2, p=2, optimized
!   LEM: 0.1
    real(prec), target :: coeffs_7A(4) = &
                (/ 0.268330095781759925_prec,  0.683368274569431595_prec, &
                   0.731669904218240075_prec,  0.316631725430568405_prec /)

!   AB scheme, real: Y 4-4
!   s=4, p=4; symmetric 3-fold Strang composition recombined (Yoshida 1990)
!   the only real s=p=4 scheme
!   LEM: 3.7
    real(prec), parameter :: omega1 = 1.0_prec/(2.0_prec-2.0_prec**(1.0_prec/3.0_prec))
    real(prec), parameter :: omega2 = 1.0_prec - 2.0_prec*omega1
    real(prec), target :: coeffs_8(7) = &
                (/ 0.5_prec*omega1,          omega1, &
                   0.5_prec*(omega1+omega2), omega2, &
                   0.5_prec*(omega1+omega2), omega1, &
                   0.5_prec*omega1 /)

!   AB scheme, real: Emb 4/3 AK p
!   s=5, p=4, palindromic optimized (Auzinger/Koch 2015)
!   LEM: 0.06
    real(prec), target :: coeffs_9(10) = &
                (/ 0.125962888700250514_prec,  0.333588446797901933_prec, &
                   0.751193431379145450_prec, -0.338296598434303506_prec, &
                   0.127551831557005609_prec,  0.127551831557005609_prec, &
                  -0.338296598434303506_prec,  0.751193431379145450_prec, &
                   0.333588446797901933_prec,  0.125962888700250514_prec /)

!   Associated scheme ( Emb 4/3 AK p)
!   s=5, p=3, optimized
!   stages 2 & 3 may be recombined
!   LEM: 0.04
    real(prec), target :: coeffs_9A(10) = &
                (/ 0.125962888700250514_prec,  0.333588446797901933_prec, &
                   0.751193431379145450_prec, -0.338296598434303506_prec, &
                   0.0_prec,                   0.261153550449697153_prec, &
                  -0.242703571757396124_prec,  0.596114052266110425_prec, &
                   0.365547251678000160_prec,  0.147440548920593995_prec /)
!   AB scheme, real: Emb 4/3 AK s
!   s=5, p=4, symmetric (Auzinger/Koch 2015)
!   best symmetric 5-stage 4th order scheme
!   LEM: 0.2
    real(prec), target :: coeffs_10(9) = &
                (/ 0.267171359000977615_prec,  -0.361837907604416033_prec, &
                  -0.0338279096695056672_prec,  0.861837907604416033_prec, &
                   0.5333131013370561044_prec,  0.861837907604416033_prec, &
                  -0.0338279096695056672_prec, -0.361837907604416033_prec, &
                   0.267171359000977615_prec /) 

!   Associated scheme (Emb 4/3 AK s)
!   s=5, p=3, optimized
!   LEM: 0.03
    real(prec), target :: coeffs_10A(10) = &
                (/ 0.267171359000977615_prec,  -0.361837907604416033_prec, &
                  -0.0338279096695056672_prec,  0.861837907604416033_prec, &
                   0.5333131013370561044_prec,  0.395088376480991403_prec, &
                   0.267171359000977615_prec,  -0.361837907604416033_prec, &
                  -0.0338279096695056672_prec,  0.466749531123424630_prec /)
!   AB scheme, real: Emb 4/3 M/AK
!   s=6, p=4, symmetric, optimized
!   see also McLachlan (1995): similar coefficients
!   LEM: 0.02
    real(prec), target :: coeffs_11(11) = &
                (/ 0.093500348726330576_prec,   0.439051727817158558_prec, &
                  -0.069094369881095037_prec,  -0.136536314071511211_prec, &
                   0.4755940211547644618_prec,  0.394969172508705306_prec, &
                   0.4755940211547644618_prec, -0.136536314071511211_prec, &
                  -0.069094369881095037_prec,   0.439051727817158558_prec, &
                   0.093500348726330576_prec /) 

!   Associated scheme (Emb 4/3 M/AK)
!   s=6, p=3, optimized (Auzinger/Koch 2015)
!   LEM: 0.01
    real(prec), target :: coeffs_11A(11) = &
                (/ 0.093500348726330576_prec,  0.439051727817158558_prec, &
                  -0.069094369881095037_prec, -0.136536314071511211_prec, &
                   0.4755940211547644618_prec, 0.394969172508705306_prec, &
                   0.3441299989650415828_prec, 0.093500348726330576_prec, &
                  -0.069094369881095037_prec, -0.136536314071511211_prec, &
                   0.224964370916053455_prec /)

!   AB scheme, real: Emb 4/3 BM PRK/A
!   s=7, p=4, symmetric, 'PRK', optimized (Blanes/Moan 2002)
!   LEM: 0.01           
    real(prec), target :: coeffs_12(13) = &
                (/ 0.0792036964311954608_prec,  0.209515106613361891_prec, &
                   0.353172906049773948_prec,  -0.143851773179818077_prec, &
                  -0.0420650803577191948_prec,  0.434336666566456186_prec, &
                   0.2193769557534995720_prec,  0.434336666566456186_prec, &
                  -0.0420650803577191948_prec, -0.143851773179818077_prec, &
                   0.353172906049773948_prec,   0.209515106613361891_prec, &
                   0.0792036964311954608_prec /)

!    Associated scheme (Emb 4/3 BM PRK/A)
!    s=4, p=3, optimized (Auzinger 2015)
!    LEM: 0.1
    real(prec), target :: coeffs_12A(8) = &
                (/ 0.0792036964311954608_prec, 0.209515106613361891_prec, &
                   0.353172906049773948_prec,  0.634279607840446390_prec, &
                  -0.351497633364616918_prec, -0.0576123055740887430_prec, &
                   0.919121030883647509_prec,  0.213817591120280462_prec /)

!   AB scheme, real: PP 3/4 A
!   s=3, p=3, optimized; palindromic (Auzinger/Ketcheson 2014)
!   same as in Emb 3/2 AKS
!   LEM: 0.15
!   Let (I) = given scheme, and (II) = (I) applied with A,B exchanged
!   Then the additive scheme ((I)+(II))/2 has order p=4; LEM: 0.4
!   ((I)-(II))/2 is local error estimator for (I)
    real(prec), target :: coeffs_13(6) = &
                (/ 0.268330095781759925_prec,  0.919661523017399857_prec, &
                  -0.187991618799159782_prec, -0.187991618799159782_prec, &
                   0.919661523017399857_prec,  0.268330095781759925_prec /)
!   AB scheme, real: Emb 5/4 AK (i)
!   s=7, p=5 (Auzinger/Ketcheson 2013,2015)
!   LEM: 0.3
    real(prec), target :: coeffs_14(14) = &
                (/ 0.475018345144539497_prec, -0.402020995028838599_prec, &
                   0.021856594741098449_prec,  0.345821780864741783_prec, &
                  -0.334948298035883491_prec,  0.400962967485371350_prec, &
                   0.512638174652696736_prec,  0.980926531879316517_prec, &
                  -0.011978701020553904_prec, -1.362064898669775624_prec, &
                  -0.032120004263046859_prec,  0.923805029000837468_prec, &
                   0.369533888781149572_prec,  0.112569584468347105_prec /)

!   Associated scheme (Emb 5/4 AK (i))
!   s=7, p=4
!   LEM: 0.1
    real(prec), target :: coeffs_14A(13) = &
                (/ 0.475018345144539497_prec, -0.402020995028838599_prec, &
                   0.021856594741098449_prec,  0.345821780864741783_prec, &
                  -0.334948298035883491_prec,  0.392182083210588712_prec, &
                   0.541870352983573775_prec, -0.128645858055930737_prec, &
                  -0.107236018614181158_prec,  0.554422771123937193_prec, &
                   0.307677712342444599_prec,  0.238240217885501648_prec, &
                   0.095761311438408329_prec /)

!   AB scheme, real: Emb 5/4 AK (ii)
!   s=7, p=5 (Auzinger/Ketcheson 2013,2015)
!   LEM: 0.3
    real(prec), target :: coeffs_15(14) = &
                (/ 0.475018345144539497_prec, -0.402020995028838599_prec, &
                   0.021856594741098449_prec,  0.345821780864741783_prec, &
                  -0.334948298035883491_prec,  0.400962967485371350_prec, &
                   0.512638174652696736_prec,  0.980926531879316517_prec, &
                  -0.011978701020553904_prec, -1.362064898669775624_prec, &
                  -0.032120004263046859_prec,  0.923805029000837468_prec, &
                   0.369533888781149572_prec,  0.112569584468347105_prec /)

!   Associated scheme (Emb 5/4 AK (ii))
!   s=8, p=4, optimized
!   LEM: 0.02
    real(prec), target :: coeffs_15A(16) = &
                (/ 0.475018345144539497_prec,   -0.402020995028838599_prec, &
                   0.021856594741098449_prec,    0.345821780864741783_prec, &
                  -0.334948298035883491_prec,    0.400962967485371350_prec, &
                   0.512638174652696736_prec,    0.980926531879316517_prec, &
                  -0.00596874002994121298_prec, -1.28006227193004091_prec, &
                  -0.0418443254169191836_prec,   0.709119365029440952_prec, &
                   0.0616955797702736790_prec,   0.132683037231661766_prec, &
                   0.31155266917413552658_prec,  0.112569584468347105_prec /)

!   AB scheme, real: Emb 5/4 A
!   s=8, p=5, best palindromic found (Auzinger 2015)
!   LEM: 0.15
    real(prec), target :: coeffs_16(16) = &
                (/ 0.201651044312324230_prec,   0.578800656272664932_prec, &
                   0.562615975356569200_prec,   0.273128836056524479_prec, &
                   0.253874038247554845_prec,  -0.102733803148432142_prec, &
                  -0.835351693190370636_prec,   0.068014946093165092_prec, &
                   0.068014946093165092_prec,  -0.835351693190370636_prec, & 
                  -0.102733803148432142_prec,   0.253874038247554845_prec, &
                   0.273128836056524479_prec,   0.562615975356569200_prec, &
                   0.578800656272664932_prec,   0.201651044312324230_prec /)

!   Associated scheme (Emb 5/4 A)
!   s=7, p=4, optimized
!   LEM: 0.04  
    real(prec), target :: coeffs_16A(14) = &
                (/ 0.201651044312324230_prec,   0.578800656272664932_prec, &
                   0.562615975356569200_prec,   0.273128836056524479_prec, &
                   0.253874038247554845_prec,  -0.102733803148432142_prec, &
                  -0.726560425753743689_prec,  -0.800888361904193390_prec, &
                  -0.0944858667147550235_prec,  0.259957296127666422_prec, &
                   0.246962230324312930_prec,   0.594143062900597984_prec, &
                   0.5559430042277375075_prec,  0.197592313695171715_prec /)

!   AB scheme, real: PP 5/6 A
!   s=8, p=5, best palindromic found (Auzinger 2015)
!   same as in Emb 5/4 A
!   LEM: 0.15
!   Let (I) = given scheme, and (II) = (I) applied with A,B exchanged
!   Then the additive scheme ((I)+(II))/2 has order p=6; LEM: 0.7
!   ((I)-(II))/2 is local error estimator for (I)
    real(prec), target :: coeffs_17(16) = &
                (/ 0.201651044312324230_prec,   0.578800656272664932_prec, &
                   0.562615975356569200_prec,   0.273128836056524479_prec, &
                   0.253874038247554845_prec,  -0.102733803148432142_prec, &
                  -0.835351693190370636_prec,   0.068014946093165092_prec, &
                   0.068014946093165092_prec,  -0.835351693190370636_prec, & 
                  -0.102733803148432142_prec,   0.253874038247554845_prec, &
                   0.273128836056524479_prec,   0.562615975356569200_prec, &
                   0.578800656272664932_prec,   0.201651044312324230_prec /)

!   AB scheme, real: Y 8-6
!   s=8, p=6; 7-fold symmetric Strang composition recombined (Yoshida 1990)
!   LEM: 16.
    real(prec), parameter :: w1 = -1.177679984178871007
    real(prec), parameter :: w2 =  0.2355732133593581337
    real(prec), parameter :: w3 = 0.7845136104775572638
    real(prec), parameter :: w0 = 1.0_prec-2.0_prec*(w1+w2+w3)
    real(prec), target :: coeffs_18(15) = &
                (/ 0.5_prec*w3,      w3, &
                   0.5_prec*(w2+w3), w2, &
                   0.5_prec*(w1+w2), w1, &
                   0.5_prec*(w0+w1), w0, &
                   0.5_prec*(w0+w1), w1, &
                   0.5_prec*(w1+w2), w2, &
                   0.5_prec*(w2+w3), w3, &
                   0.5_prec*w3 /)

!   AB scheme, real: BM 11-6 PRK
!   s=11, p=6, symmetric, 'PRK', optimized (Blanes/Moan 2002)
!   LEM: 0.05                  
    real(prec), target :: coeffs_19(21) = &
                (/ 0.0502627644003922_prec,  0.148816447901042_prec, &
                   0.413514300428344_prec,  -0.132385865767784_prec, &
                   0.0450798897943977_prec,  0.067307604692185_prec, &
                  -0.188054853819569_prec,   0.432666402578175_prec, &
                   0.541960678450780_prec,  -0.016404589403618_prec, &
                  -0.7255255585086898_prec, -0.016404589403618_prec, &
                   0.541960678450780_prec,   0.432666402578175_prec, &
                  -0.188054853819569_prec,   0.067307604692185_prec, &
                   0.0450798897943977_prec, -0.132385865767784_prec, &
                   0.413514300428344_prec,   0.148816447901042_prec, &
                   0.0502627644003922_prec /)

!   ABC scheme, real: PP 3/4 A 3
!   s=6, p=3, best palindromic found (Auzinger 2015)
!   LEM: 0.7
!   Let (I) = given scheme, and (II) = (I) applied with A,B,C replaced by C,B,A
!   Then the additive scheme ((I)+(II))/2 has order p=4
!   ((I)-(II))/2 is local error estimator for (I)
    real(prec), target :: coeffs_ABC_6(18) = &

                (/ 0.461601939364879971_prec, -0.266589223588183997_prec,  -0.360420727960349671_prec, &
                 -0.0678710530507800801_prec, 0.0924576733143338354_prec,   0.579154058410941403_prec, &
                 -0.0958868852260720250_prec,  0.674131550273850162_prec,   0.483422668461380403_prec, &
                   0.483422668461380403_prec,  0.674131550273850162_prec, -0.0958868852260720250_prec, &
                   0.579154058410941403_prec, 0.0924576733143338354_prec, -0.0678710530507800801_prec, &
                  -0.360420727960349671_prec, -0.266589223588183997_prec,   0.461601939364879971_prec /)

!!!!!! COMPLEX SCHEMES !!!!!!!!!!!!!!

!   AB scheme, complex: Milne 2/2 c (i)
!   (I) s=2, p=2 (Auzinger 2015)
!   LEM: 0.15
    complex(prec), target :: coeffs_c1(4) = &
       (/ cmplx(12.0_prec/37.0_prec, -2.0_prec/37.0_prec, prec), &
          cmplx(25.0_prec/34.0_prec, -1.0_prec/17.0_prec, prec), &
          cmplx(25.0_prec/37.0_prec, +2.0_prec/37.0_prec, prec), &
          cmplx(9.0_prec/34.0_prec,  +1.0_prec/17.0_prec, prec) /)

!   Associated scheme (Milne 2/2 c (i))
!   (II) s=2, p=2
!   LEM: 2.
!   'Milne device': With κ=14/15 - i/20, κ*((I)-(II)) is local error estimator for (I)          
    complex(prec), target :: coeffs_c1A(4) = &
       (/ cmplx(4.0_prec/5.0_prec, 2.0_prec/5.0_prec, prec), &
          cmplx(1.0_prec/2.0_prec, 1.0_prec, prec), &
          cmplx(1.0_prec/5.0_prec, -2.0_prec/5.0_prec, prec), &
          cmplx(1.0_prec/2.0_prec,-1.0_prec, prec) /)

!   AB scheme, complex: Milne 2/2 c (ii)
!   (I) s=2, p=2, optimized (Auzinger 2015)
!   LEM: 0.1      
    complex(prec), target :: coeffs_c2(4) = &
       (/ cmplx(4.0_prec/13.0_prec,  -1.0_prec/26.0_prec, prec), &
          cmplx(18.0_prec/25.0_prec, -1.0_prec/25.0_prec, prec), &
          cmplx(9.0_prec/13.0_prec,  +1.0_prec/26.0_prec, prec), &
          cmplx(7.0_prec/25.0_prec,  +1.0_prec/25.0_prec, prec) /)

!   Associated scheme(Milne 2/2 c (ii))
!   (II) s=2, p=2, a2 purely imaginary
!   LEM: 2.
!   'Milne device': With κ=1/20 + i/60, κ*((I)-(II)) is local error estimator for (I)          
    complex(prec), target :: coeffs_c2A(4) = &
       (/ cmplx(1.0_prec, +1.0_prec/2.0_prec, prec), &
          cmplx(0.0_prec, 1.0_prec, prec), &
          cmplx(0.0_prec, -1.0_prec/2.0_prec, prec), &
          cmplx(1.0_prec,  -1.0_prec, prec) /)

!   AB scheme, complex: Emb 3/2 A c
!   s=3, p=3, optimal?; palindromic (Auzinger 2015)
!   LEM: 0.3      
    complex(prec), target :: coeffs_c4(6) = &
     (/ cmplx(0.201639688260407656_prec, + 0.105972321241365172_prec, prec), & 
        cmplx(0.387747410753696807_prec, + 0.100071120693574555_prec, prec), &
        cmplx(0.410612900985895537_prec, - 0.206043441934939727_prec, prec), &
        cmplx(0.410612900985895537_prec, - 0.206043441934939727_prec, prec), &
        cmplx(0.387747410753696807_prec, + 0.100071120693574555_prec, prec), &
        cmplx(0.201639688260407656_prec, + 0.105972321241365172_prec, prec) /) 

!   Associated scheme(Emb 3/2 A c)
!   s=3, p=2
!   LEM: 0.15
    complex(prec), target :: coeffs_c4A(5) = &
     (/ cmplx(0.201639688260407656_prec, + 0.105972321241365172_prec, prec), & 
        cmplx(0.387747410753696807_prec, + 0.100071120693574555_prec, prec), &
        cmplx(0.502190133760050497_prec, - 0.0910042566309608042_prec, prec), &
        cmplx(0.612252589246303193_prec, - 0.100071120693574555_prec, prec), &
        cmplx(0.296170177979541847_prec, - 0.0149680646104043678_prec, prec) /)




#ifdef _QUADPRECISION_
end module tssmq_splitting_schemes
#else
end module tssm_splitting_schemes
#endif
