#ifdef _QUADPRECISION_

#include <stdint.h>
#include <quadmath.h> 

#if 1
typedef union
{
    __float128 value;

    struct{
        uint64_t u0;
        uint64_t u1;
    } words64;

} myfloat128;

typedef union
{
    __complex128 value;

    struct{
        uint64_t u0;
        uint64_t u1;
        uint64_t u2;
        uint64_t u3;
    } words64;

} mycomplex128;


#define F(x) (x.value)

#else

typedef __float128  myfloat128;
typedef __complex128  mycomplex128;

#define F(x) (x)

#endif

#ifdef _REAL_
  #define _COMPLEX_OR_REAL_ __float128
  #define _WRAPPED_COMPLEX_OR_REAL_ myfloat128
#else
  #define _COMPLEX_OR_REAL_ __complex128
  #define _WRAPPED_COMPLEX_OR_REAL_ mycomplex128
#endif

#ifdef _ROTSYM_
  #define _DIM_ 1
  #ifdef _REAL_
    #ifdef _QUADPRECISION_
     #define S(x) tssmq_ ## x ## _bessel_rotsym_real_1d
     #define W(x) tssmq_ ## x ## _bessel_rotsym_real_1d_wf128
    #else
     #define S(x) tssm_ ## x ## _bessel_rotsym_real_1d
     #define W(x) tssm_ ## x ## _bessel_rotsym_real_1d_wf128
    #endif 
  #else
    #ifdef _QUADPRECISION_
     #define S(x) tssmq_ ## x ## _bessel_rotsym_1d
     #define W(x) tssmq_ ## x ## _bessel_rotsym_1d_wf128
    #else
     #define S(x) tssm_ ## x ## _bessel_rotsym_1d
     #define W(x) tssm_ ## x ## _bessel_rotsym_1d_wf128
    #endif
  #endif
#else
  #define _DIM_ 2
  #ifdef _REAL_
    #ifdef _QUADPRECISION_
     #define S(x) tssmq_ ## x ## _fourier_bessel_real_2d
     #define W(x) tssmq_ ## x ## _fourier_bessel_real_2d_wf128
    #else
     #define S(x) tssm_ ## x ## _fourier_bessel_real_2d
     #define W(x) tssm_ ## x ## _fourier_bessel_real_2d_wf128
    #endif 
  #else
    #ifdef _QUADPRECISION_
     #define S(x) tssmq_ ## x ## _fourier_bessel_2d
     #define W(x) tssmq_ ## x ## _fourier_bessel_2d_wf128
    #else
     #define S(x) tssm_ ## x ## _fourier_bessel_2d
     #define W(x) tssm_ ## x ## _fourier_bessel_2d_wf128
    #endif
  #endif
#endif

#ifdef _ROTSYM_
void* S(new)(int nr, __float128 r_max, int boundary_conditions);
void* W(new)(int nr, myfloat128 r_max, int boundary_conditions)
{
    return S(new)(nr,  F(r_max), boundary_conditions);
}

#else
void* S(new)(int M, int nr, int nfr, __float128 r_max, 
               int boundary_conditions, int quadrature_formula,
               char* filename, int filename_length);
void* W(new)(int M, int nr, int nfr, myfloat128 r_max, 
               int boundary_conditions, int quadrature_formula,
               char* filename, int filename_length);
{
    return S(new)(M, nr, nfr,  F(r_max),
                boundary_conditions, quadrature_formula,
                filename, filename_length);
}
#endif


void S(finalize)(void *m);
void W(finalize)(void *m)
{
    S(finalize)(m);
}

#ifndef _ROTSYM_
void S(save_coeffs)(void *m, char* filename, int filename_length);
void W(save_coeffs)(void *m, char* filename, int filename_length)
{
    S(save_coeffs)(m, filename, filename_length);
}    
#endif


void* S(new_wf)(void *psi); 
void* W(new_wf)(void *psi) 
{
    return S(new_wf)(psi);
}


void S(finalize_wf)(void *psi);
void W(finalize_wf)(void *psi)
{
    S(finalize_wf)(psi);
}


int S(is_real_space_wf)(void *psi);
int W(is_real_space_wf)(void *psi)
{
    return  S(is_real_space_wf)(psi);
}


void S(to_real_space_wf)(void *psi);
void W(to_real_space_wf)(void *psi)
{
    S(to_real_space_wf)(psi);
}


void S(to_frequency_space_wf)(void *psi);
void W(to_frequency_space_wf)(void *psi)
{
    S(to_frequency_space_wf)(psi);
}


void S(propagate_A_wf)(void *psi, __float128 dt);
void W(propagate_A_wf)(void *psi, myfloat128 dt)
{
    S(propagate_A_wf)(psi, F(dt));
}    


void S(propagate_B_wf)(void *psi, __float128 dt);
void W(propagate_B_wf)(void *psi, myfloat128 dt)
{
    S(propagate_B_wf)(psi, F(dt));
}    


void S(propagate_C_wf)(void *psi, __float128 dt);
void W(propagate_C_wf)(void *psi, myfloat128 dt)
{
    S(propagate_C_wf)(psi, F(dt));
}    


void S(add_apply_A_wf)(void *this, void *other,
                  _COMPLEX_OR_REAL_ coefficient);
void W(add_apply_A_wf)(void *this, void *other,
                  _WRAPPED_COMPLEX_OR_REAL_ coefficient)
{
    S(add_apply_A_wf)(this, other, F(coefficient));
}    


void S(save_wf)(void *psi, char* filename, int filename_length);
void W(save_wf)(void *psi, char* filename, int filename_length)
{
    S(save_wf)(psi, filename, filename_length);
}    


void S(load_wf)(void *psi, char* filename, int filename_length);
void W(load_wf)(void *psi, char* filename, int filename_length)
{
    S(load_wf)(psi, filename, filename_length);
} 


void S(copy_wf)(void *psi, void *source);
void W(copy_wf)(void *psi, void *source)
{
    S(copy_wf)(psi, source);
}    


__float128 S(norm_wf)(void *psi);
myfloat128 W(norm_wf)(void *psi)
{
    myfloat128 res;
    F(res) = S(norm_wf)(psi);
    return res;
}    


__float128 S(norm_in_frequency_space_wf)(void *psi);
myfloat128 W(norm_in_frequency_space_wf)(void *psi)
{
    myfloat128 res;
    F(res) = S(norm_in_frequency_space_wf)(psi);
    return res;
}


__float128 S(normalize_wf)(void *psi);
myfloat128 W(normalize_wf)(void *psi)
{
    myfloat128 res;
    F(res) = S(normalize_wf)(psi);
    return res;
}  


__float128 S(distance_wf)(void* psi1, void* psi2);
myfloat128 W(distance_wf)(void* psi1, void* psi2)
{
    myfloat128 res;
    F(res) = S(distance_wf)(psi1, psi2);
    return res;
}    

#ifndef _ROTSYM_
__float128 S(inner_product_wf)(void *psi, void *other);
myfloat128 W(inner_product_wf)(void *psi, void *other)
{
    myfloat128 res;
    F(res) = S(inner_product_wf)(psi, other);
    return res;
}  
#endif

void S(scale_wf)(void *psi, _COMPLEX_OR_REAL_ factor);
void W(scale_wf)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ factor)
{
     S(scale_wf)(psi, F(factor));
}     


void S(axpy_wf)(void *this, void *other, _COMPLEX_OR_REAL_ factor);
void W(axpy_wf)(void *this, void *other, _WRAPPED_COMPLEX_OR_REAL_ factor)
{
    S(axpy_wf)(this, other, F(factor));
}


#ifndef _ROTSYM
_COMPLEX_OR_REAL_ S(evaluate_wf)(void *psi, __float128 x, __float128 y);
_WRAPPED_COMPLEX_OR_REAL_ W(evaluate_wf)(void *psi, myfloat128 x, myfloat128 y)
{
    _WRAPPED_COMPLEX_OR_REAL_ res;
    F(res) = S(evaluate_wf)(psi, F(x), F(y));
    return res;
}    
#endif


void* S(get_data_wf)(void *psi, int dim[_DIM_]);
void* W(get_data_wf)(void *psi, int dim[_DIM_])
{
    return S(get_data_wf)(psi, dim);
}


/* different to the fourier version */
void* S(get_eigenvalues)(void *m, int dim[_DIM_]);
void* W(get_eigenvalues)(void *m, int dim[_DIM_])
{
    return S(get_eigenvalues)(m, dim);
}


void* S(get_nodes)(void* m, int dim[1], int which);
void* W(get_nodes)(void* m, int dim[1], int which)
{
    return S(get_nodes)(m, dim, which);
}    

/* not present in the fourier version */
void* S(get_weights)(void* m, int dim[1]);
void* W(get_weights)(void* m, int dim[1])
{
    return S(get_weights)(m, dim);
}    

void* S(get_L)(void* m, int dim[_DIM_+1]);
void* W(get_L)(void* m, int dim[_DIM_+1])
{
    return S(get_L)(m, dim);
}

int S(get_nr)(void *m); 
int W(get_nr)(void *m) 
{
    return S(get_nr)(m);
}    

__float128 S(get_rmax)(void *m);
myfloat128 W(get_rmax)(void *m)
{
    myfloat128 res;
    F(res) = S(get_rmax)(m);
    return res;
}

#ifndef _ROTSYM_
int S(get_ntheta)(void *m); 
int W(get_ntheta)(void *m) 
{
    return S(get_ntheta)(m);
}   

int S(get_nfr)(void *m); 
int W(get_nfr)(void *m) 
{
    return S(get_nfr)(m);
}  
#endif

#ifdef _ROTSYM_

#ifndef _REAL_ 
void S(rset_wf)(void *psi, __float128 (*f)(__float128));
void W(rset_wf)(void *psi, myfloat128 (*f)(myfloat128))
{
    __float128 ff(__float128 x)
    {
        myfloat128 res, xx;
        F(xx) = x;
        res = f(xx);
        return F(res);
    }
    S(rset_wf)(psi, ff);
}    

void S(rset_t_wf)(void *psi, __float128 (*f)(__float128, __float128));
void W(rset_t_wf)(void *psi, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 t)
    {
        myfloat128 res, xx, tt;
        F(xx) = x;
        F(tt) = t;
        res = f(xx, tt);
        return F(res);
    }
    S(rset_t_wf)(psi, ff);
}
#endif   

void  S(set_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128));
void  W(set_wf)(void *psi,
    _WRAPPED_COMPLEX_OR_REAL_ (*f)(myfloat128))
{
    _COMPLEX_OR_REAL_ ff(__float128 x)
    {
        _WRAPPED_COMPLEX_OR_REAL_ res;
        myfloat128 xx;
        F(xx) = x;
        res = f(xx);
        return F(res);
    }
    S(set_wf)(psi, ff);
}    

void  S(set_t_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128));
void  W(set_t_wf)(void *psi,
    _WRAPPED_COMPLEX_OR_REAL_ (*f)(myfloat128, myfloat128))
{
    _COMPLEX_OR_REAL_ ff(__float128 x, __float128 t)
    {
        _WRAPPED_COMPLEX_OR_REAL_ res;
        myfloat128 xx, tt;
        F(xx) = x;
        F(tt) = t;
        res = f(xx, tt);
        return F(res);
    }
    S(set_t_wf)(psi, ff);
}

#else

#ifndef _REAL_ 
void S(rset_wf)(void *psi, __float128 (*f)(__float128, __float128));
void W(rset_wf)(void *psi, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 y)
    {
        myfloat128 res, xx, yy;
        F(xx) = x;
        F(yy) = y;
        res = f(xx, yy);
        return F(res);
    }
    S(rset_wf)(psi, ff);
}    

void S(rset_t_wf)(void *psi, __float128 (*f)(__float128, 
                                  __float128, __float128));
void W(rset_t_wf)(void *psi, myfloat128 (*f)(myfloat128, 
                                  myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 y, __float128 t)
    {
        myfloat128 res, xx, yy, tt;
        F(xx) = x;
        F(yy) = y;
        F(tt) = t;
        res = f(xx, yy, tt);
        return F(res);
    }
    S(rset_t_wf)(psi, ff);
}
#endif   

void  S(set_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128));
void  W(set_wf)(void *psi,
    _WRAPPED_COMPLEX_OR_REAL_ (*f)(myfloat128, myfloat128))
{
    _COMPLEX_OR_REAL_ ff(__float128 x, __float128 y)
    {
        _WRAPPED_COMPLEX_OR_REAL_ res;
        myfloat128 xx, yy;
        F(xx) = x;
        F(yy) = y;
        res = f(xx, yy);
        return F(res);
    }
    S(set_wf)(psi, ff);
}    

void  S(set_t_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128));
void  W(set_t_wf)(void *psi,
    _WRAPPED_COMPLEX_OR_REAL_ (*f)(myfloat128, myfloat128, myfloat128))
{
    _COMPLEX_OR_REAL_ ff(__float128 x, __float128 y, __float128 t)
    {
        _WRAPPED_COMPLEX_OR_REAL_ res;
        myfloat128 xx, yy, tt;
        F(xx) = x;
        F(yy) = y;
        F(tt) = t;
        res = f(xx, yy, tt);
        return F(res);
    }
    S(set_t_wf)(psi, ff);
}

#endif



#endif /* _QUADPRECISION_ */
