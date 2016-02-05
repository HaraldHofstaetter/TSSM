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
#define S0(x,y)  x ## _fourier_real_ ## y ## d 
#define W0(x,y)  x ## _fourier_real_ ## y ## d_wf128 
#else
#define S0(x,y)  x ## _fourier_ ## y ## d 
#define W0(x,y)  x ## _fourier_ ## y ## d_wf128
#endif
#define S1(x,y) S0(x,y)
#define S(x) S1(x,_DIM_)
#define W1(x,y) W0(x,y)
#define W(x) W1(x,_DIM_)
#ifdef _REAL_
 #define _WAVE_FUNCTION_ real_wave_function
 #define _COMPLEX_OR_REAL_ __float128
 #define _WRAPPED_COMPLEX_OR_REAL_ myfloat128
#else
 #define _WAVE_FUNCTION_ wave_function
 #define _COMPLEX_OR_REAL_ __complex128
 #define _WRAPPED_COMPLEX_OR_REAL_ mycomplex128
#endif

void* S(c_new)(int nx, __float128 xmin, __float128 xmax,
#if(_DIM_>=2)
                        int ny, __float128 ymin, __float128 ymax, 
#endif
#if(_DIM_>=3)
                        int nz, __float128 zmin, __float128 zmax, 
#endif
                        int boundary_conditions);
__attribute__((visibility("default")))
void* W(c_new)(int nx, myfloat128 xmin, myfloat128 xmax, 
#if(_DIM_>=2)
                        int ny, myfloat128 ymin, myfloat128 ymax, 
#endif
#if(_DIM_>=3)
                        int nz, myfloat128 zmin, myfloat128 zmax, 
#endif
                        int boundary_conditions)
{
    return   S(c_new)(nx, F(xmin), F(xmax), 
#if(_DIM_>=2)
                              ny, F(ymin), F(ymax), 
#endif
#if(_DIM_>=3)
                              nz, F(zmin), F(zmax), 
#endif
                              boundary_conditions);
}


void S(c_finalize)(void *m);
void W(c_finalize)(void *m)
{
    S(c_finalize)(m);
}


void* S(c_new_wf)(void *psi); 
void* W(c_new_wf)(void *psi) 
{
    return S(c_new_wf)(psi);
}


void S(c_finalize_wf)(void *psi);
void W(c_finalize_wf)(void *psi)
{
    S(c_finalize_wf)(psi);
}


int S(c_is_real_space_wf)(void *psi);
int W(c_is_real_space_wf)(void *psi)
{
    return  S(c_is_real_space_wf)(psi);
}


void S(c_to_real_space_wf)(void *psi);
void W(c_to_real_space_wf)(void *psi)
{
    S(c_to_real_space_wf)(psi);
}


void S(c_to_frequency_space_wf)(void *psi);
void W(c_to_frequency_space_wf)(void *psi)
{
    S(c_to_frequency_space_wf)(psi);
}


void S(c_propagate_A_wf)(void *psi, __float128 dt);
void W(c_propagate_A_wf)(void *psi, myfloat128 dt)
{
    S(c_propagate_A_wf)(psi, F(dt));
}    


void S(c_propagate_B_wf)(void *psi, __float128 dt);
void W(c_propagate_B_wf)(void *psi, myfloat128 dt)
{
    S(c_propagate_B_wf)(psi, F(dt));
}    


void S(c_propagate_C_wf)(void *psi, __float128 dt);
void W(c_propagate_C_wf)(void *psi, myfloat128 dt)
{
    S(c_propagate_C_wf)(psi, F(dt));
}    


void S(c_add_apply_A_wf)(void *this, void *other,
                  _COMPLEX_OR_REAL_ coefficient);
void W(c_add_apply_A_wf)(void *this, void *other,
                  _WRAPPED_COMPLEX_OR_REAL_ coefficient)
{
    S(c_add_apply_A_wf)(this, other, F(coefficient));
}    


void S(c_save_wf)(void *psi, char* filename, int filename_length);
void W(c_save_wf)(void *psi, char* filename, int filename_length)
{
    S(c_save_wf)(psi, filename, filename_length);
}    


void S(c_load_wf)(void *psi, char* filename, int filename_length);
void W(c_load_wf)(void *psi, char* filename, int filename_length)
{
    S(c_load_wf)(psi, filename, filename_length);
} 


void S(c_copy_wf)(void *psi, void *source);
void W(c_copy_wf)(void *psi, void *source)
{
    W(c_copy_wf)(psi, source);
}    


__float128 S(c_norm2_wf)(void *psi);
myfloat128 W(c_norm2_wf)(void *psi)
{
    myfloat128 res;
    F(res) = S(c_norm2_wf)(psi);
    return res;
}    


__float128 S(c_norm2_in_frequency_space_wf)(void *psi);
myfloat128 W(c_norm2_in_frequency_space_wf)(void *psi)
{
    myfloat128 res;
    F(res) = S(c_norm2_in_frequency_space_wf)(psi);
    return res;
}


__float128 S(c_normalize_wf)(void *psi);
myfloat128 W(c_normalize_wf)(void *psi)
{
    myfloat128 res;
    F(res) = S(c_normalize_wf)(psi);
    return res;
}  


__float128 S(c_distance_wf)(void* psi1, void* psi2);
myfloat128 W(c_distance_wf)(void* psi1, void* psi2)
{
    myfloat128 res;
    F(res) = S(c_distance_wf)(psi1, psi2);
    return res;
}    


void S(c_scale_wf)(void *psi, _COMPLEX_OR_REAL_ factor);
void W(c_scale_wf)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ factor)
{
     S(c_scale_wf)(psi, F(factor));
}     


void S(c_axpy_wf)(void *this, void *other, _COMPLEX_OR_REAL_ factor);
void W(c_axpy_wf)(void *this, void *other, _WRAPPED_COMPLEX_OR_REAL_ factor)
{
    S(c_axpy_wf)(this, other, F(factor));
}


void* S(c_get_data_wf)(void *psi, int dim[_DIM_]);
void* W(c_get_data_wf)(void *psi, int dim[_DIM_])
{
    return S(c_get_data_wf)(psi, dim);
}


void* S(c_get_eigenvalues)(void *m, int dim[1], int which);
void* W(c_get_eigenvalues)(void *m, int dim[1], int which)
{
    return S(c_get_eigenvalues)(m, dim, which);
}


void* S(c_get_nodes)(void* m, int dim[1], int which);
void* W(c_get_nodes)(void* m, int dim[1], int which)
{
    return S(c_get_nodes)(m, dim, which);
}    


int S(c_get_nx)(void *m); 
int W(c_get_nx)(void *m) 
{
    return S(c_get_nx)(m);
}    


__float128 S(c_get_xmin)(void *m);
myfloat128 W(c_get_xmin)(void *m)
{
    myfloat128 res;
    F(res) = S(c_get_xmin)(m);
    return res;
}


__float128 S(c_get_xmax)(void *m);
myfloat128 W(c_get_xmax)(void *m)
{
    myfloat128 res;
    F(res) = S(c_get_xmax)(m);
    return res;
}


#if(_DIM_>=2)
int S(c_get_ny)(void *m);
int W(c_get_ny)(void *m) 
{
    return S(c_get_ny)(m);
}    


__float128 S(c_get_ymin)(void *m);
myfloat128 W(c_get_ymin)(void *m)
{
    myfloat128 res;
    F(res) = S(c_get_ymin)(m);
    return res;
}


__float128 S(c_get_ymax)(void *m);
myfloat128 W(c_get_ymax)(void *m)
{
    myfloat128 res;
    F(res) = S(c_get_ymax)(m);
    return res;
}
#endif


#if(_DIM_>=3)
int S(c_get_nz)(void *m);
int W(c_get_nz)(void *m) 
{
    return S(c_get_nz)(m);
}    


__float128 S(c_get_zmin)(void *m);
myfloat128 W(c_get_zmin)(void *m)
{
    myfloat128 res;
    F(res) = S(c_get_zmin)(m);
    return res;
}


__float128 S(c_get_zmax)(void *m);
myfloat128 W(c_get_zmax)(void *m)
{
    myfloat128 res;
    F(res) = S(c_get_zmax)(m);
    return res;
}
#endif


#ifndef _REAL_ 
#if(_DIM_==1)
void S(c_rset_wf)(void *psi, __float128 (*f)(__float128));
void W(c_rset_wf)(void *psi, myfloat128 (*f)(myfloat128))
{
    /* Note: Here and below we use nested functions, these are supported as
       extension in GNU C, but not in ISO standard C, see 
       https://gcc.gnu.org/onlinedocs/gcc/Nested-Functions.html
    */   
    __float128 ff(__float128 x)
    {
        myfloat128 res, xx;
        F(xx) = x;
        res = f(xx);
        return F(res);
    }
    S(c_rset_wf)(psi, ff);
}    

void S(c_rset_t_wf)(void *psi, __float128 (*f)(__float128, __float128));
void W(c_rset_t_wf)(void *psi, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 t)
    {
        myfloat128 res, xx, tt;
        F(xx) = x;
        F(tt) = t;
        res = f(xx, tt);
        return F(res);
    }
    S(c_rset_t_wf)(psi, ff);
}

#elif(_DIM_==2)
void S(c_rset_wf)(void *psi, __float128 (*f)(__float128, __float128));
void W(c_rset_wf)(void *psi, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 y)
    {
        myfloat128 res, xx, yy;
        F(xx) = x;
        F(yy) = y;
        res = f(xx, yy);
        return F(res);
    }
    S(c_rset_wf)(psi, ff);
}    

void S(c_rset_t_wf)(void *psi, __float128 (*f)(__float128, 
                                  __float128, __float128));
void W(c_rset_t_wf)(void *psi, myfloat128 (*f)(myfloat128, 
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
    S(c_rset_t_wf)(psi, ff);
}

#elif(_DIM_==3)
void S(c_rset_wf)(void *psi, __float128 (*f)(__float128, 
                                __float128, __float128));
void W(c_rset_wf)(void *psi, myfloat128 (*f)(myfloat128, 
                                myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 y, __float128 z)
    {
        myfloat128 res, xx, yy, zz;
        F(xx) = x;
        F(yy) = y;
        F(zz) = z;
        res = f(xx, yy, zz);
        return F(res);
    }
    S(c_rset_wf)(psi, ff);
}    

void S(c_rset_t_wf)(void *psi, __float128 (*f)(__float128, 
                                  __float128, __float128, __float128));
void W(c_rset_t_wf)(void *psi, myfloat128 (*f)(myfloat128, 
                                  myfloat128, myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 y, __float128 z, __float128 t)
    {
        myfloat128 res, xx, yy, zz, tt;
        F(xx) = x;
        F(yy) = y;
        F(zz) = z;
        F(tt) = t;
        res = f(xx, yy, zz, tt);
        return F(res);
    }
    S(c_rset_t_wf)(psi, ff);
}

#endif   
#endif   


#if(_DIM_==1)
void  S(c_set_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128));
void  W(c_set_wf)(void *psi,
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
    S(c_set_wf)(psi, ff);
}    

void  S(c_set_t_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128));
void  W(c_set_t_wf)(void *psi,
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
    S(c_set_t_wf)(psi, ff);
}

#elif(_DIM_==2)
void  S(c_set_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128));
void  W(c_set_wf)(void *psi,
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
    S(c_set_wf)(psi, ff);
}    

void  S(c_set_t_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128));
void  W(c_set_t_wf)(void *psi,
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
    S(c_set_t_wf)(psi, ff);
}

#elif(_DIM_==3)
void  S(c_set_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128));
void  W(c_set_wf)(void *psi,
    _WRAPPED_COMPLEX_OR_REAL_ (*f)(myfloat128, myfloat128, myfloat128))
{
    _COMPLEX_OR_REAL_ ff(__float128 x, __float128 y, __float128 z)
    {
        _WRAPPED_COMPLEX_OR_REAL_ res;
        myfloat128 xx, yy, zz;
        F(xx) = x;
        F(yy) = y;
        F(zz) = z;
        res = f(xx, yy, zz);
        return F(res);
    }
    S(c_set_wf)(psi, ff);
}    

void  S(c_set_t_wf)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128, __float128));
void  W(c_set_t_wf)(void *psi,
    _WRAPPED_COMPLEX_OR_REAL_ (*f)(myfloat128, myfloat128, myfloat128, myfloat128))
{
    _COMPLEX_OR_REAL_ ff(__float128 x, __float128 y, __float128 z, __float128 t)
    {
        _WRAPPED_COMPLEX_OR_REAL_ res;
        myfloat128 xx, yy, zz, tt;
        F(xx) = x;
        F(yy) = y;
        F(zz) = z;
        F(tt) = t;
        res = f(xx, yy, zz, tt);
        return F(res);
    }
    S(c_set_t_wf)(psi, ff);
}

#endif   




#endif
