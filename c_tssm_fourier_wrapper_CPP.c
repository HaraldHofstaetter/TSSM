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
#define S0(x,y)  x ## _real_ ## y ## D 
#define W0(x,y)  x ## _real_ ## y ## D_w 
#else
#define S0(x,y)  x ## _ ## y ## D 
#define W0(x,y)  x ## _ ## y ## D_w
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

void* W(c_new_fourier)(int nx, __float128 xmin, __float128 xmax,
#if(_DIM_>=2)
                        int ny, __float128 ymin, __float128 ymax, 
#endif
#if(_DIM_>=3)
                        int nz, __float128 zmin, __float128 zmax, 
#endif
                        int boundary_conditions);

void* S(c_new_fourier)(int nx, myfloat128 xmin, myfloat128 xmax, 
#if(_DIM_>=2)
                        int ny, myfloat128 ymin, myfloat128 ymax, 
#endif
#if(_DIM_>=3)
                        int nz, myfloat128 zmin, myfloat128 zmax, 
#endif
                        int boundary_conditions)
{
    return   W(c_new_fourier)(nx, F(xmin), F(xmax), 
#if(_DIM_>=2)
                              ny, F(ymin), F(ymax), 
#endif
#if(_DIM_>=3)
                              nz, F(zmin), F(zmax), 
#endif
                              boundary_conditions);
}


void W(c_finalize_fourier)(void *m);
void S(c_finalize_fourier)(void *m)
{
    W(c_finalize_fourier)(m);
}


void* W(c_new_wf_fourier)(void *psi); 
void* S(c_new_wf_fourier)(void *psi) 
{
    return W(c_new_wf_fourier)(psi);
}


void W(c_finalize_wf_fourier)(void *psi);
void S(c_finalize_wf_fourier)(void *psi)
{
    W(c_finalize_wf_fourier)(psi);
}


int W(c_is_real_space_wf_fourier)(void *psi);
int S(c_is_real_space_wf_fourier)(void *psi)
{
    return  W(c_is_real_space_wf_fourier)(psi);
}


void W(c_to_real_space_wf_fourier)(void *psi);
void S(c_to_real_space_wf_fourier)(void *psi)
{
    W(c_to_real_space_wf_fourier)(psi);
}


void W(c_to_frequency_space_wf_fourier)(void *psi);
void S(c_to_frequency_space_wf_fourier)(void *psi)
{
    W(c_to_frequency_space_wf_fourier)(psi);
}


void W(c_propagate_A_wf_fourier)(void *psi, __float128 dt);
void S(c_propagate_A_wf_fourier)(void *psi, myfloat128 dt)
{
    W(c_propagate_A_wf_fourier)(psi, F(dt));
}    


void W(c_propagate_B_wf_fourier)(void *psi, __float128 dt);
void S(c_propagate_B_wf_fourier)(void *psi, myfloat128 dt)
{
    W(c_propagate_B_wf_fourier)(psi, F(dt));
}    


void W(c_propagate_C_wf_fourier)(void *psi, __float128 dt);
void S(c_propagate_C_wf_fourier)(void *psi, myfloat128 dt)
{
    W(c_propagate_C_wf_fourier)(psi, F(dt));
}    


void W(c_add_apply_A_wf_fourier)(void *this, void *other,
                  _COMPLEX_OR_REAL_ coefficient);
void S(c_add_apply_A_wf_fourier)(void *this, void *other,
                  _WRAPPED_COMPLEX_OR_REAL_ coefficient)
{
    W(c_add_apply_A_wf_fourier)(this, other, F(coefficient));
}    


void W(c_save_wf_fourier)(void *psi, char* filename, int filename_length);
void S(c_save_wf_fourier)(void *psi, char* filename, int filename_length)
{
    W(c_save_wf_fourier)(psi, filename, filename_length);
}    


void W(c_load_wf_fourier)(void *psi, char* filename, int filename_length);
void S(c_load_wf_fourier)(void *psi, char* filename, int filename_length)
{
    W(c_load_wf_fourier)(psi, filename, filename_length);
} 


void W(c_copy_wf_fourier)(void *psi, void *source);
void S(c_copy_wf_fourier)(void *psi, void *source)
{
    S(c_copy_wf_fourier)(psi, source);
}    


__float128 W(c_norm2_wf_fourier)(void *psi);
myfloat128 S(c_norm2_wf_fourier)(void *psi)
{
    myfloat128 res;
    F(res) = W(c_norm2_wf_fourier)(psi);
    return res;
}    


__float128 W(c_norm2_in_frequency_space_wf_fourier)(void *psi);
myfloat128 S(c_norm2_in_frequency_space_wf_fourier)(void *psi)
{
    myfloat128 res;
    F(res) = W(c_norm2_in_frequency_space_wf_fourier)(psi);
    return res;
}


__float128 W(c_normalize_wf_fourier)(void *psi);
myfloat128 S(c_normalize_wf_fourier)(void *psi)
{
    myfloat128 res;
    F(res) = W(c_normalize_wf_fourier)(psi);
    return res;
}  


__float128 W(c_distance_wf_fourier)(void* psi1, void* psi2);
myfloat128 S(c_distance_wf_fourier)(void* psi1, void* psi2)
{
    myfloat128 res;
    F(res) = W(c_distance_wf_fourier)(psi1, psi2);
    return res;
}    


void W(c_scale_wf_fourier)(void *psi, _COMPLEX_OR_REAL_ factor);
void S(c_scale_wf_fourier)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ factor)
{
     W(c_scale_wf_fourier)(psi, F(factor));
}     


void W(c_axpy_wf_fourier)(void *this, void *other, _COMPLEX_OR_REAL_ factor);
void S(c_axpy_wf_fourier)(void *this, void *other, _WRAPPED_COMPLEX_OR_REAL_ factor)
{
    W(c_axpy_wf_fourier)(this, other, F(factor));
}


void* W(c_get_data_wf_fourier)(void *psi, int dim[_DIM_]);
void* S(c_get_data_wf_fourier)(void *psi, int dim[_DIM_])
{
    return W(c_get_data_wf_fourier)(psi, dim);
}


void* W(c_get_eigenvalues_fourier)(void *m, int dim[1], int which);
void* S(c_get_eigenvalues_fourier)(void *m, int dim[1], int which)
{
    return W(c_get_eigenvalues_fourier)(m, dim, which);
}


void* W(c_get_nodes_fourier)(void* m, int dim[1], int which);
void* S(c_get_nodes_fourier)(void* m, int dim[1], int which)
{
    return W(c_get_nodes_fourier)(m, dim, which);
}    


int W(c_get_nx_fourier)(void *m); 
int S(c_get_nx_fourier)(void *m) 
{
    return W(c_get_nx_fourier)(m);
}    


__float128 W(c_get_xmin_fourier)(void *m);
myfloat128 S(c_get_xmin_fourier)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_xmin_fourier)(m);
    return res;
}


__float128 W(c_get_xmax_fourier)(void *m);
myfloat128 S(c_get_xmax_fourier)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_xmax_fourier)(m);
    return res;
}


#if(_DIM_>=2)
int W(c_get_ny_fourier)(void *m);
int S(c_get_ny_fourier)(void *m) 
{
    return W(c_get_ny_fourier)(m);
}    


__float128 W(c_get_ymin_fourier)(void *m);
myfloat128 S(c_get_ymin_fourier)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_ymin_fourier)(m);
    return res;
}


__float128 W(c_get_ymax_fourier)(void *m);
myfloat128 S(c_get_ymax_fourier)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_ymax_fourier)(m);
    return res;
}
#endif


#if(_DIM_>=3)
int W(c_get_nz_fourier)(void *m);
int S(c_get_nz_fourier)(void *m) 
{
    return W(c_get_nz_fourier)(m);
}    


__float128 W(c_get_zmin_fourier)(void *m);
myfloat128 S(c_get_zmin_fourier)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_zmin_fourier)(m);
    return res;
}


__float128 W(c_get_zmax_fourier)(void *m);
myfloat128 S(c_get_zmax_fourier)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_zmax_fourier)(m);
    return res;
}
#endif


#ifndef _REAL_ 
#if(_DIM_==1)
void W(c_rset_wf_fourier)(void *psi, __float128 (*f)(__float128));
void S(c_rset_wf_fourier)(void *psi, myfloat128 (*f)(myfloat128))
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
    W(c_rset_wf_fourier)(psi, ff);
}    

void W(c_rset_t_wf_fourier)(void *psi, __float128 (*f)(__float128, __float128));
void S(c_rset_t_wf_fourier)(void *psi, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 t)
    {
        myfloat128 res, xx, tt;
        F(xx) = x;
        F(tt) = t;
        res = f(xx, tt);
        return F(res);
    }
    W(c_rset_t_wf_fourier)(psi, ff);
}

#elif(_DIM_==2)
void W(c_rset_wf_fourier)(void *psi, __float128 (*f)(__float128, __float128));
void S(c_rset_wf_fourier)(void *psi, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 y)
    {
        myfloat128 res, xx, yy;
        F(xx) = x;
        F(yy) = y;
        res = f(xx, yy);
        return F(res);
    }
    W(c_rset_wf_fourier)(psi, ff);
}    

void W(c_rset_t_wf_fourier)(void *psi, __float128 (*f)(__float128, 
                                  __float128, __float128));
void S(c_rset_t_wf_fourier)(void *psi, myfloat128 (*f)(myfloat128, 
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
    W(c_rset_t_wf_fourier)(psi, ff);
}

#elif(_DIM_==3)
void W(c_rset_wf_fourier)(void *psi, __float128 (*f)(__float128, 
                                __float128, __float128));
void S(c_rset_wf_fourier)(void *psi, myfloat128 (*f)(myfloat128, 
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
    W(c_rset_wf_fourier)(psi, ff);
}    

void W(c_rset_t_wf_fourier)(void *psi, __float128 (*f)(__float128, 
                                  __float128, __float128, __float128));
void S(c_rset_t_wf_fourier)(void *psi, myfloat128 (*f)(myfloat128, 
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
    W(c_rset_t_wf_fourier)(psi, ff);
}

#endif   
#endif   


#if(_DIM_==1)
void  W(c_set_wf_fourier)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128));
void  S(c_set_wf_fourier)(void *psi,
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
    W(c_set_wf_fourier)(psi, ff);
}    

void  W(c_set_t_wf_fourier)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128));
void  S(c_set_t_wf_fourier)(void *psi,
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
    W(c_set_t_wf_fourier)(psi, ff);
}

#elif(_DIM_==2)
void  W(c_set_wf_fourier)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128));
void  S(c_set_wf_fourier)(void *psi,
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
    W(c_set_wf_fourier)(psi, ff);
}    

void  W(c_set_t_wf_fourier)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128));
void  S(c_set_t_wf_fourier)(void *psi,
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
    W(c_set_t_wf_fourier)(psi, ff);
}

#elif(_DIM_==3)
void  W(c_set_wf_fourier)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128));
void  S(c_set_wf_fourier)(void *psi,
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
    W(c_set_wf_fourier)(psi, ff);
}    

void  W(c_set_t_wf_fourier)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128, __float128));
void  S(c_set_t_wf_fourier)(void *psi,
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
    W(c_set_t_wf_fourier)(psi, ff);
}

#endif   




#endif