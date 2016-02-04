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

#ifdef _HERMITE_
#ifdef _REAL_
#define S0(x,y)  x ## _hermite_real_ ## y ## d 
#define W0(x,y)  x ## _hermite_real_ ## y ## d_w 
#else
#define S0(x,y)  x ## _hermite_ ## y ## d 
#define W0(x,y)  x ## _hermite_ ## y ## d_w 
#endif
#else
#ifdef _REAL_
#define S0(x,y)  x ## _real_ ## y ## d 
#define W0(x,y)  x ## _real_ ## y ## d_w 
#else
#define S0(x,y)  x ## _ ## y ## d 
#define W0(x,y)  x ## _ ## y ## d_w
#endif
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


#ifdef _HERMITE_

#else
#if(_DIM==1)                        
void* W(c_new_schroedinger)(
#ifdef _HERMITE_
                        int nx, __float128 omega_x,
#else
                        int nx, __float128 xmin, __float128 xmax,
#endif                        
                        __float128 hbar, __float128 mass, 
                        __float128 (*potential)(__float128),
                        int with_potential, __float128 cubic_coupling,
                        int boundary_conditions);

void* S(c_new_schroedinger)(
#ifdef _HERMITE_
                        int nx, myfloat128 _omega_x, 
#else
                        int nx, myfloat128 xmin, myfloat128 xmax, 
#endif                        
                        myfloat hbar, myfloat mass, 
                        myfloat (*potential)(myfloat),
                        int with_potential, myfloat cubic_coupling,
                        int boundary_conditions)
{
    __float128 pp(__float128 x)
    {
        myfloat128 res, xx;
        F(xx) = x;
        res = potential(xx);
        return F(res);
    }        
    return   W(c_new_schroedinger)(
#ifdef _HERMITE_
                              nx, F(omega_x),
#else
                              nx, F(xmin), F(xmax), 
#endif                        
                              F(hbar), F(mass), pp,
                              with_potential, F(cubic_coupling)
                              boundary_conditions);
}
#elif(_DIM==2)                        
void* W(c_new_schroedinger)(
#ifdef _HERMITE_
                        int nx, __float128 omega_x,
                        int ny, __float128 omega_y,
#else
                        int nx, __float128 xmin, __float128 xmax,
                        int ny, __float128 ymin, __float128 ymax, 
#endif             
                        __float128 hbar, __float128 mass, 
                        __float128 (*potential)(__float128,__float128),
                        int with_potential, __float128 cubic_coupling,
                        int boundary_conditions);

void* S(c_new_schroedinger)(
#ifdef _HERMITE_
                        int nx, myfloat128 _omega_x, 
                        int ny, myfloat128 _omega_y, 
#else
                        int nx, myfloat128 xmin, myfloat128 xmax, 
                        int ny, myfloat128 ymin, myfloat128 ymax, 
#endif                              
                        myfloat hbar, myfloat mass, 
                        myfloat (*potential)(myfloat,myfloat),
                        int with_potential, myfloat cubic_coupling,
                        int boundary_conditions)
{
    __float128 pp(__float128 x, __float128 y)
    {
        myfloat128 res, xx, yy;
        F(xx) = x;
        F(yy) = y;
        res = potential(xx, yy);
        return F(res);
    }        
    return   W(c_new_schroedinger)(
#ifdef _HERMITE_
                              nx, F(omega_x),
                              ny, F(omega_y),
#else
                              nx, F(xmin), F(xmax), 
                              ny, F(ymin), F(ymax), 
#endif                        
                              F(hbar), F(mass), pp,
                              with_potential, F(cubic_coupling)
                              boundary_conditions);
}
#elif(_DIM==3)                        
void* W(c_new_schroedinger)(
#ifdef _HERMITE_
                        int nx, __float128 omega_x,
                        int ny, __float128 omega_y,
                        int nz, __float128 omega_z,
#else
                        int nx, __float128 xmin, __float128 xmax,
                        int ny, __float128 ymin, __float128 ymax, 
                        int nz, __float128 zmin, __float128 zmax, 
#endif             
                        __float128 hbar, __float128 mass, 
                        __float128 (*potential)(__float128,__float128,__float128),
                        int with_potential, __float128 cubic_coupling,
                        int boundary_conditions);

void* S(c_new_schroedinger)(
#ifdef _HERMITE_
                        int nx, myfloat128 _omega_x, 
                        int ny, myfloat128 _omega_y, 
                        int nz, myfloat128 _omega_z, 
#else
                        int nx, myfloat128 xmin, myfloat128 xmax, 
                        int ny, myfloat128 ymin, myfloat128 ymax, 
                        int nz, myfloat128 zmin, myfloat128 zmax, 
#endif                              
                        myfloat hbar, myfloat mass, 
                        myfloat (*potential)(myfloat,myfloat,myfloat),
                        int with_potential, myfloat cubic_coupling,
                        int boundary_conditions)
{
    __float128 pp(__float128 x, __float128 y, __float128 z)
    {
        myfloat128 res, xx, yy, zz;
        F(xx) = x;
        F(yy) = y;
        F(zz) = z;
        res = potential(xx, yy, zz);
        return F(res);
    }        
    return   W(c_new_schroedinger)(
#ifdef _HERMITE_
                              nx, F(omega_x),
                              ny, F(omega_y),
                              nz, F(omega_z),
#else
                              nx, F(xmin), F(xmax), 
                              ny, F(ymin), F(ymax), 
                              nz, F(zmin), F(zmax), 
#endif                        
                              F(hbar), F(mass), pp,
                              with_potential, F(cubic_coupling)
                              boundary_conditions);
}

#endif
#endif


void W(c_finalize_schroedinger)(void *m);
void S(c_finalize_schroedinger)(void *m)
{
    W(c_finalize_schroedinger)(m);
}


void* W(c_new_wf_schroedinger)(void *psi); 
void* S(c_new_wf_schroedinger)(void *psi) 
{
    return W(c_new_wf_schroedinger)(psi);
}


void W(c_finalize_wf_schroedinger)(void *psi);
void S(c_finalize_wf_schroedinger)(void *psi)
{
    W(c_finalize_wf_schroedinger)(psi);
}


int W(c_is_real_space_wf_schroedinger)(void *psi);
int S(c_is_real_space_wf_schroedinger)(void *psi)
{
    return  W(c_is_real_space_wf_schroedinger)(psi);
}


void W(c_to_real_space_wf_schroedinger)(void *psi);
void S(c_to_real_space_wf_schroedinger)(void *psi)
{
    W(c_to_real_space_wf_schroedinger)(psi);
}


void W(c_to_frequency_space_wf_schroedinger)(void *psi);
void S(c_to_frequency_space_wf_schroedinger)(void *psi)
{
    W(c_to_frequency_space_wf_schroedinger)(psi);
}

void W(c_propagate_A_wf_schroedinger)(void *psi, _COMPLEX_OR_REAL_ dt);
void S(c_propagate_A_wf_schroedinger)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ dt)
{
    W(c_propagate_A_wf_schroedinger)(psi, F(dt));
}    

void W(c_propagate_B_wf_schroedinger)(void *psi, _COMPLEX_OR_REAL_ dt);
void S(c_propagate_B_wf_schroedinger)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ dt)
{
    W(c_propagate_B_wf_schroedinger)(psi, F(dt));
}    

void W(c_propagate_C_wf_schroedinger)(void *psi, _COMPLEX_OR_REAL_ dt);
void S(c_propagate_C_wf_schroedinger)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ dt)
{
    W(c_propagate_C_wf_schroedinger)(psi, F(dt));
}    


void W(c_add_apply_A_wf_schroedinger)(void *this, void *other,
                  _COMPLEX_OR_REAL_ coefficient);
void S(c_add_apply_A_wf_schroedinger)(void *this, void *other,
                  _WRAPPED_COMPLEX_OR_REAL_ coefficient)
{
    W(c_add_apply_A_wf_schroedinger)(this, other, F(coefficient));
}    


void W(c_save_wf_schroedinger)(void *psi, char* filename, int filename_length);
void S(c_save_wf_schroedinger)(void *psi, char* filename, int filename_length)
{
    W(c_save_wf_schroedinger)(psi, filename, filename_length);
}    


void W(c_load_wf_schroedinger)(void *psi, char* filename, int filename_length);
void S(c_load_wf_schroedinger)(void *psi, char* filename, int filename_length)
{
    W(c_load_wf_schroedinger)(psi, filename, filename_length);
} 


void W(c_copy_wf_schroedinger)(void *psi, void *source);
void S(c_copy_wf_schroedinger)(void *psi, void *source)
{
    S(c_copy_wf_schroedinger)(psi, source);
}    


__float128 W(c_norm2_wf_schroedinger)(void *psi);
myfloat128 S(c_norm2_wf_schroedinger)(void *psi)
{
    myfloat128 res;
    F(res) = W(c_norm2_wf_schroedinger)(psi);
    return res;
}    


__float128 W(c_norm2_in_frequency_space_wf_schroedinger)(void *psi);
myfloat128 S(c_norm2_in_frequency_space_wf_schroedinger)(void *psi)
{
    myfloat128 res;
    F(res) = W(c_norm2_in_frequency_space_wf_schroedinger)(psi);
    return res;
}


__float128 W(c_normalize_wf_schroedinger)(void *psi);
myfloat128 S(c_normalize_wf_schroedinger)(void *psi)
{
    myfloat128 res;
    F(res) = W(c_normalize_wf_schroedinger)(psi);
    return res;
}  


__float128 W(c_distance_wf_schroedinger)(void* psi1, void* psi2);
myfloat128 S(c_distance_wf_schroedinger)(void* psi1, void* psi2)
{
    myfloat128 res;
    F(res) = W(c_distance_wf_schroedinger)(psi1, psi2);
    return res;
}    


void W(c_scale_wf_schroedinger)(void *psi, _COMPLEX_OR_REAL_ factor);
void S(c_scale_wf_schroedinger)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ factor)
{
     W(c_scale_wf_schroedinger)(psi, F(factor));
}     


void W(c_axpy_wf_schroedinger)(void *this, void *other, _COMPLEX_OR_REAL_ factor);
void S(c_axpy_wf_schroedinger)(void *this, void *other, _WRAPPED_COMPLEX_OR_REAL_ factor)
{
    W(c_axpy_wf_schroedinger)(this, other, F(factor));
}


void* W(c_get_data_wf_schroedinger)(void *psi, int dim[_DIM_]);
void* S(c_get_data_wf_schroedinger)(void *psi, int dim[_DIM_])
{
    return W(c_get_data_wf_schroedinger)(psi, dim);
}


void* W(c_get_eigenvalues_schroedinger)(void *m, int dim[1], int which);
void* S(c_get_eigenvalues_schroedinger)(void *m, int dim[1], int which)
{
    return W(c_get_eigenvalues_schroedinger)(m, dim, which);
}


void* W(c_get_nodes_schroedinger)(void* m, int dim[1], int which);
void* S(c_get_nodes_schroedinger)(void* m, int dim[1], int which)
{
    return W(c_get_nodes_schroedinger)(m, dim, which);
}    

#ifdef _HERMITE_
void* W(c_get_weights_schroedinger)(void* m, int dim[1], int which);
void* S(c_get_weights_schroedinger)(void* m, int dim[1], int which)
{
    return W(c_get_weights_schroedinger)(m, dim, which);
}    

void* W(c_get_H)(void* m, int dim[2], int which);
void* S(c_get_H)(void* m, int dim[2], int which)
{
    return W(c_get_H)(m, dim, which);
}   
#endif

int W(c_get_nx_schroedinger)(void *m); 
int S(c_get_nx_schroedinger)(void *m) 
{
    return W(c_get_nx_schroedinger)(m);
}    

#ifdef _HERMITE_
__float128 W(c_get_omega_x_schroedinger)(void *m);
myfloat128 S(c_get_omega_x_schroedinger)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_omega_x_schroedinger)(m);
    return res;
}

#else
__float128 W(c_get_xmin_schroedinger)(void *m);
myfloat128 S(c_get_xmin_schroedinger)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_xmin_schroedinger)(m);
    return res;
}

__float128 W(c_get_xmax_schroedinger)(void *m);
myfloat128 S(c_get_xmax_schroedinger)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_xmax_schroedinger)(m);
    return res;
}
#endif


#if(_DIM_>=2)
int W(c_get_ny_schroedinger)(void *m);
int S(c_get_ny_schroedinger)(void *m) 
{
    return W(c_get_ny_schroedinger)(m);
}    

#ifdef _HERMITE_
__float128 W(c_get_omega_y_schroedinger)(void *m);
myfloat128 S(c_get_omega_y_schroedinger)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_omega_y_schroedinger)(m);
    return res;
}
#else
__float128 W(c_get_ymin_schroedinger)(void *m);
myfloat128 S(c_get_ymin_schroedinger)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_ymin_schroedinger)(m);
    return res;
}

__float128 W(c_get_ymax_schroedinger)(void *m);
myfloat128 S(c_get_ymax_schroedinger)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_ymax_schroedinger)(m);
    return res;
}
#endif
#endif


#if(_DIM_>=3)
int W(c_get_nz_schroedinger)(void *m);
int S(c_get_nz_schroedinger)(void *m) 
{
    return W(c_get_nz_schroedinger)(m);
}    

#ifdef _HERMITE_
__float128 W(c_get_omega_z_schroedinger)(void *m);
myfloat128 S(c_get_omega_z_schroedinger)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_omega_z_schroedinger)(m);
    return res;
}
#else
__float128 W(c_get_zmin_schroedinger)(void *m);
myfloat128 S(c_get_zmin_schroedinger)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_zmin_schroedinger)(m);
    return res;
}


__float128 W(c_get_zmax_schroedinger)(void *m);
myfloat128 S(c_get_zmax_schroedinger)(void *m)
{
    myfloat128 res;
    F(res) = W(c_get_zmax_schroedinger)(m);
    return res;
}
#endif
#endif


#ifndef _REAL_ 
#if(_DIM_==1)
void W(c_rset_wf_schroedinger)(void *psi, __float128 (*f)(__float128));
void S(c_rset_wf_schroedinger)(void *psi, myfloat128 (*f)(myfloat128))
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
    W(c_rset_wf_schroedinger)(psi, ff);
}    

void W(c_rset_t_wf_schroedinger)(void *psi, __float128 (*f)(__float128, __float128));
void S(c_rset_t_wf_schroedinger)(void *psi, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 t)
    {
        myfloat128 res, xx, tt;
        F(xx) = x;
        F(tt) = t;
        res = f(xx, tt);
        return F(res);
    }
    W(c_rset_t_wf_schroedinger)(psi, ff);
}

#elif(_DIM_==2)
void W(c_rset_wf_schroedinger)(void *psi, __float128 (*f)(__float128, __float128));
void S(c_rset_wf_schroedinger)(void *psi, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 y)
    {
        myfloat128 res, xx, yy;
        F(xx) = x;
        F(yy) = y;
        res = f(xx, yy);
        return F(res);
    }
    W(c_rset_wf_schroedinger)(psi, ff);
}    

void W(c_rset_t_wf_schroedinger)(void *psi, __float128 (*f)(__float128, 
                                  __float128, __float128));
void S(c_rset_t_wf_schroedinger)(void *psi, myfloat128 (*f)(myfloat128, 
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
    W(c_rset_t_wf_schroedinger)(psi, ff);
}

#elif(_DIM_==3)
void W(c_rset_wf_schroedinger)(void *psi, __float128 (*f)(__float128, 
                                __float128, __float128));
void S(c_rset_wf_schroedinger)(void *psi, myfloat128 (*f)(myfloat128, 
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
    W(c_rset_wf_schroedinger)(psi, ff);
}    

void W(c_rset_t_wf_schroedinger)(void *psi, __float128 (*f)(__float128, 
                                  __float128, __float128, __float128));
void S(c_rset_t_wf_schroedinger)(void *psi, myfloat128 (*f)(myfloat128, 
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
    W(c_rset_t_wf_schroedinger)(psi, ff);
}

#endif   
#endif   


#if(_DIM_==1)
void  W(c_set_wf_schroedinger)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128));
void  S(c_set_wf_schroedinger)(void *psi,
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
    W(c_set_wf_schroedinger)(psi, ff);
}    

void  W(c_set_t_wf_schroedinger)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128));
void  S(c_set_t_wf_schroedinger)(void *psi,
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
    W(c_set_t_wf_schroedinger)(psi, ff);
}

#elif(_DIM_==2)
void  W(c_set_wf_schroedinger)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128));
void  S(c_set_wf_schroedinger)(void *psi,
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
    W(c_set_wf_schroedinger)(psi, ff);
}    

void  W(c_set_t_wf_schroedinger)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128));
void  S(c_set_t_wf_schroedinger)(void *psi,
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
    W(c_set_t_wf_schroedinger)(psi, ff);
}

#elif(_DIM_==3)
void  W(c_set_wf_schroedinger)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128));
void  S(c_set_wf_schroedinger)(void *psi,
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
    W(c_set_wf_schroedinger)(psi, ff);
}    

void  W(c_set_t_wf_schroedinger)(void *psi,
    _COMPLEX_OR_REAL_ (*f)(__float128, __float128, __float128, __float128));
void  S(c_set_t_wf_schroedinger)(void *psi,
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
    W(c_set_t_wf_schroedinger)(psi, ff);
}

#endif   



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! schroedinger specific methods !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#if(_DIM_==1)
void W(c_set_potential_schroedinger)(void *m, __float128 (*f)(__float128));
void S(c_set_potential_schroedinger)(void *m, myfloat128 (*f)(myfloat128))
{
    __float128 ff(__float128 x)
    {
        myfloat128 res, xx;
        F(xx) = x;
        res = f(xx);
        return F(res);
    }
    W(c_set_potential_schroedinger)(m, ff);
}

#elif(_DIM_==2)
void W(c_set_potential_schroedinger)(void *m, __float128 (*f)(__float128, __float128));
void S(c_set_potential_schroedinger)(void *m, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 y)
    {
        myfloat128 res, xx, yy;
        F(xx) = x;
        F(yy) = y;
        res = f(xx, yy);
        return F(res);
    }
    W(c_set_potential_schroedinger)(m, ff);
}    

#elif(_DIM_==3)
void W(c_set_potential_schroedinger)(void *m, __float128 (*f)(__float128, 
                                __float128, __float128));
void S(c_set_potential_schroedinger)(void *m, myfloat128 (*f)(myfloat128, 
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
    W(c_set_potential_schroedinger)(m, ff);
}    

#endif

void* W(c_get_potential_schroedinger)(void *m, int dim[_DIM_]);
void* S(c_get_potential_schroedinger)(void *m, int dim[_DIM_])
{
    return W(c_get_potential_schroedinger)(m, dim);
}

void W(c_load_potential_schroedinger)(void *m, char* filename, int filename_length);
void S(c_load_potential_schroedinger)(void *m, char* filename, int filename_length)
{
    W(c_load_potential_schroedinger)(m, filename, filename_length);
} 

void W(c_save_potential_schroedinger)(void *m, char* filename, int filename_length);
void S(c_save_potential_schroedinger)(void *m, char* filename, int filename_length)
{
    W(c_save_potential_schroedinger)(m, filename, filename_length);
}

void W(c_imaginary_time_propagate_A_wf_schroedinger)(void *psi, _COMPLEX_OR_REAL_ dt);
void S(c_imaginary_time_propagate_A_wf_schroedinger)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ dt)
{
    W(c_imaginary_time_propagate_A_wf_schroedinger)(psi, F(dt));
}    

void W(c_imaginary_time_propagate_B_wf_schroedinger)(void *psi, _COMPLEX_OR_REAL_ dt);
void S(c_imaginary_time_propagate_B_wf_schroedinger)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ dt)
{
    W(c_imaginary_time_propagate_B_wf_schroedinger)(psi, F(dt));
}    

void W(c_imaginary_time_propagate_C_wf_schroedinger)(void *psi, _COMPLEX_OR_REAL_ dt);
void S(c_imaginary_time_propagate_C_wf_schroedinger)(void *psi, _WRAPPED_COMPLEX_OR_REAL_ dt)
{
    W(c_imaginary_time_propagate_C_wf_schroedinger)(psi, F(dt));
}    

void W(c_add_apply_B_wf_schroedinger)(void *this, void *other,
                  _COMPLEX_OR_REAL_ coefficient);
void S(c_add_apply_B_wf_schroedinger)(void *this, void *other,
                  _WRAPPED_COMPLEX_OR_REAL_ coefficient)
{
    W(c_add_apply_B_wf_schroedinger)(this, other, F(coefficient));
}

__float128 W(c_kinetic_energy_wf_schroedinger)(void *psi);
myfloat128 S(c_kinetic_energy_wf_schroedinger)(void *psi)
{
    myfloat128 res;
    F(res) = W(c_kinetic_energy_wf_schroedinger)(psi);
    return res;
} 

__float128 W(c_potential_energy_wf_schroedinger)(void *psi);
myfloat128 S(c_potential_energy_wf_schroedinger)(void *psi)
{
    myfloat128 res;
    F(res) = W(c_potential_energy_wf_schroedinger)(psi);
    return res;
}

__float128 W(c_interaction_energy_wf_schroedinger)(void *psi);
myfloat128 S(c_interaction_energy_wf_schroedinger)(void *psi)
{
    myfloat128 res;
    F(res) = W(c_interaction_energy_wf_schroedinger)(psi);
    return res;
}


#if(_DIM_==1)
__float128 W(c_observable_wf_schroedinger)(void *psi, __float128 (*f)(__float128));
myfloat128 S(c_observable_wf_schroedinger)(void *psi, myfloat128 (*f)(myfloat128))
{
    __float128 ff(__float128 x)
    {
        myfloat128 res, xx;
        F(xx) = x;
        res = f(xx);
        return F(res);
    }
    myfloat128 res;
    F(res) = W(c_observable_wf_schroedinger)(psi, ff);
    return res;
}
#elif(_DIM_==2)
__float128 W(c_observable_wf_schroedinger)(void *psi, __float128 (*f)(__float128, __float128));
myfloat128 S(c_observable_wf_schroedinger)(void *psi, myfloat128 (*f)(myfloat128, myfloat128))
{
    __float128 ff(__float128 x, __float128 y)
    {
        myfloat128 res, xx, yy;
        F(xx) = x;
        F(yy) = y;
        res = f(xx, yy);
        return F(res);
    }
    myfloat128 res;
    F(res) = W(c_observable_wf_schroedinger)(psi, ff);
    return res;
}

#elif(_DIM_==3)
__float128 W(c_observable_wf_schroedinger)(void *psi, __float128 (*f)(__float128, 
                                __float128, __float128));
myfloat128 S(c_observable_wf_schroedinger)(void *psi, myfloat128 (*f)(myfloat128, 
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
    myfloat128 res;
    F(res) = W(c_observable_wf_schroedinger)(psi, ff);
    return res;
}    

#endif

void W(c_get_energy_expectation_deviation_wf_schroedinger)(void *psi, int dim[2]);
void S(c_get_energy_expectation_deviation_wf_schroedinger)(void *psi, int dim[2])
{
   W(c_get_energy_expectation_deviation_wf_schroedinger)(psi, dim);
}

void W(c_get_realspace_observables_wf_schroedinger)(void *psi, int dim[2+2*_DIM_]);
void S(c_get_realspace_observables_wf_schroedinger)(void *psi, int dim[2+2*_DIM_])
{
   W(c_get_realspace_observables_wf_schroedinger)(psi, dim);
}


void W(c_selfconsistent_nonlinear_step_wf_schroedinger)(void *psi, 
           _COMPLEX_OR_REAL_ dt, _COMPLEX_OR_REAL_ dt1, 
           __float128 eps, int max_iters);
void S(c_selfconsistent_nonlinear_step_wf_schroedinger)(void *psi, 
           _WRAPPED_COMPLEX_OR_REAL_ dt, _COMPLEX_OR_REAL_ dt1, 
           myfloat128 eps, int max_iters)
{
    S(c_selfconsistent_nonlinear_step_wf_schroedinger)(psi, dt, dt1, eps, max_iters);
}

#endif  /* QUADPRECISION */
