#pragma once
#include <assert.h>
#include "vec3.h"

typedef struct {
    double xx;
    double yy;
    double zz;
    double xy;
    double yz;
    double zx;
} Inertia;

/* 
Inertia tensor is defined with the off-diagonal minus signs absorbed into the terms.

      [ Ixx Ixy Izx ]
I_0 = [ Ixy Iyy Iyz ]
      [ Izx Iyz Izz ]
*/
typedef struct {
    double m;
    double sx; // m * c_x
    double sy; // m * c_y
    double sz; // m * c_z
    double Ixx;
    double Iyy;
    double Izz;
    double Ixy; // minus sign absorbed into term.
    double Iyz; // minus sign absorbed into term.
    double Izx; // minus sign absorbed into term.
} SpatialInertia;

typedef struct {
    double m;
    double cg[3];
    Inertia I;
} Body;

SpatialInertia bd_translateInertia( Body B, Vec3 v, double m ) {
    SpatialInertia I0 = { 0 };
    double x = v.x;
    double y = v.y;
    double z = v.z;

    I0.m = m;
    I0.Ixx = B.I.xx + m*( z*z + y*y );
    I0.Iyy = B.I.yy + m*( x*x + z*z );
    I0.Izz = B.I.zz + m*( y*y + x*x );

    I0.Ixy = B.I.xy - m*x*y;
    I0.Iyz = B.I.yz - m*y*z;
    I0.Izx = B.I.zx - m*z*x;

    I0.sx = m * x;
    I0.sy = m * y;
    I0.sz = m * z;

    return I0;
}

// Computes vOut = I * vIn
void bd_m2f( SpatialInertia *I, Svec *vIn, Svec *vOut ) {
    double *f = vOut->v;
    double *v = vIn->v;
    double wx = v[0];   // Omega_x
    double wy = v[1];   // Omega_y
    double wz = v[2];   // Omega_z
    double vx = v[0];   // v0_x
    double vy = v[1];   // v0_y
    double vz = v[2];   // v0_z
    double m = I->m;
    double mx = I->sx;
    double my = I->sy;
    double mz = I->sz;
    double Ixx = I->Ixx;
    double Iyy = I->Iyy;
    double Izz = I->Izz;
    double Ixy = I->Ixy;
    double Iyz = I->Iyz;
    double Izx = I->Izx;

    f[ 0 ] =  Ixx*wx + Ixy*wy + Izx*wz   +    0  - mz*vy + my*vz  ;
    f[ 1 ] =  Ixy*wx + Iyy*wy + Iyz*wz   + mz*vx +    0  - mx*vz  ;
    f[ 2 ] =  Izx*wx + Iyz*wy + Izz*wz   - my*vx + mx*vy +    0   ;

    f[ 3 ] =    0*wx +  mz*wy -  my*wz   +  m*vx +    0  +    0   ;
    f[ 4 ] =  -mz*wx +   0*wy +  mx*wz   +    0  +  m*vy +    0   ;
    f[ 5 ] =   my*wx -  mx*wy +   0*wz   +    0  +    0  +  m*vz  ;
}

void si_addTo( SpatialInertia *a, SpatialInertia *b ) {
    a->Ixx += b->Ixx;
    a->Iyy += b->Iyy;
    a->Izz += b->Izz;
    a->Ixy += b->Ixy;
    a->Iyz += b->Iyz;
    a->Izx += b->Izx;
    a->m   += b->m;
    a->sx  += b->sx;
    a->sy  += b->sy;
    a->sz  += b->sz;
}

// Computes vOut = I * vIn
void si_multBySvec( SpatialInertia *I, Svec *vIn, Svec *vOut ) {
    double *vi = vIn->v;
    double *vo = vOut->v;

    // Revolute Joint
    if ( vi[ 2 ] == 1 ) {
        assert( vi[ 0 ] == 0 );
        assert( vi[ 1 ] == 0 );
        assert( vi[ 3 ] == 0 );
        assert( vi[ 4 ] == 0 );
        assert( vi[ 5 ] == 0 );

        // Linear part is m*(w x c) with w = z: (-sy, sx, 0). These signs
        // were originally flipped (+sy, -sx) -- caught by src/validate after
        // the same bug was copied into the code generator (cf. bd_m2f above,
        // which always had them right).
        vo[ 0 ] = I->Izx;
        vo[ 1 ] = I->Iyz;
        vo[ 2 ] = I->Izz;
        vo[ 3 ] = -I->sy;
        vo[ 4 ] = I->sx;
        vo[ 5 ] = 0;
        return;
    }
    assert( false );
}

// Computes fOut = fIn^T * S
double si_multForceByJointSubspace( Svec *fIn, Svec *S ) {
    double *f = fIn->v;
    double *s = S->v;

    // Revolute Joint
    if ( s[ 2 ] == 1 ) {
        assert( s[ 0 ] == 0 );
        assert( s[ 1 ] == 0 );
        assert( s[ 3 ] == 0 );
        assert( s[ 4 ] == 0 );
        assert( s[ 5 ] == 0 );

        return f[ 2 ];
    }
    assert( false );
}