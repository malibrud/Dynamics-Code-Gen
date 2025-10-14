#pragma once
#include "vec3.h"

typedef struct {
    double xx;
    double yy;
    double zz;
    double xy;
    double yz;
    double zx;
} Inertia;

typedef struct {
    double m;
    double mcx;
    double mcy;
    double mcz;
    double mcxcy;
    double mcycz;
    double mczcx;
    double Ixx;
    double Iyy;
    double Izz;
    double Ixy;
    double Iyz;
    double Izx;
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

    I0.Ixy = B.I.xy + m*x*y;
    I0.Iyz = B.I.yz + m*y*z;
    I0.Izx = B.I.zx + m*z*x;

    I0.mcx = m * x;
    I0.mcy = m * y;
    I0.mcz = m * z;

    I0.mcxcy = m * x * y;
    I0.mcycz = m * y * z;
    I0.mczcx = m * z * x;

    return I0;
}

// Computes vOut = I * vIn
void bd_m2f( SpatialInertia *I, Svec *vIn, Svec *vOut ) {
    double *f = vOut->v;
    double *v = vIn->v;
    double m = I->m;
    double mx = I->mcx;
    double my = I->mcy;
    double mz = I->mcz;
    double mxy = I->mcxcy;
    double myz = I->mcycz;
    double mzx = I->mczcx;
    double Ixx = I->Ixx;
    double Iyy = I->Iyy;
    double Izz = I->Izz;
    double Ixy = I->Ixy;
    double Iyz = I->Iyz;
    double Izx = I->Izx;

    f[ 0 ] =  Ixx*v[0] - Ixy*v[1] - Izx*v[2]   +    0    - mz*v[4] + my*v[5];
    f[ 1 ] = -Ixy*v[0] + Iyy*v[1] - Iyz*v[2]   + mz*v[3] +    0    - mx*v[5];
    f[ 2 ] = -Izx*v[0] - Iyz*v[1] + Izz*v[2]   - my*v[3] + mx*v[4] +    0   ;

    f[ 3 ] =    0*v[0] +  mz*v[1] -  my*v[2]   +  m*v[3] +    0    +    0   ;
    f[ 4 ] =  -mz*v[0] +   0*v[1] +  mx*v[2]   +    0    +  m*v[4] +    0   ;
    f[ 5 ] =   my*v[0] -  mx*v[1] +   0*v[2]   +    0    +    0    +  m*v[5];
}
