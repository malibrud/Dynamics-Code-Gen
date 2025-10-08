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

    I0.Ixy = B.I.xy - m*x*y;
    I0.Iyz = B.I.yz - m*y*z;
    I0.Izx = B.I.zx - m*z*x;

    I0.mcx = m * x;
    I0.mcy = m * y;
    I0.mcz = m * z;

    return I0;
}
