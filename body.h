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

SpatialInertia bd_translateInertia( Body I, Vec3 v, double m ) {
    Inertia I, Vec3 v, double m 
    SpatialInertia I0 = {};
    double x = v.x;
    double y = v.y;
    double z = v.z;

    I0.m = m;
    I0.Ixx = I.xx + m*( z*z + y*y );
    I0.Iyy = I.yy + m*( x*x + z*z );
    I0.Izz = I.zz + m*( y*y + x*x );

    I0.Ixy = I.xy - m*x*y;
    I0.Iyz = I.yz - m*y*z;
    I0.Izx = I.zx - m*z*x;

    I0.mcx = m * x;
    I0.mcy = m * y;
    I0.mcz = m * z;

    return I0;
}
