#pragma once

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
