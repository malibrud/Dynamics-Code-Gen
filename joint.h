#pragma once
#include <inttypes.h>

#define J_REVOLUTE 0
#define J_PRISMATIC 1

typedef struct {
    uint8_t type;
    uint8_t dof;
    double *q;
    double *qd;
    double *qdd;
} Joint;