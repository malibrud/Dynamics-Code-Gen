#pragma once
#include <inttypes.h>
#include "xform.h"

#define J_REVOLUTE 0
#define J_PRISMATIC 1

typedef struct {
    uint8_t type;
    uint8_t dof;
    Xform *xtree;
    Xform *xjoint;
    double *q;
    double *qd;
    double *qdd;
} Joint;

void j_updateXform( Joint *j ) {
    switch ( j->type ) {
        case J_REVOLUTE:
            x_update( j->xjoint, j->q );
            break;
        case J_PRISMATIC:
            break;
    }
}