#pragma once
#include <inttypes.h>

#define J_REVOLUTE 0
#define J_PRISMATIC 1
#define J_FLOATING 2   // genuine 6-DOF joint: quaternion orientation, S = identity(6)

// typedef struct {
//     uint8_t type;
//     uint8_t dof;
//     Xform *xtree;
//     Xform *xjoint;
//     double *q;
//     double *qd;
//     double *qdd;
// } Joint;
// 
// void j_updateXform( Joint *j ) {
//     switch ( j->type ) {
//         case J_REVOLUTE:
//             x_update( j->xjoint, j->q );
//             break;
//         case J_PRISMATIC:
//             break;
//     }
// }

// void j_calc( )