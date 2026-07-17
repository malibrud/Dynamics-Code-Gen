#pragma once
#include <inttypes.h>

#define J_REVOLUTE 0
#define J_PRISMATIC 1
#define J_FLOATING 2   // genuine 6-DOF joint: quaternion orientation, S = identity(6)

// Loop-joint types (RBDA ch. 8): these close kinematic loops and impose
// constraints instead of adding coordinates, so they never appear as tree
// joints. nc = 6 minus the joint's motion freedoms; see loopNc() in
// gen_loops.h.
#define J_CYLINDRICAL 3
#define J_SPHERICAL   4
#define J_POINTONLINE 5  // sphere-in-cylinder: 2 constraints, closes planar loops without redundancy
#define J_DISTANCE    6  // rigid rod between two anchor points: 1 constraint

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