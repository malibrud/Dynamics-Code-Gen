#pragma once

#include "system.h"

void inv_dyn_rnea( System *sys, double *q, double qd, double *qdd, double *out_tau ) {
    int NB = sys->nBod;

    // update transforms
    for ( int i = 1 ; i < NB ; i++ ) {
        Joint *j = &sys->joints[ i ];
        int iParent = sys->parent[ i ];
        j_updateXform( j );
        x_compose( j->xtree, j->xjoint, sys->xFromParent[ i ] );
        x_compose( sys->xFromBase[ i-1 ], sys->xFromParent[ i ], sys->xFromBase[ i ] );
    }
}