#pragma once

#include "system.h"

void fwd_dyn_crba( System *sys, double *q, double *qd, double *tau, double *out_qdd ) {
    int NB = sys->nBod;

    // This section computes C = ID( model, q, qd, 0, f ) .. See eq (6.2) in RBDA, pp 103
    // The algorithm is in Table 5.1, RBDA pp 96.
    for ( int i = 0 ; i < NB ; i++ ) {
        int iParent = sys->parent[ i ];
        // Calculate the joint values, i.e. jcalc() in RBDA.
        Svec S = {};
        Svec vJ = {};
        if ( sys->jtype[ i ] == J_REVOLUTE ) {
            // Get X_J from q.
            x_update( sys->xJoint[ i ], &q[ i ] );

            // Get S_i.
            S.v[ 2 ] = 1;

            // Get v_J.
            vJ.v[ 2 ] = qd[ i ];
        } else {
            assert( false );
        }
        
        // Compute 
        x_compose( sys->xTree[ i ], sys->xJoint[ i ], sys->xFromParent[ i ] );
        Xform *xFromBaseParent = ( i == 0 ) 
            ? &x_identity 
            : sys->xFromBase[ iParent ];
        x_compose( xFromBaseParent, sys->xFromParent[ i ], sys->xFromBase[ i ] );
    }
}