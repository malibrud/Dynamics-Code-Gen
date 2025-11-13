#pragma once

#include "system.h"
#include <stdlib.h>

void fwd_dyn_crba( System *sys, double *q, double *qd, double *tau, double *out_qdd ) {
    int NB = sys->nBod;

    // This section computes C = ID( model, q, qd, 0, f ) .. See eq (6.2) in RBDA, pp 103
    // The algorithm is in Table 5.1, RBDA pp 96.
    for ( int i = 0 ; i < NB ; i++ ) {
        int iParent = sys->parent[ i ];
        // Calculate the joint values, i.e. jcalc() in RBDA.
        Svec vJ = {};
        Svec aJ = {};
        if ( sys->jtype[ i ] == J_REVOLUTE ) {
            // Get X_J from q.
            x_update( sys->xJoint[ i ], &q[ i ] );

            // Get S_i.
            sys->S.v[ 2 ] = 1;

            // Get v_J = S_i * qdd_i.
            vJ.v[ 2 ] = qd[ i ];

            // Get a_J = S_i * qdd_i.
            aJ.v[ 2 ] = qd[ i ];
        } else {
            assert( false );
        }
        
        // This section computes the parent to body transform and the base to body transforms.
        // See RBDA Table 5.1.
        // Note that the body number scheme differs from RBDA in that the first link is numbered 0.
        // The base is then denoted as -1 = \lambda(0)

        // i_X_\lambda(i) = X_J * X_T(i)
        x_compose( sys->xTree[ i ], sys->xJoint[ i ], sys->xFromParent[ i ] );

        // i_X_\lambda(0) = i_X_\lambda(i) * lambda(i)_X_\lambda(0)
        Xform *xFromBaseParent = ( i == 0 ) 
            ? &x_identity 
            : sys->xFromBase[ iParent ];
        x_compose( xFromBaseParent, sys->xFromParent[ i ], sys->xFromBase[ i ] );

        // This section computes the velocity of the body in body coordinates.
        // v_i = i_X_L(i) * v_L(i) + v_J
        if ( i == 0 ) {
            sys->vel[ i ] = vJ;
        } else {
            x_mapply( sys->xFromParent[ i ], &sys->vel[ iParent ], &sys->vel[ i ] );
            sv_addTo( &sys->vel[ i ], &vJ );
        }

        // This section computes the acceleration of the body in body coordinates.
        // a_i = i_X_L(i) * a_L(i) + v_i x v_J
        // Note: Because we are computing C, we assume aJ = 0.
        if ( i == 0 ) {
            x_mapply( sys->xFromParent[ i ], &sys->grav, &sys->acc[ i ] );
        } else {
            x_mapply( sys->xFromParent[ i ], &sys->acc[ iParent ], &sys->acc[ i ] );
        }
        Svec vi_cross_vJ;
        sv_mcross( &sys->vel[ i ], &vJ, &vi_cross_vJ );
        sv_addTo( &sys->acc[ i ], &vi_cross_vJ );

        // This section computes the force on the body in body coordinates.
        // f_i = I_i * a_i + v_i x* I_i * v_i - i_X*_0 * f_ext_i
        // Note: Because we are computing C, we assume aJ = 0.
        Svec *fi = &sys->force[ i ]; // I0_i * a_i

        // first term
        bd_m2f( &sys->I0[ i ], &sys->acc[ i ], fi );

        // second term
        Svec Iv; // Iv = I0_i * v_i
        bd_m2f( &sys->I0[ i ], &sys->vel[ i ], &Iv );
        Svec vIv; // vi x* Iv
        sv_fcross( &sys->vel[ i ], &Iv, &vIv );
        sv_addTo( fi, &vIv );

        // third term
        Svec fext;
        x_fapply( sys->xFromBase[ i ], &sys->fext[ iParent ], &fext );
        sv_subFrom( fi, &fext );
    }

    for ( int i = NB-1 ; i >= 0 ; i-- ) {
        // Get the joint force as the element of C, i.e. eq (6.2).
        // C_i = S_i^T * f_i
        if ( sys->jtype[ i ] == J_REVOLUTE ) {
            sys->C[ i ] = sys->force[ i ].v[ 2 ];
        } else { 
            assert( false );
        }

        // Accumulate forces back up the chain from child to parent.
        // f_L(i) += L(i)_X_i * f_i
        int iParent = sys->parent[ i ];
        Svec childForce;
        x_finvapply( sys->xFromParent[ i ], &sys->force[ i ], &childForce );
        sv_addTo( &sys->force[ iParent ] , &childForce );
    }

    // Calculate H
    // See RBDA, Table 6.2, pp 107

    // H = 0
    memset( sys->H, 0, NB*NB*sizeof(double) );

    for ( int i = 0 ; i < NB ; i++ ) {
        sys->Ic[i] = sys->I0[i];
    }

    // Generate composite rigid body inertias.
    // Ic_i = I0_i + Sum(j = mu(i)) { i_X*_j * Ic_j * j_X_i }
    for ( int i = 0 ; i < NB ; i++ ) sys->Ic[ i ] = sys->I0[ i ];
    for ( int i = NB-1 ; i >= 0 ; i-- ) {
        int iParent = sys->parent[ i ];
        if ( i != 0 ) {
            SpatialInertia Ip; // Ic in parent frame;
            x_Iinvapply( *sys->xFromParent, &sys->Ic[ i ], &Ip );
            si_addTo( &sys->Ic[ iParent ], &Ip);
        }
        Svec F;
        si_multBySvec( &sys->Ic[ i ], &sys->S[ i ], &F );
        sys->H[ i*NB + i ] = si_multforcebyjointsubspace( &F, &sys->S[ i ] );

        int j = i;
        while ( j != 0 ) {
            
        }
    }
}