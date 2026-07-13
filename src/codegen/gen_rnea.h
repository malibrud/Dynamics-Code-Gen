#pragma once

#include <stdio.h>
#include "system.h"
#include "prln.h"
#include "svec.h"
#include "inertia.h"

// Computes the joint-space bias force C = ID(model, q, qd, 0, fext) via the
// RNEA outward pass (velocity, acceleration, force) plus inward accumulation
// (RBDA Table 5.1, run with aJ = 0 -- see CRBA.pdf slide 3: H*qdd + C = tau =
// ID(model,q,qd,qdd,fext), so C = ID(model,q,qd,0,fext)). External forces
// are applied only to bodies whose sys->hasFext[i] flag is set; others cost
// nothing. Must run after gen_fwd_dyn_crba() has populated
// sys->xFromParent[]/xFromBase[]/S[]/vJ[], and sys->bodies[]/I[]/grav have
// been filled in by the caller.
void gen_rnea_C( FILE *f, System *sys ) {
    prln_file = f;
    int NB = sys->nBod;

    prln( "" );
    prln( "// Bias force C = ID(model, q, qd, 0).  Table 5.1 in RBDA, outward pass." );
    prln( "" );

    for ( int i = 0 ; i < NB ; i++ ) {
        int iParent = sys->parent[ i ];
        auto X6 = xformTo6x6( sys->xFromParent[ i ] );

        // v_i = i_X_lambda(i) * v_lambda(i) + vJ_i  (v_lambda(i) = 0 at the base)
        if ( iParent == -1 ) {
            sys->v[ i ] = sys->vJ[ i ];
        } else {
            Svec vParentInI = smat6MulVec( X6, sys->v[ iParent ] );
            sys->v[ i ] = svecAdd( "", vParentInI, sys->vJ[ i ] );
        }
        sys->v[ i ].baseName = format( "v{}_", i );
        sys->v[ i ].emit();

        // a_i = i_X_lambda(i) * a_lambda(i) + v_i x vJ_i  (aJ = 0; a_lambda(0) = -a_grav)
        Svec aParent = smat6MulVec( X6, iParent == -1 ? sys->grav : sys->a[ iParent ] );
        Svec viCrossVJ = svCrossMotion( "", sys->v[ i ], sys->vJ[ i ] );
        sys->a[ i ] = svecAdd( format( "a{}_", i ), aParent, viCrossVJ );
        sys->a[ i ].emit();

        // f_i = I_i * a_i + v_i x* (I_i * v_i) - i_X*_0 * fext_i
        Svec Ia = spatialInertiaMulVec( format( "Ia{}_", i ), sys->I[ i ], sys->a[ i ] );
        Ia.emit();
        Svec Iv = spatialInertiaMulVec( format( "Iv{}_", i ), sys->I[ i ], sys->v[ i ] );
        Iv.emit();
        Svec vIv = svCrossForce( format( "vIv{}_", i ), sys->v[ i ], Iv );
        vIv.emit();

        Svec fBias = svecAdd( format( "f{}_0_", i ), Ia, vIv );

        if ( sys->hasFext[ i ] ) {
            // fext_i arrives in absolute coordinates and must move to body-i
            // coordinates: i_X*_0 * fext = (0_X_i)^T * fext, where 0_X_i is
            // the inverse of our xFromBase[i] (spatial_v2's
            // apply_external_forces does the same via Xa{i}' \ f_ext{i}).
            Svec fx( format( "fext{}_", i ) );
            for ( int k = 0 ; k < 6 ; k++ ) {
                prln( "double fext%d_%d = fext[ %d ];", i, k, i*6 + k );
                fx.v[ k ] = format( "fext{}_{}", i, k );
            }

            Xform xToBase = invertXform( format( "Xinv{}", i ), sys->xFromBase[ i ] );
            xToBase.emit();
            auto Xinv6  = xformTo6x6( xToBase );
            auto Xstar6 = smat6Transpose( Xinv6 );
            Svec fextBody = smat6MulVec( Xstar6, fx );
            fextBody.baseName = format( "fextB{}_", i );
            fextBody.emit();

            fBias = svecSub( format( "f{}_0_", i ), fBias, fextBody );
        }

        sys->force[ i ] = fBias;
        sys->force[ i ].emit();
        prln( "" );
    }

    for ( int i = NB - 1 ; i >= 0 ; i-- ) {
        prln( "double C%d = %s;", i, sys->force[ i ].v[2].c_str() );

        int iParent = sys->parent[ i ];
        if ( iParent != -1 ) {
            // f_lambda(i) += lambda(i)_X*_i * f_i, i.e. X_i^T * f_i (same
            // transpose trick as the Ic/H propagation in gen_crba.h).
            auto X6  = xformTo6x6( sys->xFromParent[ i ] );
            auto Xt6 = smat6Transpose( X6 );
            Svec childForce = smat6MulVec( Xt6, sys->force[ i ] );
            childForce.baseName = format( "fUp{}_", i );
            childForce.emit();

            sys->force[ iParent ] = svecAdd( format( "f{}_via{}_", iParent, i ), sys->force[ iParent ], childForce );
            sys->force[ iParent ].emit();
        }
    }
    prln( "" );
}
