#pragma once

#include <stdio.h>
#include "system.h"
#include "prln.h"
#include "inertia.h"

// Computes each body's own spatial inertia I_i (about its own frame) from
// sys->bodies[], populating sys->I[]. Both gen_rnea_C() and gen_crba_H() need
// this, so it runs once, before either. Must run after gen_fwd_dyn_crba() and
// after sys->bodies[] has been filled in by the caller.
void gen_body_inertias( FILE *f, System *sys ) {
    prln_file = f;
    int NB = sys->nBod;

    prln( "" );
    prln( "// Each body's own spatial inertia, about its own frame." );
    prln( "" );
    for ( int i = 0 ; i < NB ; i++ ) {
        sys->I[ i ] = bodySpatialInertia( format( "I{}", i ), sys->bodies[ i ] );
        sys->I[ i ].emit();
    }
    prln( "" );
}

// Computes the joint-space inertia matrix H via the composite-rigid-body
// algorithm (RBDA Table 6.4, pseudocode on slide 8 of CRBA.pdf). Must run
// after gen_body_inertias() has populated sys->I[], and gen_fwd_dyn_crba()
// has populated sys->xFromParent[]/S[]. Appends to -- but does not close --
// the C function gen_fwd_dyn_crba() opened; gen_solve_qdd() closes it.
void gen_crba_H( FILE *f, System *sys ) {
    prln_file = f;
    int NB = sys->nBod;

    prln( "" );
    prln( "// Composite rigid body algorithm: H = CRBA(model, q).  Table 6.4 in RBDA." );
    prln( "" );

    // Ic_i = I_i for every body -- start each composite inertia as just the
    // body's own inertia; children get folded in below as we walk inward.
    for ( int i = 0 ; i < NB ; i++ ) sys->Ic[ i ] = sys->I[ i ];

    for ( int i = NB - 1 ; i >= 0 ; i-- ) {
        int iParent = sys->parent[ i ];

        if ( iParent != -1 ) {
            prln( "// Fold body %d's composite inertia into its parent, body %d.", i, iParent );
            SpatialInertia contribution = transformInertia( format( "IcUp{}", i ), sys->xFromParent[ i ], sys->Ic[ i ] );
            contribution.emit();

            // Named after the child (i), not just the parent, so that a
            // parent with more than one child gets a fresh temp for each
            // fold instead of colliding on the same variable name.
            SpatialInertia updated = addInertia( format( "Ic{}_via{}", iParent, i ), sys->Ic[ iParent ], contribution );
            updated.emit();
            sys->Ic[ iParent ] = updated;
            prln( "" );
        }

        Svec F = spatialInertiaMulVec( format( "F{}_0_", i ), sys->Ic[ i ], sys->S[ i ] );
        F.emit();
        string Hii = svecDot( F, sys->S[ i ] );
        prln( "double H%d%d = %s;", i, i, Hii.c_str() );
        prln( "" );

        // Walk up the ancestor chain, propagating F and picking off the
        // off-diagonal H entries. For a simple chain (no branches) this
        // walks all the way to the root every time; branch-induced sparsity
        // only shows up -- as permanently-zero H entries between bodies on
        // separate branches -- once the tree actually branches. Each step
        // gets a unique F{i}_{step}_ temp name so nothing gets redeclared.
        int j = i;
        int step = 0;
        while ( sys->parent[ j ] != -1 ) {
            auto X6  = xformTo6x6( sys->xFromParent[ j ] );
            auto Xt6 = smat6Transpose( X6 );
            F = smat6MulVec( Xt6, F );
            step++;
            F.baseName = format( "F{}_{}_", i, step );
            F.emit();

            j = sys->parent[ j ];

            string Hij = svecDot( F, sys->S[ j ] );
            prln( "double H%d%d = %s;", i, j, Hij.c_str() );
            prln( "double H%d%d = H%d%d;", j, i, i, j );
            prln( "" );
        }
    }
}
