#pragma once

#include <stdio.h>
#include "system.h"
#include "prln.h"

// True iff the CRBA emitted an H(i,j) temp: the ancestor walk only visits
// pairs where one body supports the other. Pairs on separate branches are
// branch-induced structural zeros (RBDA sec 6.4 / CRBA.pdf slide 10) and no
// temp exists for them.
static bool crbaHasEntry( System *sys, int i, int j ) {
    if ( i == j ) return true;
    for ( int k = sys->parent[ i ] ; k != -1 ; k = sys->parent[ k ] ) if ( k == j ) return true;
    for ( int k = sys->parent[ j ] ; k != -1 ; k = sys->parent[ k ] ) if ( k == i ) return true;
    return false;
}

// Step 3 of the RBDA forward-dynamics recipe (CRBA.pdf slide 2): solve
// H*qdd = tau - C for qdd, using the H_r_c/C(k) temps gen_crba_H() and
// gen_rnea_C() already emitted (indexed by velocity coordinate, since a
// floating joint occupies a 6-wide block), and write the result to out_qdd.
// Closes the C function gen_fwd_dyn_crba() opened.
void gen_solve_qdd( FILE *f, System *sys ) {
    prln_file = f;
    int NB = sys->nBod;
    int NV = sys->nv;

    // Which body owns each velocity coordinate, for the structural-zero test.
    vector<int> bodyOf( NV );
    for ( int i = 0 ; i < NB ; i++ ) {
        int w = sys->jtype[ i ] == J_FLOATING ? 6 : 1;
        for ( int k = 0 ; k < w ; k++ ) bodyOf[ sys->vidx[ i ] + k ] = i;
    }

    prln( "" );
    prln( "// Solve H*qdd = tau - C for qdd. Entries written as 0.0 are" );
    prln( "// branch-induced structural zeros: those body pairs share no" );
    prln( "// common subtree, so the CRBA never computes them." );
    prln( "double H[ %d ] = {", NV * NV );
    prln_indent_level++;
    for ( int r = 0 ; r < NV ; r++ ) {
        string row;
        for ( int c = 0 ; c < NV ; c++ )
            row += crbaHasEntry( sys, bodyOf[ r ], bodyOf[ c ] ) ? format( "H_{}_{}, ", r, c ) : "0.0, "s;
        prln( "%s", row.c_str() );
    }
    prln_indent_level--;
    prln( "};" );

    prln( "double rhs[ %d ] = {", NV );
    prln_indent_level++;
    for ( int k = 0 ; k < NV ; k++ ) prln( "tau[ %d ] - C%d,", k, k );
    prln_indent_level--;
    prln( "};" );

    prln( "double qdd[ %d ];", NV );
    prln( "solveLinearSystem( %d, H, rhs, qdd );", NV );
    for ( int k = 0 ; k < NV ; k++ ) prln( "out_qdd[ %d ] = qdd[ %d ];", k, k );

    prln_indent_level--;
    prln( "}" );
}
