#pragma once

#include <stdio.h>
#include "system.h"
#include "prln.h"

// Step 3 of the RBDA forward-dynamics recipe (CRBA.pdf slide 2): solve
// H*qdd = tau - C for qdd, using the H(i,j)/C(i) temps gen_crba_H() and
// gen_rnea_C() already emitted, and write the result to out_qdd. Closes the
// C function gen_fwd_dyn_crba() opened.
void gen_solve_qdd( FILE *f, System *sys ) {
    prln_file = f;
    int NB = sys->nBod;

    prln( "" );
    prln( "// Solve H*qdd = tau - C for qdd." );
    prln( "double H[ %d ] = {", NB * NB );
    prln_indent_level++;
    for ( int i = 0 ; i < NB ; i++ ) {
        string row;
        for ( int j = 0 ; j < NB ; j++ ) row += format( "H{}{}, ", i, j );
        prln( "%s", row.c_str() );
    }
    prln_indent_level--;
    prln( "};" );

    prln( "double rhs[ %d ] = {", NB );
    prln_indent_level++;
    for ( int i = 0 ; i < NB ; i++ ) prln( "tau[ %d ] - C%d,", i, i );
    prln_indent_level--;
    prln( "};" );

    prln( "double qdd[ %d ];", NB );
    prln( "solveLinearSystem( %d, H, rhs, qdd );", NB );
    for ( int i = 0 ; i < NB ; i++ ) prln( "out_qdd[ %d ] = qdd[ %d ];", i, i );

    prln_indent_level--;
    prln( "}" );
}
