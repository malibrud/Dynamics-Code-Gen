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

    if ( sys->loops.empty() ) {
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
        return;
    }

    // Closed loops: method 2 of RBDA sec 8.5. One elimination of H serves
    // both the tree solve qdd0 = H^-1 (tau - C) and Y = H^-1 K^T; then the
    // nc x nc system (K Y) lambda = kful - K qdd0 gives the constraint
    // forces, and qdd = qdd0 + Y lambda. K's structural zeros (tree joints
    // not on any loop path, plus mechanism-specific cancellations) are baked
    // as 0.0 literals.
    int NC = sys->nc;
    int M = NC + 1;

    prln( "" );
    prln( "// Loop solve: B = [ tau - C | K^T ], then H^-1 applied to all columns." );
    prln( "double B[ %d ] = {", NV * M );
    prln_indent_level++;
    for ( int r = 0 ; r < NV ; r++ ) {
        string row = format( "tau[ {} ] - C{}, ", r, r );
        for ( int c = 0 ; c < NC ; c++ ) {
            const string &e = sys->Kentry[ c ][ r ];
            row += ( e == "0" ? "0.0"s : e ) + ", ";
        }
        prln( "%s", row.c_str() );
    }
    prln_indent_level--;
    prln( "};" );
    prln( "solveLinearSystemM( %d, H, %d, B );", NV, M );

    prln( "double Km[ %d ] = {", NC * NV );
    prln_indent_level++;
    for ( int r = 0 ; r < NC ; r++ ) {
        string row;
        for ( int c = 0 ; c < NV ; c++ ) {
            const string &e = sys->Kentry[ r ][ c ];
            row += ( e == "0" ? "0.0"s : e ) + ", ";
        }
        prln( "%s", row.c_str() );
    }
    prln_indent_level--;
    prln( "};" );

    prln( "double kful[ %d ] = {", NC );
    prln_indent_level++;
    for ( int r = 0 ; r < NC ; r++ ) prln( "%s,", sys->kfulName[ r ].c_str() );
    prln_indent_level--;
    prln( "};" );

    prln( "// Ac = K Y, bc = kful - K qdd0 (qdd0 is column 0 of B, Y the rest)." );
    prln( "double Ac[ %d ], bc[ %d ], lambda[ %d ];", NC * NC, NC, NC );
    prln( "for ( int r = 0 ; r < %d ; r++ ) {", NC );
    prln_indent_level++;
    prln( "bc[ r ] = kful[ r ];" );
    prln( "for ( int c = 0 ; c < %d ; c++ ) Ac[ r*%d + c ] = 0.0;", NC, NC );
    prln( "for ( int k = 0 ; k < %d ; k++ ) {", NV );
    prln_indent_level++;
    prln( "double krk = Km[ r*%d + k ];", NV );
    prln( "if ( krk == 0.0 ) continue;" );
    prln( "bc[ r ] -= krk * B[ k*%d ];", M );
    prln( "for ( int c = 0 ; c < %d ; c++ ) Ac[ r*%d + c ] += krk * B[ k*%d + 1 + c ];", NC, NC, M );
    prln_indent_level--;
    prln( "}" );
    prln_indent_level--;
    prln( "}" );
    prln( "solveLinearSystem( %d, Ac, bc, lambda );", NC );

    prln( "for ( int r = 0 ; r < %d ; r++ ) {", NV );
    prln_indent_level++;
    prln( "double s = B[ r*%d ];", M );
    prln( "for ( int c = 0 ; c < %d ; c++ ) s += B[ r*%d + 1 + c ] * lambda[ c ];", NC, M );
    prln( "out_qdd[ r ] = s;" );
    prln_indent_level--;
    prln( "}" );

    prln_indent_level--;
    prln( "}" );
}
