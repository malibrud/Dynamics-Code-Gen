#pragma once

#include <stdio.h>
#include <assert.h>
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

    // The expanded parent array (RBDA Table 6.4): lam[r] is the velocity
    // coordinate that precedes r on the path to the base. Within a floating
    // body's 6-wide block the coordinates form a chain; the block's first
    // coordinate points at the parent body's LAST one.
    vector<int> lam( NV );
    for ( int i = 0 ; i < NB ; i++ ) {
        int v0 = sys->vidx[ i ], w = sys->jtype[ i ] == J_FLOATING ? 6 : 1;
        int par = sys->parent[ i ];
        lam[ v0 ] = par == -1 ? -1
                  : sys->vidx[ par ] + ( sys->jtype[ par ] == J_FLOATING ? 6 : 1 ) - 1;
        for ( int k = 1 ; k < w ; k++ ) lam[ v0 + k ] = v0 + k - 1;
    }

    // lam must reproduce H's sparsity pattern exactly, or the factorization
    // would skip a nonzero: H_rc is nonzero iff one of r, c is reachable from
    // the other by walking lam. That is the same relation crbaHasEntry() tests
    // at body level, so the two must agree on every pair. Costs O(nv^2) here
    // and nothing at runtime. Also checks Featherstone's lam[r] < r, which the
    // factorization's reverse-order sweep depends on.
    for ( int r = 0 ; r < NV ; r++ ) {
        assert( lam[ r ] < r );
        for ( int c = 0 ; c < NV ; c++ ) {
            bool viaLam = false;
            for ( int k = r ; k != -1 ; k = lam[ k ] ) if ( k == c ) viaLam = true;
            for ( int k = c ; k != -1 ; k = lam[ k ] ) if ( k == r ) viaLam = true;
            assert( viaLam == crbaHasEntry( sys, bodyOf[ r ], bodyOf[ c ] ) );
        }
    }

    prln( "" );
    prln( "// Solve H*qdd = tau - C for qdd, via the sparse L^T D L factorization" );
    prln( "// of RBDA sec 6.5 (see ltdlFactor/ltdlSolveM above)." );
    prln( "//" );
    prln( "// lam is the expanded parent array over velocity coordinates. It is" );
    prln( "// also H's sparsity pattern: a pair of coordinates has a nonzero H" );
    prln( "// entry exactly when one is reachable from the other along lam, so" );
    prln( "// every entry NOT emitted below is a branch-induced structural zero" );
    prln( "// (those bodies share no common subtree) that the factorization" );
    prln( "// provably never touches." );
    prln( "static const int Hlam[ %d ] = {", NV );
    prln_indent_level++;
    {
        string row;
        for ( int r = 0 ; r < NV ; r++ ) row += format( "{}, ", lam[ r ] );
        prln( "%s", row.c_str() );
    }
    prln_indent_level--;
    prln( "};" );

    prln( "// Only the lower triangle is ever read -- the factorization and both" );
    prln( "// substitutions stay below the diagonal -- so H is left uninitialized" );
    prln( "// and only that half is written." );
    prln( "double H[ %d ];", NV * NV );
    for ( int r = 0 ; r < NV ; r++ ) {
        string row;
        for ( int c = 0 ; c <= r ; c++ )
            if ( crbaHasEntry( sys, bodyOf[ r ], bodyOf[ c ] ) )
                row += format( "H[ {} ] = H_{}_{}; ", r * NV + c, r, c );
        if ( !row.empty() ) prln( "%s", row.c_str() );
    }

    if ( sys->loops.empty() ) {
        prln( "double rhs[ %d ] = {", NV );
        prln_indent_level++;
        for ( int k = 0 ; k < NV ; k++ ) {
            string rhs = format( "tau[ {} ] - C{}", k, k );
            if ( !sys->Qentry.empty() && sys->Qentry[ k ] != "0" ) rhs += " + " + sys->Qentry[ k ];
            prln( "%s,", rhs.c_str() );
        }
        prln_indent_level--;
        prln( "};" );

        prln( "ltdlFactor( %d, Hlam, H );", NV );
        prln( "ltdlSolveM( %d, Hlam, H, 1, rhs );", NV );
        for ( int k = 0 ; k < NV ; k++ ) prln( "out_qdd[ %d ] = rhs[ %d ];", k, k );

        prln_indent_level--;
        prln( "}" );
        return;
    }

    // Closed loops: method 2 of RBDA sec 8.5. One factorization of H serves
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
        string row = format( "tau[ {} ] - C{}", r, r );
        if ( !sys->Qentry.empty() && sys->Qentry[ r ] != "0" ) row += " + " + sys->Qentry[ r ];
        row += ", ";
        for ( int c = 0 ; c < NC ; c++ ) {
            const string &e = sys->Kentry[ c ][ r ];
            row += ( e == "0" ? "0.0"s : e ) + ", ";
        }
        prln( "%s", row.c_str() );
    }
    prln_indent_level--;
    prln( "};" );
    prln( "ltdlFactor( %d, Hlam, H );", NV );
    prln( "ltdlSolveM( %d, Hlam, H, %d, B );", NV, M );

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
