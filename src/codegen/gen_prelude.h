#pragma once

#include <stdio.h>
#include "prln.h"

// Emits the file-scope preamble: includes, and a small generic dense linear
// solve used by gen_solve_qdd() to solve H*qdd = tau - C. Must run before any
// other generator, since it's file-scope (not inside the generated function).
// withLoops additionally emits a multiple-right-hand-side variant: the
// method-2 loop solve (RBDA sec 8.5) needs H^-1 applied to tau - C and to
// every column of K^T, and one elimination pass serves all of them.
void gen_prelude( FILE *f, int withLoops ) {
    prln_file = f;

    prln( "#include <math.h>" );
    prln( "#include <string.h>" );
    prln( "" );

    if ( withLoops ) {
        prln( "// Solves A X = B in place for X (left in B). A is n x n row-major and is" );
        prln( "// destroyed; B is n x m row-major. Gaussian elimination with partial" );
        prln( "// pivoting, one elimination shared by all m right-hand sides." );
        prln( "static void solveLinearSystemM( int n, double *A, int m, double *B ) {" );
        prln_indent_level++;
        prln( "for ( int col = 0 ; col < n ; col++ ) {" );
        prln_indent_level++;
        prln( "int pivot = col;" );
        prln( "double best = fabs( A[ col*n + col ] );" );
        prln( "for ( int row = col + 1 ; row < n ; row++ ) {" );
        prln_indent_level++;
        prln( "double v = fabs( A[ row*n + col ] );" );
        prln( "if ( v > best ) { best = v; pivot = row; }" );
        prln_indent_level--;
        prln( "}" );
        prln( "if ( pivot != col ) {" );
        prln_indent_level++;
        prln( "for ( int k = 0 ; k < n ; k++ ) {" );
        prln_indent_level++;
        prln( "double tmp = A[ col*n + k ]; A[ col*n + k ] = A[ pivot*n + k ]; A[ pivot*n + k ] = tmp;" );
        prln_indent_level--;
        prln( "}" );
        prln( "for ( int k = 0 ; k < m ; k++ ) {" );
        prln_indent_level++;
        prln( "double tmp = B[ col*m + k ]; B[ col*m + k ] = B[ pivot*m + k ]; B[ pivot*m + k ] = tmp;" );
        prln_indent_level--;
        prln( "}" );
        prln_indent_level--;
        prln( "}" );
        prln( "for ( int row = col + 1 ; row < n ; row++ ) {" );
        prln_indent_level++;
        prln( "double factor = A[ row*n + col ] / A[ col*n + col ];" );
        prln( "for ( int k = col ; k < n ; k++ ) A[ row*n + k ] -= factor * A[ col*n + k ];" );
        prln( "for ( int k = 0 ; k < m ; k++ ) B[ row*m + k ] -= factor * B[ col*m + k ];" );
        prln_indent_level--;
        prln( "}" );
        prln_indent_level--;
        prln( "}" );
        prln( "for ( int row = n - 1 ; row >= 0 ; row-- ) {" );
        prln_indent_level++;
        prln( "for ( int k = 0 ; k < m ; k++ ) {" );
        prln_indent_level++;
        prln( "double sum = B[ row*m + k ];" );
        prln( "for ( int c = row + 1 ; c < n ; c++ ) sum -= A[ row*n + c ] * B[ c*m + k ];" );
        prln( "B[ row*m + k ] = sum / A[ row*n + row ];" );
        prln_indent_level--;
        prln( "}" );
        prln_indent_level--;
        prln( "}" );
        prln_indent_level--;
        prln( "}" );
        prln( "" );
    }
    prln( "// Solves A x = b in place for x. A is n x n, row-major, and is" );
    prln( "// destroyed; b is destroyed too. Gaussian elimination with partial" );
    prln( "// pivoting -- fine for the small, dense systems produced here." );
    prln( "// Once bodies actually branch, replace this with the sparse LTDL" );
    prln( "// factorization that exploits branch-induced sparsity (see the" );
    prln( "// generated H's zero pattern, RBDA Table 6.5 / CRBA.pdf slide 13)." );
    prln( "static void solveLinearSystem( int n, double *A, double *b, double *x ) {" );
    prln_indent_level++;
    prln( "for ( int col = 0 ; col < n ; col++ ) {" );
    prln_indent_level++;
    prln( "int pivot = col;" );
    prln( "double best = fabs( A[ col*n + col ] );" );
    prln( "for ( int row = col + 1 ; row < n ; row++ ) {" );
    prln_indent_level++;
    prln( "double v = fabs( A[ row*n + col ] );" );
    prln( "if ( v > best ) { best = v; pivot = row; }" );
    prln_indent_level--;
    prln( "}" );
    prln( "if ( pivot != col ) {" );
    prln_indent_level++;
    prln( "for ( int k = 0 ; k < n ; k++ ) {" );
    prln_indent_level++;
    prln( "double tmp = A[ col*n + k ]; A[ col*n + k ] = A[ pivot*n + k ]; A[ pivot*n + k ] = tmp;" );
    prln_indent_level--;
    prln( "}" );
    prln( "double tmp = b[ col ]; b[ col ] = b[ pivot ]; b[ pivot ] = tmp;" );
    prln_indent_level--;
    prln( "}" );
    prln( "for ( int row = col + 1 ; row < n ; row++ ) {" );
    prln_indent_level++;
    prln( "double factor = A[ row*n + col ] / A[ col*n + col ];" );
    prln( "for ( int k = col ; k < n ; k++ ) A[ row*n + k ] -= factor * A[ col*n + k ];" );
    prln( "b[ row ] -= factor * b[ col ];" );
    prln_indent_level--;
    prln( "}" );
    prln_indent_level--;
    prln( "}" );
    prln( "for ( int row = n - 1 ; row >= 0 ; row-- ) {" );
    prln_indent_level++;
    prln( "double sum = b[ row ];" );
    prln( "for ( int k = row + 1 ; k < n ; k++ ) sum -= A[ row*n + k ] * x[ k ];" );
    prln( "x[ row ] = sum / A[ row*n + row ];" );
    prln_indent_level--;
    prln( "}" );
    prln_indent_level--;
    prln( "}" );
    prln( "" );
}
