#pragma once

#include <stdio.h>
#include "prln.h"

// Emits the file-scope preamble: includes, and the linear solves used by
// gen_solve_qdd(). Must run before any other generator, since it's file-scope
// (not inside the generated function).
//
// H*qdd = tau - C goes through Featherstone's sparse L^T L factorization
// (RBDA sec 6.5). withLoops additionally emits the small dense solver, needed
// only by the method-2 loop path for its nc x nc constraint system.
void gen_prelude( FILE *f, int withLoops ) {
    prln_file = f;

    prln( "#include <math.h>" );
    prln( "#include <string.h>" );
    prln( "" );

    // ---- Sparse L^T D L factorization of H (RBDA sec 6.5, Tables 6.3 / 6.5).
    //
    // H is symmetric positive definite, so it needs no pivoting; and given the
    // expanded parent array lam (Table 6.4) the factorization has ZERO
    // FILL-IN -- L carries exactly H's branch-induced sparsity pattern, the
    // same pattern crbaHasEntry() tests. So the factorization never touches a
    // structural zero, and its cost is sum_i depth(i)^2 / 2 rather than the
    // n^3/3 of a dense elimination.
    //
    // Table 6.3 gives both L^T L and L^T D L; this is the latter, whose L is
    // unit lower-triangular so the factorization needs no square roots. That
    // is not a micro-optimization: with n sqrt() calls on the critical path,
    // L^T L measured *slower* than the dense elimination it replaced on every
    // model small enough to have no branch sparsity to win back (elbow,
    // stanford, rrrs and friends, all nv <= 9).
    //
    // The book indexes from 1 with lam(i) = 0 for "no parent"; this is the
    // 0-based form, with -1 for no parent. Only the lower triangle is read or
    // written, and lam[i] < i always holds (model.h requires parent < i).
    prln( "// Sparse L^T D L factorization of H, in place (RBDA Table 6.3): L is" );
    prln( "// unit lower-triangular, stored below the diagonal, and D on it." );
    prln( "// lam is the expanded parent array over velocity coordinates (RBDA" );
    prln( "// Table 6.4), which is exactly H's sparsity pattern: on row k, the" );
    prln( "// nonzeros below the diagonal are columns lam[k], lam[lam[k]], ..." );
    prln( "// H is symmetric positive definite, so no pivoting is needed, and" );
    prln( "// the factorization has zero fill-in -- the branch-induced zeros are" );
    prln( "// never touched. Only the lower triangle is accessed." );
    prln( "//" );
    prln( "// A singular H (a massless body, say) divides by zero and puts an inf" );
    prln( "// or NaN in out_qdd rather than reporting it: this runs in the inner" );
    prln( "// loop of a real-time integrator and must not do I/O. Validation" );
    prln( "// catches it immediately -- a NaN fails every tolerance." );
    prln( "static void ltdlFactor( int n, const int *lam, double *H ) {" );
    prln_indent_level++;
    prln( "for ( int k = n - 1 ; k >= 0 ; k-- ) {" );
    prln_indent_level++;
    prln( "for ( int i = lam[ k ] ; i != -1 ; i = lam[ i ] ) {" );
    prln_indent_level++;
    prln( "double a = H[ k*n + i ] / H[ k*n + k ];" );
    prln( "for ( int j = i ; j != -1 ; j = lam[ j ] ) H[ i*n + j ] -= a * H[ k*n + j ];" );
    prln( "H[ k*n + i ] = a;" );
    prln_indent_level--;
    prln( "}" );
    prln_indent_level--;
    prln( "}" );
    prln_indent_level--;
    prln( "}" );
    prln( "" );

    prln( "// Solves H X = B in place for X (left in B), given H already factored" );
    prln( "// by ltdlFactor. H = L^T D L, so X = L^-1 D^-1 ( L^-T B ): the two" );
    prln( "// in-situ back-substitutions of RBDA Table 6.5 with L_ii replaced by 1" );
    prln( "// (the unit-triangular form the table calls out), each walking only" );
    prln( "// the nonzero pattern via lam, with the diagonal solve folded into the" );
    prln( "// forward sweep. B is n x m row-major, and all m right-hand sides" );
    prln( "// share the one factorization -- which is what makes the method-2 loop" );
    prln( "// solve (RBDA sec 8.5) affordable for its 1 + nc columns." );
    prln( "static void ltdlSolveM( int n, const int *lam, const double *L, int m, double *B ) {" );
    prln_indent_level++;
    prln( "for ( int i = n - 1 ; i >= 0 ; i-- ) {            // B = L^-T B" );
    prln_indent_level++;
    prln( "for ( int j = lam[ i ] ; j != -1 ; j = lam[ j ] ) {" );
    prln_indent_level++;
    prln( "double lij = L[ i*n + j ];" );
    prln( "for ( int c = 0 ; c < m ; c++ ) B[ j*m + c ] -= lij * B[ i*m + c ];" );
    prln_indent_level--;
    prln( "}" );
    prln_indent_level--;
    prln( "}" );
    prln( "for ( int i = 0 ; i < n ; i++ ) {                 // B = L^-1 D^-1 B" );
    prln_indent_level++;
    prln( "double d = L[ i*n + i ];" );
    prln( "for ( int c = 0 ; c < m ; c++ ) B[ i*m + c ] /= d;" );
    prln( "for ( int j = lam[ i ] ; j != -1 ; j = lam[ j ] ) {" );
    prln_indent_level++;
    prln( "double lij = L[ i*n + j ];" );
    prln( "for ( int c = 0 ; c < m ; c++ ) B[ i*m + c ] -= lij * B[ j*m + c ];" );
    prln_indent_level--;
    prln( "}" );
    prln_indent_level--;
    prln( "}" );
    prln_indent_level--;
    prln( "}" );
    prln( "" );

    if ( !withLoops ) return;

    prln( "// Solves A x = b in place for x. A is n x n, row-major, and is" );
    prln( "// destroyed; b is destroyed too. Gaussian elimination with partial" );
    prln( "// pivoting. Used only for the method-2 loop path's nc x nc constraint" );
    prln( "// system A = K H^-1 K^T, which has no tree structure to exploit and is" );
    prln( "// small (nc <= 33 across every model here); H itself goes through the" );
    prln( "// sparse factorization above." );
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
