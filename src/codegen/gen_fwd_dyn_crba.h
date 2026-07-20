#pragma once

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include "system.h"
#include "prln.h"
#include "gen_kinematics.h"

// #define VAR_SIZE 64
// typedef struct {
//     char s[ VAR_SIZE ];
// } Var;
// 
// Var _Var( const char *f, ... ) {
//     Var v;
//     va_list args;
//     va_start( args, f );
//     vsnprintf( v.s, VAR_SIZE, f, args );
//     va_end( args );
//     return v;
// }
// 
// typedef struct {
//     Var v[ 6 ];
// } Svec_;
// 
// typedef struct {
//     Var E[ 3 ][ 3 ];
//     Var r[ 3 ];
// } xForm;
// 
// Svec_ genS( int iBody, int jType ) {
//     Svec_ S;
//     switch( jType ) {
//         case J_REVOLUTE:
//             S.v[ 0 ] = _Var( "0" );
//             S.v[ 1 ] = _Var( "0" );
//             S.v[ 2 ] = _Var( "1" );
//             S.v[ 3 ] = _Var( "0" );
//             S.v[ 4 ] = _Var( "0" );
//             S.v[ 5 ] = _Var( "0" );
//             break;
//         default:
//             assert( false );
//     }
//     return S;
// }
// 
// Svec_ genVJ( int iBody, Svec_ *S, Var *qd ) {
//     Svec_ vJ;
//     for ( int i = 0 ; i < 6 ; i++ ) {
//         bool isElement0 = !strcmp( S->v[ i ].s, "0" );
//         vJ.v[ i ] = isElement0 ? _Var( "0") : *qd;
//     }
//     return vJ;
// }

void gen_fwd_dyn_crba( FILE *f, System *sys ) {
    prln_file = f;

    int NB = sys->nBod;
    prln( "// Coordinate counts baked into this generated function, so harnesses can" );
    prln( "// detect a dynamics.c that is stale relative to the mechanism file. nq and" );
    prln( "// nv differ when the mechanism has floating bodies: each takes 7 q slots" );
    prln( "// (unit quaternion w,x,y,z -- normalized on input -- then position x,y,z)" );
    prln( "// but 6 qd/tau/out_qdd slots (body-frame spatial velocity: wx,wy,wz,vx,vy,vz)." );
    prln( "const int fwd_dyn_crba_nbod = %d;", NB );
    prln( "const int fwd_dyn_crba_nq = %d;", sys->nq );
    prln( "const int fwd_dyn_crba_nv = %d;", sys->nv );
    prln( "const int fwd_dyn_crba_nact = %d;", sys->nact );
    prln( "" );
    prln( "// fext: external spatial forces, one per body, expressed in absolute (base)" );
    prln( "// coordinates: fext[ i*6 + k ] with k = 0..2 the moment about the base origin" );
    prln( "// and k = 3..5 the linear force (RBDA convention). Slots for bodies that were" );
    prln( "// generated without external-force support are never read, so fext may be" );
    prln( "// NULL when no body has it enabled." );
    prln( "// act: force-element inputs, fwd_dyn_crba_nact slots in element" );
    prln( "// declaration order (src/model/model.h's 'input' keyword); elements" );
    prln( "// without 'input' never read a slot, so act may be NULL when" );
    prln( "// fwd_dyn_crba_nact is 0 (same contract as fext)." );
    prln( "void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *act, double *out_qdd ) {" );
    prln_indent_level++;

    prln( "// The outward pass of the computation of C.  This is the first for loop of Table 5.1 in RBDA. ");
    prln( "" );
    gen_kinematics( f, sys );

    // Note: does not close the C function -- gen_crba_H() (see gen_crba.h)
    // appends the H-matrix computation into this same function body using
    // the xFromParent[]/S[] just populated above, and closes the brace.
}