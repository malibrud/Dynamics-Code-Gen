#pragma once

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include "system.h"
#include "prln.h"

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
    prln( "" );
    prln( "// fext: external spatial forces, one per body, expressed in absolute (base)" );
    prln( "// coordinates: fext[ i*6 + k ] with k = 0..2 the moment about the base origin" );
    prln( "// and k = 3..5 the linear force (RBDA convention). Slots for bodies that were" );
    prln( "// generated without external-force support are never read, so fext may be" );
    prln( "// NULL when no body has it enabled." );
    prln( "void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *out_qdd ) {" );
    prln_indent_level++;

    prln( "// The outward pass of the computation of C.  This is the first for loop of Table 5.1 in RBDA. ");
    prln( "" );
    for ( int i = 0 ; i < NB ; i++ ) {
        prln( "// Body %d", i );
        int jType = sys->jtype[ i ];
        char ax = sys->jaxis[ i ];
        int qi = sys->qidx[ i ];
        int vi = sys->vidx[ i ];

        if ( jType == J_FLOATING ) {
            // Genuine 6-DOF joint (RBDA sec 9.4): orientation is a quaternion
            // -- normalized here, so any nonzero quaternion is accepted and
            // there is no gimbal lock anywhere -- position is a plain
            // 3-vector, S is the 6x6 identity (left implicit; gen_rnea and
            // gen_crba special-case J_FLOATING), and vJ is the body-frame
            // spatial velocity read straight from qd.
            prln( "double quw%d = q[ %d ];", i, qi );
            prln( "double qux%d = q[ %d ];", i, qi + 1 );
            prln( "double quy%d = q[ %d ];", i, qi + 2 );
            prln( "double quz%d = q[ %d ];", i, qi + 3 );
            prln( "double qnorm%d = 1.0 / sqrt( quw%d*quw%d + qux%d*qux%d + quy%d*quy%d + quz%d*quz%d );",
                  i, i, i, i, i, i, i, i, i );
            prln( "double Qw%d = quw%d * qnorm%d;", i, i, i );
            prln( "double Qx%d = qux%d * qnorm%d;", i, i, i );
            prln( "double Qy%d = quy%d * qnorm%d;", i, i, i );
            prln( "double Qz%d = quz%d * qnorm%d;", i, i, i );

            // Rotation matrix in this codebase's E convention (the transpose
            // of spatial_v2's rq.m, matching how XFORM_RZ relates to rz.m).
            prln( "double Efb%d_00 = 2.0*( Qw%d*Qw%d + Qx%d*Qx%d - 0.5 );", i, i, i, i, i );
            prln( "double Efb%d_01 = 2.0*( Qx%d*Qy%d - Qw%d*Qz%d );", i, i, i, i, i );
            prln( "double Efb%d_02 = 2.0*( Qx%d*Qz%d + Qw%d*Qy%d );", i, i, i, i, i );
            prln( "double Efb%d_10 = 2.0*( Qx%d*Qy%d + Qw%d*Qz%d );", i, i, i, i, i );
            prln( "double Efb%d_11 = 2.0*( Qw%d*Qw%d + Qy%d*Qy%d - 0.5 );", i, i, i, i, i );
            prln( "double Efb%d_12 = 2.0*( Qy%d*Qz%d - Qw%d*Qx%d );", i, i, i, i, i );
            prln( "double Efb%d_20 = 2.0*( Qx%d*Qz%d - Qw%d*Qy%d );", i, i, i, i, i );
            prln( "double Efb%d_21 = 2.0*( Qy%d*Qz%d + Qw%d*Qx%d );", i, i, i, i, i );
            prln( "double Efb%d_22 = 2.0*( Qw%d*Qw%d + Qz%d*Qz%d - 0.5 );", i, i, i, i, i );
            prln( "double rfb%d_0 = q[ %d ];", i, qi + 4 );
            prln( "double rfb%d_1 = q[ %d ];", i, qi + 5 );
            prln( "double rfb%d_2 = q[ %d ];", i, qi + 6 );

            Xform XJ;
            XJ.baseName = format( "Xjfb{}", i );
            for ( int r = 0 ; r < 3 ; r++ )
            for ( int c = 0 ; c < 3 ; c++ )
                XJ.E[r][c] = format( "Efb{}_{}{}", i, r, c );
            for ( int k = 0 ; k < 3 ; k++ )
                XJ.r[k] = format( "rfb{}_{}", i, k );
            sys->xJoint[ i ] = XJ;

            Svec vJ( format( "vJ{}_", i ) );
            for ( int k = 0 ; k < 6 ; k++ ) {
                prln( "double vJ%d_%d = qd[ %d ];", i, k, vi + k );
                vJ.v[ k ] = format( "vJ{}_{}", i, k );
            }
            sys->vJ[ i ] = vJ;

            sys->xFromParent[ i ] = compose( format( "Xlink{}", i ), sys->xJoint[ i ], sys->xTree[ i ] );
            prln( "" );
            sys->xFromParent[ i ].emit();
            prln( "" );
        } else {
            string q = format( "q{}", i );
            string qd = format( "qd{}", i );

            sys->xFromParent[ i ] = compose( format( "Xlink{}", i ), sys->xJoint[ i ], sys->xTree[ i ] );

            // Motion subspace: revolute about axis a -> S = (a; 0); prismatic
            // along a -> S = (0; a).
            int sIdx = ( jType == J_REVOLUTE ? 0 : 3 ) + ( ax - 'x' );
            sys->S[ i ][ sIdx ] = "1";
            sys->S[ i ].baseName = format( "S{}_", i );

            sys->vJ[ i ] = sys->S[ i ] * qd;
            sys->vJ[ i ].baseName = format( "vJ{}_", i );

            prln( "double %s = q[ %d ];", q.c_str(), qi );
            prln( "double %s = qd[ %d ];", qd.c_str(), vi );
            if ( jType == J_REVOLUTE ) {
                prln( "double sin_%s = sin( %s );", q.c_str(), q.c_str() );
                prln( "double cos_%s = cos( %s );", q.c_str(), q.c_str() );
            }

            prln( "" );
            sys->S[ i ].emit();
            prln( "" );
            sys->vJ[ i ].emit();
            prln( "" );
            sys->xFromParent[ i ].emit();
            prln( "" );
        }

        // xFromBase is composed only after xFromParent has been emitted, so it
        // references the Xlink temps rather than re-expanding their expressions.
        // Compose with the actual parent's base transform -- not body i-1's,
        // which is the same thing only for unbranched chains. For a body whose
        // parent is the base, i_X_base is just Xlink_i -- already emitted, no
        // new temps needed.
        int iParent = sys->parent[ i ];
        if ( iParent != -1 ) {
            sys->xFromBase[ i ] = compose( format( "Xglobal{}", i ), sys->xFromParent[ i ], sys->xFromBase[ iParent ] );
            sys->xFromBase[ i ].emit();
            prln( "" );
        } else {
            sys->xFromBase[ i ] = sys->xFromParent[ i ];
        }

    }

    // Note: does not close the C function -- gen_crba_H() (see gen_crba.h)
    // appends the H-matrix computation into this same function body using
    // the xFromParent[]/S[] just populated above, and closes the brace.
}