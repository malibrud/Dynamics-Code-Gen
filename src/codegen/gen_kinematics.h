#pragma once

#include <stdio.h>
#include "system.h"
#include "prln.h"
#include "svec.h"
#include "xform.h"
#include "inertia.h"

// Per-body kinematics and velocity emission, factored out of
// gen_fwd_dyn_crba() and gen_rnea_C() so a second generated function
// (gen_observe(), stage 5) can build the same S/vJ/xFromParent/xFromBase/v
// quantities from the same q/qd without duplicating the dynamics (CRBA H,
// RNEA bias force) it doesn't need. Both emit into whatever function body is
// currently open via prln_file; neither opens or closes a function.

// Body i's joint transform, motion subspace S, joint velocity vJ, and the
// composed xFromParent/xFromBase transforms -- exactly the per-body loop
// gen_fwd_dyn_crba() used to run inline, unchanged, so callers see identical
// generated code.
static void gen_kinematics_body( FILE *f, System *sys, int i ) {
    prln_file = f;
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

        // The body's ordinary right-handed attitude matrix R(Q), mapping
        // body coordinates to base coordinates -- the standard quaternion
        // rotation matrix, so Q is the usual body attitude quaternion.
        prln( "double Rfb%d_00 = 2.0*( Qw%d*Qw%d + Qx%d*Qx%d - 0.5 );", i, i, i, i, i );
        prln( "double Rfb%d_01 = 2.0*( Qx%d*Qy%d - Qw%d*Qz%d );", i, i, i, i, i );
        prln( "double Rfb%d_02 = 2.0*( Qx%d*Qz%d + Qw%d*Qy%d );", i, i, i, i, i );
        prln( "double Rfb%d_10 = 2.0*( Qx%d*Qy%d + Qw%d*Qz%d );", i, i, i, i, i );
        prln( "double Rfb%d_11 = 2.0*( Qw%d*Qw%d + Qy%d*Qy%d - 0.5 );", i, i, i, i, i );
        prln( "double Rfb%d_12 = 2.0*( Qy%d*Qz%d - Qw%d*Qx%d );", i, i, i, i, i );
        prln( "double Rfb%d_20 = 2.0*( Qx%d*Qz%d - Qw%d*Qy%d );", i, i, i, i, i );
        prln( "double Rfb%d_21 = 2.0*( Qy%d*Qz%d + Qw%d*Qx%d );", i, i, i, i, i );
        prln( "double Rfb%d_22 = 2.0*( Qw%d*Qw%d + Qz%d*Qz%d - 0.5 );", i, i, i, i, i );
        prln( "double rfb%d_0 = q[ %d ];", i, qi + 4 );
        prln( "double rfb%d_1 = q[ %d ];", i, qi + 5 );
        prln( "double rfb%d_2 = q[ %d ];", i, qi + 6 );

        Xform XJ;
        XJ.baseName = format( "Xjfb{}", i );
        // E is the coordinate transform (base-from-body), i.e. R^T -- a
        // pure index swap on the temps above, free at runtime. Same
        // reasoning as Xform::setJointRot() for the single-DOF joints.
        for ( int r = 0 ; r < 3 ; r++ )
        for ( int c = 0 ; c < 3 ; c++ )
            XJ.E[r][c] = format( "Rfb{}_{}{}", i, c, r );
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

// Runs gen_kinematics_body() for every body, in order.
static void gen_kinematics( FILE *f, System *sys ) {
    prln_file = f;
    int NB = sys->nBod;
    for ( int i = 0 ; i < NB ; i++ ) gen_kinematics_body( f, sys, i );
}

// Body i's spatial velocity v_i = i_X_lambda(i) * v_lambda(i) + vJ_i
// (v_lambda(i) = 0 at the base), from the kinematics gen_kinematics_body()
// already emitted (sys->xFromParent[i], sys->vJ[i]). Populates sys->v[i].
// Exactly the velocity half of the RNEA outward pass gen_rnea_C() used to
// run inline; the acceleration/force half stays in gen_rnea_C() since only
// the dynamics pass needs it.
static void gen_velocity_body( FILE *f, System *sys, int i ) {
    prln_file = f;
    int iParent = sys->parent[ i ];
    auto X6 = xformTo6x6( sys->xFromParent[ i ] );

    if ( iParent == -1 ) {
        sys->v[ i ] = sys->vJ[ i ];
    } else {
        Svec vParentInI = smat6MulVec( X6, sys->v[ iParent ] );
        sys->v[ i ] = svecAdd( "", vParentInI, sys->vJ[ i ] );
    }
    sys->v[ i ].baseName = format( "v{}_", i );
    sys->v[ i ].emit();
}
