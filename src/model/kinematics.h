// Shared spatial-algebra and forward-kinematics layer: everything needed to
// turn a (Model, q, qd) into per-body poses and velocities, with none of the
// dynamics.
//
// This was extracted verbatim from src/validate/dense.h, which is still its
// biggest consumer. It moved here because it has a second consumer that has
// nothing to do with validation: the track contact engine
// (src/contact/track_contact.h) needs to know where every shoe and roller is
// and how fast it is moving, and a shipped force producer must not depend on
// the validation reference.
//
// "dense" in denseKinematics() still means what it always did -- explicit 6x6
// Plucker matrices composed by plain matrix multiplication, as opposed to the
// generator's simplified straight-line code -- so the name is kept, and with
// it the property that anything cross-checked against it shares no code with
// src/codegen.
//
// Conventions (the ones that bite; see the note at the top of model.h):
//   - X = [ E, 0; -E*r~, E ] with E the COORDINATE transform (child-from-
//     parent), which is the TRANSPOSE of the child's physical rotation R.
//     Callers build the ordinary right-handed R and transpose it explicitly.
//     That is what makes q the true physical joint angle and qd = dq/dt.
//   - Spatial velocities are [ w ; vO ], angular first, vO the velocity of the
//     body-fixed point instantaneously at the frame origin.

#pragma once
#include <math.h>
#include "model.h"

// ---------------------------------------------------------------------------
// Generic dense linear algebra.
// ---------------------------------------------------------------------------
static void m6_mul( const double A[6][6], const double B[6][6], double C[6][6] ) {
    for ( int r = 0 ; r < 6 ; r++ )
    for ( int c = 0 ; c < 6 ; c++ ) {
        double s = 0;
        for ( int k = 0 ; k < 6 ; k++ ) s += A[r][k] * B[k][c];
        C[r][c] = s;
    }
}

static void m6_mulvec( const double A[6][6], const double x[6], double y[6] ) {
    for ( int r = 0 ; r < 6 ; r++ ) {
        double s = 0;
        for ( int c = 0 ; c < 6 ; c++ ) s += A[r][c] * x[c];
        y[r] = s;
    }
}

// ---------------------------------------------------------------------------
// Spatial-vector building blocks, straight from RBDA.
// ---------------------------------------------------------------------------
static void skew( const double v[3], double S[3][3] ) {
    S[0][0] = 0;     S[0][1] = -v[2]; S[0][2] =  v[1];
    S[1][0] =  v[2]; S[1][1] = 0;     S[1][2] = -v[0];
    S[2][0] = -v[1]; S[2][1] =  v[0]; S[2][2] = 0;
}

static void transpose3( const double A[3][3], double At[3][3] ) {
    for ( int i = 0 ; i < 3 ; i++ )
    for ( int j = 0 ; j < 3 ; j++ )
        At[i][j] = A[j][i];
}

// Plucker motion transform from the COORDINATE transform E (child-from-parent)
// and translation r (RBDA eq 2.24 with our generator's parametrization):
// X = [ E, 0; -E*rx, E ]. Note E is the transpose of the child's physical
// rotation -- callers build the ordinary R and transpose it (see jointXS).
static void plux( const double E[3][3], const double r[3], double X[6][6] ) {
    double rx[3][3], Erx[3][3];
    skew( r, rx );
    for ( int i = 0 ; i < 3 ; i++ )
    for ( int j = 0 ; j < 3 ; j++ ) {
        double s = 0;
        for ( int k = 0 ; k < 3 ; k++ ) s += E[i][k] * rx[k][j];
        Erx[i][j] = s;
    }
    for ( int i = 0 ; i < 3 ; i++ )
    for ( int j = 0 ; j < 3 ; j++ ) {
        X[ i   ][ j   ] = E[i][j];
        X[ i   ][ j+3 ] = 0;
        X[ i+3 ][ j   ] = -Erx[i][j];
        X[ i+3 ][ j+3 ] = E[i][j];
    }
}

// Unit vector for a canonical joint axis.
static void axisVec( char axis, double a[3] ) {
    a[0] = a[1] = a[2] = 0;
    if      ( axis == 'x' ) a[0] = 1;
    else if ( axis == 'y' ) a[1] = 1;
    else                    a[2] = 1;
}

// The ordinary right-handed rotation matrix R(q) about a canonical axis --
// the physical rotation the joint applies to its child frame, matching
// XFORM_RX/RY/RZ in the generator. This is NOT the Plucker E block; callers
// transpose it (see jointXS).
static void rotAxis( char axis, double q, double R[3][3] ) {
    double c = cos( q ), s = sin( q );
    if ( axis == 'x' ) {
        R[0][0] = 1; R[0][1] = 0; R[0][2] =  0;
        R[1][0] = 0; R[1][1] = c; R[1][2] = -s;
        R[2][0] = 0; R[2][1] = s; R[2][2] =  c;
    } else if ( axis == 'y' ) {
        R[0][0] =  c; R[0][1] = 0; R[0][2] = s;
        R[1][0] =  0; R[1][1] = 1; R[1][2] = 0;
        R[2][0] = -s; R[2][1] = 0; R[2][2] = c;
    } else {
        R[0][0] = c; R[0][1] = -s; R[0][2] = 0;
        R[1][0] = s; R[1][1] =  c; R[1][2] = 0;
        R[2][0] = 0; R[2][1] =  0; R[2][2] = 1;
    }
}

// The body-to-base rotation R of a floating body's quaternion slice (w,x,y,z),
// normalized here. Split out of jointXS because callers outside the dynamics
// need the physical attitude directly: the integrator's rdot = R(Q)*v_body
// (src/sim/integrate.h) and any consumer converting a body-frame vector to
// world. jointXS transposes it into the Plucker E block; these callers do not.
static void quatToR( const double qseg[4], double R[3][3] ) {
    double w = qseg[0], x = qseg[1], y = qseg[2], z = qseg[3];
    double n = 1.0 / sqrt( w*w + x*x + y*y + z*z );
    w *= n; x *= n; y *= n; z *= n;
    R[0][0] = 2.0*( w*w + x*x - 0.5 ); R[0][1] = 2.0*( x*y - w*z );       R[0][2] = 2.0*( x*z + w*y );
    R[1][0] = 2.0*( x*y + w*z );       R[1][1] = 2.0*( w*w + y*y - 0.5 ); R[1][2] = 2.0*( y*z - w*x );
    R[2][0] = 2.0*( x*z - w*y );       R[2][1] = 2.0*( y*z + w*x );       R[2][2] = 2.0*( w*w + z*z - 0.5 );
}

// Joint transform XJ and motion subspace S (6 x ndof, stored as columns of a
// 6x6 array) for a joint at position coordinates qseg (this body's slice of
// the packed q vector). Revolute: rotation about the axis, S = (axis; 0).
// Prismatic: translation along the axis, S = (0; axis). Floating: quaternion
// orientation plus position, S = identity(6) -- a genuine 6-DOF joint with
// no gimbal lock.
//
// In every case the physical rotation is built as the ordinary right-handed R
// and transposed into the Plucker E block (child-from-parent coordinates), so
// q is the true physical joint coordinate and S stays on the positive axis --
// which is what makes qd = dq/dt. See Xform::setJointRot() in the generator.
static int jointXS( const ModelBody *b, const double *qseg, double XJ[6][6], double S[6][6] ) {
    double a3[3], Eident[3][3] = { {1,0,0},{0,1,0},{0,0,1} }, r0[3] = { 0, 0, 0 };
    for ( int r = 0 ; r < 6 ; r++ )
    for ( int c = 0 ; c < 6 ; c++ ) S[r][c] = 0;

    if ( b->jtype == JOINT_FLOATING ) {
        // The quaternion is the usual body attitude (body-to-base).
        double R[3][3], E[3][3];
        quatToR( qseg, R );
        transpose3( R, E );
        plux( E, qseg + 4, XJ );
        for ( int k = 0 ; k < 6 ; k++ ) S[k][k] = 1;
        return 6;
    }

    axisVec( b->axis, a3 );
    if ( b->jtype == JOINT_REVOLUTE ) {
        double R[3][3], E[3][3];
        rotAxis( b->axis, qseg[0], R );
        transpose3( R, E );
        plux( E, r0, XJ );
        for ( int k = 0 ; k < 3 ; k++ ) S[k][0] = a3[k];
    } else { // JOINT_PRISMATIC (spherical/cylindrical rejected upstream)
        double rq[3] = { a3[0]*qseg[0], a3[1]*qseg[0], a3[2]*qseg[0] };
        plux( Eident, rq, XJ );
        for ( int k = 0 ; k < 3 ; k++ ) S[k+3][0] = a3[k];
    }
    return 1;
}

// Spatial inertia from mass, CG offset, and rotational inertia about the CG
// (RBDA eq 2.63, spatial_v2's mcI):
//   I = [ Icm + m*cx*cx^T,  m*cx ;
//         m*cx^T,           m*1  ]     with cx = skew(c), cx^T = -cx.
static void mcI( double m, const double c[3], const double Ic3[6], double I[6][6] ) {
    double cx[3][3];
    skew( c, cx );

    double Ibar[3][3] = {
        { Ic3[0], Ic3[3], Ic3[5] },
        { Ic3[3], Ic3[1], Ic3[4] },
        { Ic3[5], Ic3[4], Ic3[2] },
    };
    // Ibar += m * cx * cx^T
    for ( int i = 0 ; i < 3 ; i++ )
    for ( int j = 0 ; j < 3 ; j++ ) {
        double s = 0;
        for ( int k = 0 ; k < 3 ; k++ ) s += cx[i][k] * cx[j][k];
        Ibar[i][j] += m * s;
    }

    for ( int i = 0 ; i < 3 ; i++ )
    for ( int j = 0 ; j < 3 ; j++ ) {
        I[ i   ][ j   ] = Ibar[i][j];
        I[ i   ][ j+3 ] =  m * cx[i][j];
        I[ i+3 ][ j   ] = -m * cx[i][j]; // m * cx^T
        I[ i+3 ][ j+3 ] = ( i == j ) ? m : 0;
    }
}

// Origin of the frame a base-to-frame transform X maps into, in base
// coordinates: X = [E 0; B E] with B = -E*r~, so r~ = -E^T B (the same
// recovery point-on-line's delta uses, here for all three components).
static void xformOrigin( const double X[6][6], double p[3] ) {
    double rx = 0, ry = 0, rz = 0;
    for ( int k = 0 ; k < 3 ; k++ ) {
        rx -= X[k][2] * X[k+3][1]; // (r~)[2][1]
        ry -= X[k][0] * X[k+3][2]; // (r~)[0][2]
        rz -= X[k][1] * X[k+3][0]; // (r~)[1][0]
    }
    p[0] = rx; p[1] = ry; p[2] = rz;
}

// ---------------------------------------------------------------------------
// Kinematics-only pass -- Xa[i] (i_X_0) and v[i] (body-frame spatial
// velocity) for every body -- mirroring the first loop of fd_dense_full but
// stopping there (no inertia, no bias force). Feeds denseObserveTsda in
// dense.h (the cross-check for the generator's fwd_dyn_observe) and the track
// contact engine.
//
// Both outputs are indexed by ABSOLUTE body index and must be sized
// MODEL_MAX_BODIES. Callers: make them static -- Xa alone is 590 KB at
// MODEL_MAX_BODIES = 256 and will blow a default stack.
// ---------------------------------------------------------------------------
static void denseKinematics( const Model *mdl, const double *q, const double *qd,
                             double Xa[][6][6], double v[][6] ) {
    int nb = mdl->nbod;
    for ( int i = 0 ; i < nb ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];
        int par = b->parent;
        int vi = mdl->vidx[ i ];

        double Eident[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
        double XJ[6][6], XT[6][6], Xup[6][6], S[6][6];
        int ndof = jointXS( b, q + mdl->qidx[ i ], XJ, S );
        plux( Eident, b->tree, XT );
        m6_mul( XJ, XT, Xup );

        if ( par == -1 ) {
            for ( int r = 0 ; r < 6 ; r++ ) for ( int c = 0 ; c < 6 ; c++ ) Xa[i][r][c] = Xup[r][c];
        } else {
            m6_mul( Xup, Xa[par], Xa[i] );
        }

        double vJ[6] = { 0, 0, 0, 0, 0, 0 };
        for ( int r = 0 ; r < 6 ; r++ )
            for ( int k = 0 ; k < ndof ; k++ ) vJ[r] += S[r][k] * qd[ vi + k ];
        if ( par == -1 ) {
            for ( int k = 0 ; k < 6 ; k++ ) v[i][k] = vJ[k];
        } else {
            m6_mulvec( Xup, v[par], v[i] );
            for ( int k = 0 ; k < 6 ; k++ ) v[i][k] += vJ[k];
        }
    }
}
