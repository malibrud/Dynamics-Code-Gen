// Dense-matrix reference implementation of the forward dynamics, shared by
// the validator (src/validate/main.c) and the benchmark (src/bench/main.c).
//
// This implements the same computation as the code generator
// (src/codegen -> build/dynamics.c) but with full 6x6 spatial matrices and
// 6-vectors, straight from the equations in RBDA:
//
//   - motion transforms built as explicit 6x6 Plucker matrices and composed
//     by plain matrix multiplication (RBDA sec 2.8),
//   - spatial inertias built with mcI (RBDA eq 2.63),
//   - cross products as the crm/crf 6x6 operator matrices (RBDA eq 2.31/2.32),
//   - force coordinate changes done by solving X^T x = f with a generic
//     Gaussian solver (the Xa{i}' \ f_ext{i} of spatial_v2).
//
// Nothing here shares code with the generator, so agreement between the two
// is a real check: the generator's symbolic simplification (structural zeros,
// materialized temps, hand-derived block forms) against brute-force algebra.
//
// The mechanism comes from a model file (src/model/model.h). This reference
// handles any single-DOF joint the format can express -- revolute or
// prismatic about x, y, or z, arbitrary parent trees -- though only what the
// generator also supports is cross-validated today.

#pragma once
#include <math.h>
#include "../model/model.h"

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

// C = A^T * B * A (the congruence transform used by the CRBA inertia fold).
static void m6_congruence( const double A[6][6], const double B[6][6], double C[6][6] ) {
    double BA[6][6];
    m6_mul( B, A, BA );
    for ( int r = 0 ; r < 6 ; r++ )
    for ( int c = 0 ; c < 6 ; c++ ) {
        double s = 0;
        for ( int k = 0 ; k < 6 ; k++ ) s += A[k][r] * BA[k][c];
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

// y = A^T x
static void m6_tmulvec( const double A[6][6], const double x[6], double y[6] ) {
    for ( int r = 0 ; r < 6 ; r++ ) {
        double s = 0;
        for ( int c = 0 ; c < 6 ; c++ ) s += A[c][r] * x[c];
        y[r] = s;
    }
}

static double dot6( const double a[6], const double b[6] ) {
    double s = 0;
    for ( int k = 0 ; k < 6 ; k++ ) s += a[k] * b[k];
    return s;
}

// Solves A x = b (n x n, row-major) by Gaussian elimination with partial
// pivoting. A and b are destroyed.
static void linsolve( int n, double *A, double *b, double *x ) {
    for ( int col = 0 ; col < n ; col++ ) {
        int pivot = col;
        double best = fabs( A[ col*n + col ] );
        for ( int row = col + 1 ; row < n ; row++ ) {
            double v = fabs( A[ row*n + col ] );
            if ( v > best ) { best = v; pivot = row; }
        }
        if ( pivot != col ) {
            for ( int k = 0 ; k < n ; k++ ) {
                double t = A[ col*n + k ]; A[ col*n + k ] = A[ pivot*n + k ]; A[ pivot*n + k ] = t;
            }
            double t = b[ col ]; b[ col ] = b[ pivot ]; b[ pivot ] = t;
        }
        for ( int row = col + 1 ; row < n ; row++ ) {
            double f = A[ row*n + col ] / A[ col*n + col ];
            for ( int k = col ; k < n ; k++ ) A[ row*n + k ] -= f * A[ col*n + k ];
            b[ row ] -= f * b[ col ];
        }
    }
    for ( int row = n - 1 ; row >= 0 ; row-- ) {
        double s = b[ row ];
        for ( int k = row + 1 ; k < n ; k++ ) s -= A[ row*n + k ] * x[ k ];
        x[ row ] = s / A[ row*n + row ];
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
        // Unit quaternion (w,x,y,z), normalized here. R is the standard
        // quaternion rotation matrix, so the quaternion is the usual body
        // attitude (body-to-base).
        double w = qseg[0], x = qseg[1], y = qseg[2], z = qseg[3];
        double n = 1.0 / sqrt( w*w + x*x + y*y + z*z );
        w *= n; x *= n; y *= n; z *= n;
        double R[3][3], E[3][3];
        R[0][0] = 2.0*( w*w + x*x - 0.5 ); R[0][1] = 2.0*( x*y - w*z );       R[0][2] = 2.0*( x*z + w*y );
        R[1][0] = 2.0*( x*y + w*z );       R[1][1] = 2.0*( w*w + y*y - 0.5 ); R[1][2] = 2.0*( y*z - w*x );
        R[2][0] = 2.0*( x*z - w*y );       R[2][1] = 2.0*( y*z + w*x );       R[2][2] = 2.0*( w*w + z*z - 0.5 );
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

// Motion cross-product operator crm(v) (RBDA eq 2.31): [ wx, 0; v0x, wx ].
static void crm( const double v[6], double M[6][6] ) {
    double wx[3][3], v0x[3][3];
    skew( v, wx );
    skew( v + 3, v0x );
    for ( int i = 0 ; i < 3 ; i++ )
    for ( int j = 0 ; j < 3 ; j++ ) {
        M[ i   ][ j   ] = wx[i][j];
        M[ i   ][ j+3 ] = 0;
        M[ i+3 ][ j   ] = v0x[i][j];
        M[ i+3 ][ j+3 ] = wx[i][j];
    }
}

// Force cross-product operator crf(v) = -crm(v)^T (RBDA eq 2.32).
static void crf( const double v[6], double M[6][6] ) {
    double cm[6][6];
    crm( v, cm );
    for ( int i = 0 ; i < 6 ; i++ )
    for ( int j = 0 ; j < 6 ; j++ )
        M[i][j] = -cm[j][i];
}

// ---------------------------------------------------------------------------
// Loop-joint jcalc (RBDA Table 8.1 / sec 8.3): given the measured
// predecessor->successor transform X = sXp of a loop joint (canonical axis =
// z; the caller folds any x/y axis into XP/XS as a fixed rotation), returns
// the constraint-force subspace T (6 x nc columns) and the position error
// delta = T^T d, where d is the small displacement from where the successor
// frame ought to be to where it actually is.
//
// For every Table 8.1 type, T is constant in SUCCESSOR coordinates and the
// caller transforms it to base with Xs^T. point-on-line (sphere-in-cylinder,
// not in Table 8.1) is the exception: its constraint directions follow the
// predecessor's axis, so T is constant in PREDECESSOR coordinates and the
// caller must use Xp^T instead (loopTInPred below). The k_l formula (RBDA eq
// 8.47) is unchanged either way: Tdot = v x* T with v of whichever body T is
// fixed to, and the cross-product antisymmetry makes both reductions agree.
//
// point-on-line's delta is exact, not small-angle: X = [E 0; -E*rho~ E] for
// some rho (successor origin in predecessor coords), so rho~ = -E^T B
// recovers the origin position exactly, and the error is its x/y part --
// insensitive to the joint's free rotations and free z translation.
// ---------------------------------------------------------------------------
static int loopJcalc( ModelJointType t, const double X[6][6], double T[6][6], double delta[6] ) {
    for ( int r = 0 ; r < 6 ; r++ )
    for ( int c = 0 ; c < 6 ; c++ ) T[r][c] = 0;

    switch ( t ) {
    case JOINT_REVOLUTE:   // free: rotation z.  T = [e1 e2 e4 e5 e6]
        T[0][0] = 1; T[1][1] = 1; T[3][2] = 1; T[4][3] = 1; T[5][4] = 1;
        delta[0] = X[1][2];
        delta[1] = -X[0][2];
        delta[2] = X[4][2];
        delta[3] = -X[3][2];
        delta[4] = X[3][0]*X[1][0] + X[3][1]*X[1][1];
        return 5;
    case JOINT_PRISMATIC:  // free: translation z.  T = [e1 e2 e3 e4 e5]
        T[0][0] = 1; T[1][1] = 1; T[2][2] = 1; T[3][3] = 1; T[4][4] = 1;
        delta[0] = 0.5 * ( X[1][2] - X[2][1] );
        delta[1] = 0.5 * ( X[2][0] - X[0][2] );
        delta[2] = 0.5 * ( X[0][1] - X[1][0] );
        delta[3] = X[4][2];
        delta[4] = -X[3][2];
        return 5;
    case JOINT_CYLINDRICAL: // free: rotation z + translation z.  T = [e1 e2 e4 e5]
        T[0][0] = 1; T[1][1] = 1; T[3][2] = 1; T[4][3] = 1;
        delta[0] = X[1][2];
        delta[1] = -X[0][2];
        delta[2] = X[4][2];
        delta[3] = -X[3][2];
        return 4;
    case JOINT_SPHERICAL:  // free: all rotations.  T = [e4 e5 e6]
        T[3][0] = 1; T[4][1] = 1; T[5][2] = 1;
        delta[0] = X[4][0]*X[2][0] + X[4][1]*X[2][1] + X[4][2]*X[2][2];
        delta[1] = X[5][0]*X[0][0] + X[5][1]*X[0][1] + X[5][2]*X[0][2];
        delta[2] = X[3][0]*X[1][0] + X[3][1]*X[1][1] + X[3][2]*X[1][2];
        return 3;
    default: {             // JOINT_POINTONLINE. free: rotations + translation z.
        // T = [e4 e5] in PREDECESSOR coordinates (see header comment).
        T[3][0] = 1; T[4][1] = 1;
        // X = [E 0; B E] with B = -E*rho~, so rho~ = -E^T B, where rho is the
        // successor origin in predecessor coords. With
        // skew(rho) = [[0,-rz,ry],[rz,0,-rx],[-ry,rx,0]] (0-based indices):
        //   rho_x = (rho~)[2][1] = -(E^T B)[2][1] = -sum_k E[k][2]*B[k][1]
        //   rho_y = (rho~)[0][2] = -(E^T B)[0][2] = -sum_k E[k][0]*B[k][2]
        double rx = 0, ry = 0;
        for ( int k = 0 ; k < 3 ; k++ ) {
            rx -= X[k][2] * X[k+3][1];
            ry -= X[k][0] * X[k+3][2];
        }
        delta[0] = rx;
        delta[1] = ry;
        return 2;
    }
    }
}

// True if this loop type's T is fixed to the predecessor (transform with
// Xp^T) rather than the successor (Xs^T).
static int loopTInPred( ModelJointType t ) { return t == JOINT_POINTONLINE; }

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

// Coordinate rotation E aligning the canonical loop-joint z axis with the
// declared axis: E = R_axis^T where R_axis maps joint coords to body coords
// by cyclic permutation (x: jz->bx, y: jz->by, z: identity). Cyclic keeps
// the frames right-handed.
static void axisAlignE( char axis, double E[3][3] ) {
    for ( int r = 0 ; r < 3 ; r++ )
    for ( int c = 0 ; c < 3 ; c++ ) E[r][c] = 0;
    if ( axis == 'x' ) {        // R cols: jx->by, jy->bz, jz->bx
        E[0][1] = 1; E[1][2] = 1; E[2][0] = 1;
    } else if ( axis == 'y' ) { // R cols: jx->bz, jy->bx, jz->by
        E[0][2] = 1; E[1][0] = 1; E[2][1] = 1;
    } else {
        E[0][0] = 1; E[1][1] = 1; E[2][2] = 1;
    }
}

// World position of a poi's local point, from the array of base->body
// transforms fd_dense_full already built (Xa[i] = i_X_0): p = r + E^T*t,
// with r the frame's origin (xformOrigin) and E the transform's upper-left
// 3x3 block -- mirrors the generator's loopAnchorBase (gen_loops.h), but
// deliberately kept as separate, independent code (see fd_dense_full's
// header comment on why elements are applied via a different route there).
static void poiWorldPos( const ModelPoi *poi, double Xa[][6][6], double p[3] ) {
    if ( poi->body == -1 ) { p[0] = poi->p[0]; p[1] = poi->p[1]; p[2] = poi->p[2]; return; }
    double r[3];
    xformOrigin( Xa[ poi->body ], r );
    for ( int k = 0 ; k < 3 ; k++ ) {
        double s = r[k];
        for ( int i = 0 ; i < 3 ; i++ ) s += Xa[ poi->body ][i][k] * poi->p[i];
        p[k] = s;
    }
}

// A poi's body's base-coordinate spatial velocity [w; vO] (vO = velocity of
// the point instantaneously coincident with the base origin), solved rather
// than inverted -- same recipe fd_dense_full already uses for loop anchors.
// Zero for a poi fixed to the base.
static void poiWorldVel( const ModelPoi *poi, double Xa[][6][6], double v[][6], double vout[6] ) {
    if ( poi->body == -1 ) { for ( int k = 0 ; k < 6 ; k++ ) vout[k] = 0; return; }
    double M[36], b6[6];
    for ( int r = 0 ; r < 6 ; r++ ) for ( int c = 0 ; c < 6 ; c++ ) M[ r*6 + c ] = Xa[ poi->body ][r][c];
    for ( int k = 0 ; k < 6 ; k++ ) b6[k] = v[ poi->body ][k];
    linsolve( 6, M, b6, vout );
}

// ---------------------------------------------------------------------------
// Dense forward dynamics: H*qdd = tau - C, with C from the RNEA (RBDA Table
// 5.1, aJ = 0) and H from the CRBA (RBDA Table 6.4). Joints may have 1 DOF
// (revolute/prismatic) or 6 (floating), so H is nv x nv with per-joint
// blocks; q is addressed through mdl->qidx and qd/tau/qdd through mdl->vidx.
// The large working arrays are static: fd_dense is not reentrant, but stays
// comfortably off the stack even at MODEL_MAX_V.
//
// With kinematic loops (mdl->nloop > 0) this becomes RBDA ch. 8 method 2:
// build K, k (Table 8.4, in base coordinates), solve A*lambda = b with
// A = K H^-1 K^T (eq 8.38), then H qdd = tau - C + K^T lambda. K_out, if
// non-NULL, receives the nc x nv constraint matrix (row-major) so callers
// can rank-check it -- the pipeline requires properly constrained loops,
// rank(K) = nc. delta_out, if non-NULL, receives the nc position errors
// (Table 8.1's delta, stacked per loop): K is d(delta)/dq to first order,
// so the pair supports Newton-projecting a state onto the closure manifold.
//
// Force elements (TSDA/RSDA, src/model/model.h) are applied via a route
// deliberately DIFFERENT from the generator's: the generator's Q is the
// TSDA tension projected through the distance loop joint's path-walk row
// (gen_elems.h); here each TSDA becomes an explicit equal-and-opposite
// spatial force at its two anchors in base coordinates, folded into the
// per-body bias force exactly like fext (but always on, not opt-in). RSDA
// is a direct torque into its own coordinate either way, since there is no
// second formulation of "torque on your own joint" -- act is NULL-able
// exactly like fext, and only read when some element actually has 'input'.
// ---------------------------------------------------------------------------
static void fd_dense_full( const Model *mdl, const double *q, const double *qd,
                           const double *tau, const double *fext, const double *act,
                           double *qdd_out, double *K_out, double *delta_out ) {
    int nb = mdl->nbod;
    int nv = mdl->nv;
    double Xup[ MODEL_MAX_BODIES ][6][6]; // i_X_lambda(i)
    double Xa[ MODEL_MAX_BODIES ][6][6];  // i_X_0
    double I6[ MODEL_MAX_BODIES ][6][6];
    double S[ MODEL_MAX_BODIES ][6][6];   // motion subspace, ndof[i] columns
    int ndof[ MODEL_MAX_BODIES ];
    double v[ MODEL_MAX_BODIES ][6], avp[ MODEL_MAX_BODIES ][6], fvp[ MODEL_MAX_BODIES ][6];

    // -a_grav, the seed of the RNEA acceleration recursion.
    double gravNeg[6] = { 0, 0, 0, -mdl->gravity[0], -mdl->gravity[1], -mdl->gravity[2] };

    for ( int i = 0 ; i < nb ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];
        int par = b->parent;
        int vi = mdl->vidx[ i ];

        mcI( b->mass, b->cg, b->inertia, I6[i] );

        // Xup = XJ * XT, exactly as the generator composes
        // Xlink = xJoint o xTree.
        double Eident[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
        double XJ[6][6], XT[6][6];
        ndof[i] = jointXS( b, q + mdl->qidx[ i ], XJ, S[i] );
        plux( Eident, b->tree, XT );
        m6_mul( XJ, XT, Xup[i] );

        if ( par == -1 ) {
            for ( int r = 0 ; r < 6 ; r++ )
            for ( int c = 0 ; c < 6 ; c++ ) Xa[i][r][c] = Xup[i][r][c];
        } else {
            m6_mul( Xup[i], Xa[par], Xa[i] );
        }

        // v_i = Xup * v_parent + vJ,  vJ = S * qd-slice
        double vJ[6] = { 0, 0, 0, 0, 0, 0 };
        for ( int r = 0 ; r < 6 ; r++ )
            for ( int k = 0 ; k < ndof[i] ; k++ ) vJ[r] += S[i][r][k] * qd[ vi + k ];
        if ( par == -1 ) {
            for ( int k = 0 ; k < 6 ; k++ ) v[i][k] = vJ[k];
        } else {
            m6_mulvec( Xup[i], v[par], v[i] );
            for ( int k = 0 ; k < 6 ; k++ ) v[i][k] += vJ[k];
        }

        // avp_i = Xup * avp_parent + crm(v_i)*vJ   (avp_parent(base) = -a_grav)
        // For a floating root v == vJ, so the cross term is zero -- crm(v)*v
        // cancels exactly in floating point, so no special case is needed.
        double crmV[6][6], vCrossVJ[6];
        m6_mulvec( Xup[i], par == -1 ? gravNeg : avp[par], avp[i] );
        crm( v[i], crmV );
        m6_mulvec( crmV, vJ, vCrossVJ );
        for ( int k = 0 ; k < 6 ; k++ ) avp[i][k] += vCrossVJ[k];

        // fvp_i = I*avp + crf(v)*(I*v) - Xa^T \ fext_i
        double Ia[6], Iv[6], crfV[6][6], vIv[6];
        m6_mulvec( I6[i], avp[i], Ia );
        m6_mulvec( I6[i], v[i], Iv );
        crf( v[i], crfV );
        m6_mulvec( crfV, Iv, vIv );
        for ( int k = 0 ; k < 6 ; k++ ) fvp[i][k] = Ia[k] + vIv[k];

        // Bodies generated without external-force support ignore their fext
        // slots entirely (spatial_v2's empty f_ext{i} cells) -- mirror that
        // here or the two implementations diverge whenever the harness fills
        // fext for a body whose flag is off.
        if ( b->hasFext ) {
            double XaT[36], bvec[6], fb[6];
            for ( int r = 0 ; r < 6 ; r++ )
            for ( int c = 0 ; c < 6 ; c++ ) XaT[ r*6 + c ] = Xa[i][c][r];
            for ( int k = 0 ; k < 6 ; k++ ) bvec[k] = fext[ i*6 + k ];
            linsolve( 6, XaT, bvec, fb );
            for ( int k = 0 ; k < 6 ; k++ ) fvp[i][k] -= fb[k];
        }
    }

    // Force elements: TSDA becomes equal-and-opposite spatial forces at its
    // two anchors, base coordinates -- deliberately not the generator's
    // path-walk row (see the header comment above). Runs after every body's
    // Xa/v is known, since an element's anchors can be on any two bodies,
    // not just parent-child (unlike fext, which only ever needs its own
    // body's Xa and so lives inside the loop above).
    static double elemForceBase[ MODEL_MAX_BODIES ][6];
    static double elemQ[ MODEL_MAX_V ];
    for ( int i = 0 ; i < nb ; i++ ) for ( int k = 0 ; k < 6 ; k++ ) elemForceBase[i][k] = 0;
    for ( int k = 0 ; k < nv ; k++ ) elemQ[k] = 0;

    for ( int e = 0 ; e < mdl->nelem ; e++ ) {
        const ModelElem *el = &mdl->elem[ e ];

        if ( el->type == ELEM_RSDA ) {
            // Direct torque into its own coordinate: tau_e = -k*(q - q0) -
            // c*qd + act[slot]. No second formulation exists for "a torque
            // on your own joint", so this matches the generator exactly.
            int b = el->joint;
            double qq  = q[ mdl->qidx[ b ] ];
            double qqd = qd[ mdl->vidx[ b ] ];
            double tau_e = -el->k * ( qq - el->ref ) - el->c * qqd;
            if ( el->hasInput ) tau_e += act[ el->actIdx ];
            elemQ[ mdl->vidx[ b ] ] += tau_e;
            continue;
        }

        // TSDA: tension f = k*(l - l0) + c*ldot + act[slot], positive f
        // pulls the anchors together. Force on "to" points back toward
        // "from"; equal and opposite on "from".
        const ModelPoi *pf = &mdl->poi[ el->poi1 ];
        const ModelPoi *pt = &mdl->poi[ el->poi2 ];
        double pp[3], ps[3], rho[3];
        poiWorldPos( pf, Xa, pp );
        poiWorldPos( pt, Xa, ps );
        for ( int k = 0 ; k < 3 ; k++ ) rho[k] = ps[k] - pp[k];
        double l = sqrt( rho[0]*rho[0] + rho[1]*rho[1] + rho[2]*rho[2] );
        double linv = 1.0 / l;

        double vp6[6], vs6[6];
        poiWorldVel( pf, Xa, v, vp6 );
        poiWorldVel( pt, Xa, v, vs6 );
        double ppd[3], psd[3], rhod[3];
        ppd[0] = vp6[3] + vp6[1]*pp[2] - vp6[2]*pp[1];
        ppd[1] = vp6[4] + vp6[2]*pp[0] - vp6[0]*pp[2];
        ppd[2] = vp6[5] + vp6[0]*pp[1] - vp6[1]*pp[0];
        psd[0] = vs6[3] + vs6[1]*ps[2] - vs6[2]*ps[1];
        psd[1] = vs6[4] + vs6[2]*ps[0] - vs6[0]*ps[2];
        psd[2] = vs6[5] + vs6[0]*ps[1] - vs6[1]*ps[0];
        for ( int k = 0 ; k < 3 ; k++ ) rhod[k] = psd[k] - ppd[k];
        double ldot = ( rho[0]*rhod[0] + rho[1]*rhod[1] + rho[2]*rhod[2] ) * linv;

        double ff = el->k * ( l - el->ref ) + el->c * ldot;
        if ( el->hasInput ) ff += act[ el->actIdx ];

        double Fto[3] = { -ff * rho[0] * linv, -ff * rho[1] * linv, -ff * rho[2] * linv };
        if ( pt->body != -1 ) {
            double mom[3];
            mom[0] = ps[1]*Fto[2] - ps[2]*Fto[1];
            mom[1] = ps[2]*Fto[0] - ps[0]*Fto[2];
            mom[2] = ps[0]*Fto[1] - ps[1]*Fto[0];
            for ( int k = 0 ; k < 3 ; k++ ) { elemForceBase[ pt->body ][k] += mom[k]; elemForceBase[ pt->body ][3+k] += Fto[k]; }
        }
        if ( pf->body != -1 ) {
            double Ffrom[3] = { -Fto[0], -Fto[1], -Fto[2] };
            double mom[3];
            mom[0] = pp[1]*Ffrom[2] - pp[2]*Ffrom[1];
            mom[1] = pp[2]*Ffrom[0] - pp[0]*Ffrom[2];
            mom[2] = pp[0]*Ffrom[1] - pp[1]*Ffrom[0];
            for ( int k = 0 ; k < 3 ; k++ ) { elemForceBase[ pf->body ][k] += mom[k]; elemForceBase[ pf->body ][3+k] += Ffrom[k]; }
        }
    }

    // Fold elemForceBase into fvp exactly like fext (i_X*_0 = (Xa^T)^-1),
    // but always on -- elements are not opt-in per body -- and done here,
    // after every body's Xa is known, rather than inside the loop above.
    for ( int i = 0 ; i < nb ; i++ ) {
        int hasForce = 0;
        for ( int k = 0 ; k < 6 ; k++ ) if ( elemForceBase[i][k] != 0 ) hasForce = 1;
        if ( !hasForce ) continue;
        double XaT[36], bvec[6], fb[6];
        for ( int r = 0 ; r < 6 ; r++ )
        for ( int c = 0 ; c < 6 ; c++ ) XaT[ r*6 + c ] = Xa[i][c][r];
        for ( int k = 0 ; k < 6 ; k++ ) bvec[k] = elemForceBase[i][k];
        linsolve( 6, XaT, bvec, fb );
        for ( int k = 0 ; k < 6 ; k++ ) fvp[i][k] -= fb[k];
    }

    // Inward pass: C block = S_i^T fvp_i; fvp_parent += Xup^T fvp_i.
    static double C[ MODEL_MAX_V ];
    for ( int i = nb - 1 ; i >= 0 ; i-- ) {
        int par = mdl->body[i].parent;
        int vi = mdl->vidx[ i ];
        for ( int k = 0 ; k < ndof[i] ; k++ ) {
            double s = 0;
            for ( int r = 0 ; r < 6 ; r++ ) s += S[i][r][k] * fvp[i][r];
            C[ vi + k ] = s;
        }
        if ( par != -1 ) {
            double up[6];
            m6_tmulvec( Xup[i], fvp[i], up );
            for ( int k = 0 ; k < 6 ; k++ ) fvp[par][k] += up[k];
        }
    }

    // CRBA: composite inertias, then H blocks.
    double IC[ MODEL_MAX_BODIES ][6][6];
    for ( int i = 0 ; i < nb ; i++ )
        for ( int r = 0 ; r < 6 ; r++ )
        for ( int c = 0 ; c < 6 ; c++ ) IC[i][r][c] = I6[i][r][c];

    static double H[ MODEL_MAX_V ][ MODEL_MAX_V ];
    for ( int r = 0 ; r < nv ; r++ )
    for ( int c = 0 ; c < nv ; c++ ) H[r][c] = 0;

    for ( int i = nb - 1 ; i >= 0 ; i-- ) {
        int par = mdl->body[i].parent;
        int vi = mdl->vidx[ i ];
        if ( par != -1 ) {
            double folded[6][6];
            m6_congruence( Xup[i], IC[i], folded );
            for ( int r = 0 ; r < 6 ; r++ )
            for ( int c = 0 ; c < 6 ; c++ ) IC[par][r][c] += folded[r][c];
        }

        // F = IC * S_i, one 6-vector per DOF of joint i.
        double F[6][6];
        for ( int k = 0 ; k < ndof[i] ; k++ ) {
            double col[6];
            for ( int r = 0 ; r < 6 ; r++ ) col[r] = S[i][r][k];
            m6_mulvec( IC[i], col, F[k] );
        }
        // Diagonal block: H[vi+a][vi+b] = S_a . F_b.
        for ( int a = 0 ; a < ndof[i] ; a++ )
        for ( int bcol = 0 ; bcol < ndof[i] ; bcol++ ) {
            double s = 0;
            for ( int r = 0 ; r < 6 ; r++ ) s += S[i][r][a] * F[bcol][r];
            H[ vi + a ][ vi + bcol ] = s;
        }

        // Ancestor walk, propagating all of joint i's F columns together.
        int j = i;
        while ( mdl->body[j].parent != -1 ) {
            for ( int k = 0 ; k < ndof[i] ; k++ ) {
                double Fup[6];
                m6_tmulvec( Xup[j], F[k], Fup );
                for ( int r = 0 ; r < 6 ; r++ ) F[k][r] = Fup[r];
            }
            j = mdl->body[j].parent;
            int vj = mdl->vidx[ j ];
            for ( int a = 0 ; a < ndof[i] ; a++ )
            for ( int bcol = 0 ; bcol < ndof[j] ; bcol++ ) {
                double s = 0;
                for ( int r = 0 ; r < 6 ; r++ ) s += F[a][r] * S[j][r][bcol];
                H[ vi + a ][ vj + bcol ] = s;
                H[ vj + bcol ][ vi + a ] = s;
            }
        }
    }

    // Solve H qdd0 = tau - C (the tree solution; final answer if no loops).
    static double A[ MODEL_MAX_V * MODEL_MAX_V ];
    double rhs[ MODEL_MAX_V ];
    for ( int r = 0 ; r < nv ; r++ ) {
        rhs[r] = tau[r] - C[r] + elemQ[r];
        for ( int c = 0 ; c < nv ; c++ ) A[ r*nv + c ] = H[r][c];
    }
    linsolve( nv, A, rhs, qdd_out );

    if ( mdl->nloop == 0 ) return;

    // ---- Kinematic loops: RBDA Table 8.4, everything in base coordinates. --
    int nc = mdl->nc;
    double alpha = 1.0 / mdl->stab, beta = 1.0 / mdl->stab;
    static double K[ MODEL_MAX_NC ][ MODEL_MAX_V ];
    double kfull[ MODEL_MAX_NC ]; // k + kstab, i.e. the RHS of eq 8.35
    for ( int r = 0 ; r < nc ; r++ )
    for ( int c = 0 ; c < nv ; c++ ) K[r][c] = 0;

    for ( int l = 0 ; l < mdl->nloop ; l++ ) {
        const ModelLoop *lp = &mdl->loop[ l ];
        int row = mdl->cidx[ l ];

        // Xp = XP * p(k)_X_0, Xs = XS * s(k)_X_0; XP/XS fold the declared
        // axis into a fixed rotation so the canonical jcalc axis is z.
        double Eax[3][3], XP[6][6], XS[6][6], Xp[6][6], Xs[6][6];
        axisAlignE( lp->axis, Eax );
        plux( Eax, lp->ptree, XP );
        plux( Eax, lp->stree, XS );
        if ( lp->pred == -1 ) { for ( int r = 0 ; r < 6 ; r++ ) for ( int c = 0 ; c < 6 ; c++ ) Xp[r][c] = XP[r][c]; }
        else                  m6_mul( XP, Xa[ lp->pred ], Xp );
        if ( lp->succ == -1 ) { for ( int r = 0 ; r < 6 ; r++ ) for ( int c = 0 ; c < 6 ; c++ ) Xs[r][c] = XS[r][c]; }
        else                  m6_mul( XS, Xa[ lp->succ ], Xs );

        // Body velocity / velocity-product acceleration in base coordinates:
        // v0 = (i_X_0)^-1 v_i, solved rather than inverted. The fixed base
        // has v = 0 and avp = -a_grav (the RNEA seed) -- NOT zero, or the
        // -ag offset fails to cancel in as - ap.
        double vp6[6] = {0,0,0,0,0,0}, vs6[6] = {0,0,0,0,0,0}, ap6[6], as6[6];
        for ( int k = 0 ; k < 6 ; k++ ) { ap6[k] = gravNeg[k]; as6[k] = gravNeg[k]; }
        if ( lp->pred != -1 ) {
            double M[36], b6[6];
            for ( int r = 0 ; r < 6 ; r++ ) for ( int c = 0 ; c < 6 ; c++ ) M[ r*6 + c ] = Xa[ lp->pred ][r][c];
            for ( int k = 0 ; k < 6 ; k++ ) b6[k] = v[ lp->pred ][k];
            linsolve( 6, M, b6, vp6 );
            for ( int r = 0 ; r < 6 ; r++ ) for ( int c = 0 ; c < 6 ; c++ ) M[ r*6 + c ] = Xa[ lp->pred ][r][c];
            for ( int k = 0 ; k < 6 ; k++ ) b6[k] = avp[ lp->pred ][k];
            linsolve( 6, M, b6, ap6 );
        }
        if ( lp->succ != -1 ) {
            double M[36], b6[6];
            for ( int r = 0 ; r < 6 ; r++ ) for ( int c = 0 ; c < 6 ; c++ ) M[ r*6 + c ] = Xa[ lp->succ ][r][c];
            for ( int k = 0 ; k < 6 ; k++ ) b6[k] = v[ lp->succ ][k];
            linsolve( 6, M, b6, vs6 );
            for ( int r = 0 ; r < 6 ; r++ ) for ( int c = 0 ; c < 6 ; c++ ) M[ r*6 + c ] = Xa[ lp->succ ][r][c];
            for ( int k = 0 ; k < 6 ; k++ ) b6[k] = avp[ lp->succ ][k];
            linsolve( 6, M, b6, as6 );
        }

        double vrel[6];
        for ( int k = 0 ; k < 6 ; k++ ) vrel[k] = vs6[k] - vp6[k];

        int ncl;
        double Tb[6][6], delta[6]; // T columns = constraint force directions, base coords

        if ( lp->jtype == JOINT_DISTANCE ) {
            // Rod constraint |p_s - p_p| = d, carried as delta = (rho.rho -
            // d^2)/2 with rho = p_s - p_p (base coords). The constraint force
            // is a pure force rho along the line through BOTH anchors, so the
            // single spatial force T = (p_p x rho; rho) = (p_s x rho; rho)
            // serves predecessor and successor alike, and delta's gradient is
            // exactly the K row this T produces. T is fixed to neither body,
            // so instead of Table 8.4's vs x vp term the acceleration RHS
            // carries the rod's own direction rate:
            //   deltaddot = T.(as - ap) + rho.(ws x psdot - wp x ppdot)
            //             + rhodot.rhodot                    (+ K qdd terms)
            // with pdot = v0 + w x p the anchor point velocities.
            ncl = 1;
            double pp[3], ps[3], rho[3];
            xformOrigin( Xp, pp );
            xformOrigin( Xs, ps );
            for ( int k = 0 ; k < 3 ; k++ ) rho[k] = ps[k] - pp[k];

            Tb[0][0] = pp[1]*rho[2] - pp[2]*rho[1];
            Tb[1][0] = pp[2]*rho[0] - pp[0]*rho[2];
            Tb[2][0] = pp[0]*rho[1] - pp[1]*rho[0];
            for ( int k = 0 ; k < 3 ; k++ ) Tb[3+k][0] = rho[k];

            delta[0] = 0.5 * ( rho[0]*rho[0] + rho[1]*rho[1] + rho[2]*rho[2]
                             - lp->length * lp->length );

            double ppd[3], psd[3], rhod[3], wxpd[3];
            ppd[0] = vp6[3] + vp6[1]*pp[2] - vp6[2]*pp[1];
            ppd[1] = vp6[4] + vp6[2]*pp[0] - vp6[0]*pp[2];
            ppd[2] = vp6[5] + vp6[0]*pp[1] - vp6[1]*pp[0];
            psd[0] = vs6[3] + vs6[1]*ps[2] - vs6[2]*ps[1];
            psd[1] = vs6[4] + vs6[2]*ps[0] - vs6[0]*ps[2];
            psd[2] = vs6[5] + vs6[0]*ps[1] - vs6[1]*ps[0];
            for ( int k = 0 ; k < 3 ; k++ ) rhod[k] = psd[k] - ppd[k];
            wxpd[0] = vs6[1]*psd[2] - vs6[2]*psd[1] - ( vp6[1]*ppd[2] - vp6[2]*ppd[1] );
            wxpd[1] = vs6[2]*psd[0] - vs6[0]*psd[2] - ( vp6[2]*ppd[0] - vp6[0]*ppd[2] );
            wxpd[2] = vs6[0]*psd[1] - vs6[1]*psd[0] - ( vp6[0]*ppd[1] - vp6[1]*ppd[0] );

            double tcol[6], dAvp[6];
            for ( int r = 0 ; r < 6 ; r++ ) { tcol[r] = Tb[r][0]; dAvp[r] = as6[r] - ap6[r]; }
            double kvel = rho[0]*wxpd[0] + rho[1]*wxpd[1] + rho[2]*wxpd[2]
                        + rhod[0]*rhod[0] + rhod[1]*rhod[1] + rhod[2]*rhod[2];
            double kstab = 2.0 * alpha * dot6( tcol, vrel ) + beta * beta * delta[0];
            kfull[ row ] = -( dot6( tcol, dAvp ) + kvel ) - kstab;
        } else {
            // sXp of the loop joint, then [delta, T] = jcalc(type, Xs * Xp^-1).
            double XpInv[6][6], sXp[6][6], Tloc[6][6];
            {
                double M[36], col[6], e[6];
                for ( int c = 0 ; c < 6 ; c++ ) {
                    for ( int r = 0 ; r < 36 ; r++ ) M[r] = 0;
                    for ( int r = 0 ; r < 6 ; r++ ) for ( int cc = 0 ; cc < 6 ; cc++ ) M[ r*6 + cc ] = Xp[r][cc];
                    for ( int k = 0 ; k < 6 ; k++ ) e[k] = ( k == c );
                    linsolve( 6, M, e, col );
                    for ( int r = 0 ; r < 6 ; r++ ) XpInv[r][c] = col[r];
                }
            }
            m6_mul( Xs, XpInv, sXp );
            ncl = loopJcalc( lp->jtype, sXp, Tloc, delta );

            // T to base coordinates: force vectors move with X^T, using the
            // frame this joint type's T is fixed to.
            for ( int c = 0 ; c < ncl ; c++ ) {
                double tcol[6], tb[6];
                for ( int r = 0 ; r < 6 ; r++ ) tcol[r] = Tloc[r][c];
                m6_tmulvec( loopTInPred( lp->jtype ) ? Xp : Xs, tcol, tb );
                for ( int r = 0 ; r < 6 ; r++ ) Tb[r][c] = tb[r];
            }

            // kfull = -T^T (as - ap + vs x vp) - kstab_local (Table 8.4; its
            // kstab local variable is the NEGATIVE of eq 8.22's kstab, and
            // this difference IS k + kstab of eq 8.35).
            double crmVs[6][6], vsxvp[6], resid[6];
            crm( vs6, crmVs );
            m6_mulvec( crmVs, vp6, vsxvp );
            for ( int k = 0 ; k < 6 ; k++ ) resid[k] = as6[k] - ap6[k] + vsxvp[k];
            for ( int c = 0 ; c < ncl ; c++ ) {
                double tcol[6];
                for ( int r = 0 ; r < 6 ; r++ ) tcol[r] = Tb[r][c];
                double kstab = 2.0 * alpha * dot6( tcol, vrel ) + beta * beta * delta[c];
                kfull[ row + c ] = -dot6( tcol, resid ) - kstab;
            }
        }
        if ( delta_out )
            for ( int c = 0 ; c < ncl ; c++ ) delta_out[ row + c ] = delta[ c ];

        // K blocks: walk both paths to the common ancestor; -T^T 0X_i S_i on
        // the predecessor side, + on the successor side (Table 8.4). Parents
        // always have smaller indices and the base is -1, so stepping the
        // larger index meets at the nearest common ancestor.
        int i = lp->pred, j = lp->succ;
        while ( i != j ) {
            int onPred = i > j;
            int bidx = onPred ? i : j;
            double sign = onPred ? -1.0 : 1.0;
            for ( int col = 0 ; col < ndof[ bidx ] ; col++ ) {
                double scol[6], s0[6], M[36], b6[6];
                for ( int r = 0 ; r < 6 ; r++ ) { scol[r] = S[ bidx ][r][col]; }
                for ( int r = 0 ; r < 6 ; r++ ) for ( int cc = 0 ; cc < 6 ; cc++ ) M[ r*6 + cc ] = Xa[ bidx ][r][cc];
                for ( int k = 0 ; k < 6 ; k++ ) b6[k] = scol[k];
                linsolve( 6, M, b6, s0 );
                for ( int c = 0 ; c < ncl ; c++ ) {
                    double tcol[6];
                    for ( int r = 0 ; r < 6 ; r++ ) tcol[r] = Tb[r][c];
                    K[ row + c ][ mdl->vidx[ bidx ] + col ] += sign * dot6( tcol, s0 );
                }
            }
            if ( onPred ) i = mdl->body[ i ].parent;
            else          j = mdl->body[ j ].parent;
        }
    }

    if ( K_out )
        for ( int r = 0 ; r < nc ; r++ )
        for ( int c = 0 ; c < nv ; c++ ) K_out[ r*nv + c ] = K[r][c];

    // Method 2 (RBDA sec 8.5): Y = H^-1 K^T column by column, A_c = K Y,
    // b = kfull - K qdd0, solve A_c lambda = b, qdd = qdd0 + Y lambda.
    static double Y[ MODEL_MAX_V ][ MODEL_MAX_NC ];
    for ( int c = 0 ; c < nc ; c++ ) {
        double col[ MODEL_MAX_V ], sol[ MODEL_MAX_V ];
        for ( int r = 0 ; r < nv ; r++ ) {
            col[r] = K[c][r];
            for ( int cc = 0 ; cc < nv ; cc++ ) A[ r*nv + cc ] = H[r][cc];
        }
        linsolve( nv, A, col, sol );
        for ( int r = 0 ; r < nv ; r++ ) Y[r][c] = sol[r];
    }
    double Ac[ MODEL_MAX_NC * MODEL_MAX_NC ], bc[ MODEL_MAX_NC ], lambda[ MODEL_MAX_NC ];
    for ( int r = 0 ; r < nc ; r++ ) {
        bc[r] = kfull[r];
        for ( int c = 0 ; c < nv ; c++ ) bc[r] -= K[r][c] * qdd_out[c];
        for ( int c = 0 ; c < nc ; c++ ) {
            double s = 0;
            for ( int k = 0 ; k < nv ; k++ ) s += K[r][k] * Y[k][c];
            Ac[ r*nc + c ] = s;
        }
    }
    linsolve( nc, Ac, bc, lambda );
    for ( int r = 0 ; r < nv ; r++ ) {
        double s = 0;
        for ( int c = 0 ; c < nc ; c++ ) s += Y[r][c] * lambda[c];
        qdd_out[r] += s;
    }
}

static void fd_dense( const Model *mdl, const double *q, const double *qd,
                      const double *tau, const double *fext, const double *act, double *qdd_out ) {
    fd_dense_full( mdl, q, qd, tau, fext, act, qdd_out, NULL, NULL );
}

// ---------------------------------------------------------------------------
// Deterministic pseudo-random state generator, shared by the harnesses.
// ---------------------------------------------------------------------------
static unsigned long long rngState = 0x853c49e6748fea9bULL;
static double frand( void ) { // uniform in [-1, 1]
    rngState = rngState * 6364136223846793005ULL + 1442695040888963407ULL;
    return ( ( rngState >> 17 ) & 0xFFFFFF ) / (double)0x800000 - 1.0;
}

// ---------------------------------------------------------------------------
// Closure projection for loop models, shared by the validation and Chrono
// harnesses. Random states do not satisfy the loop constraints; these pull a
// state onto the closure manifold so constrained dynamics is evaluated at a
// mechanism that is actually assembled. fd_dense_full exposes delta (the
// position error) and K = d(delta)/dq to first order, so Newton with
// minimal-norm updates dq = -K^T (K K^T)^-1 delta converges to the manifold,
// and the same projector applied to K qd makes the velocity consistent.
//
// K is indexed by VELOCITY coordinate, so the Newton update lives in the
// velocity tangent space (nv slots) and must be mapped into q (nq slots).
// For single-DOF joints that is plain addition; for a floating body the six
// slots are a body-frame twist [w; v], applied through the quaternion
// tangent map -- the same rates integrateQ uses:
//     Q += 0.5 * Q x (0, w)   (then normalize)
//     r += R(Q) * v
// ---------------------------------------------------------------------------
static void tangentStep( const Model *mdl, double *q, const double *dqv, double scale ) {
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        int qi = mdl->qidx[ i ], vi = mdl->vidx[ i ];
        if ( mdl->body[ i ].jtype != JOINT_FLOATING ) {
            q[ qi ] += scale * dqv[ vi ];
            continue;
        }
        double w = q[qi], x = q[qi+1], y = q[qi+2], z = q[qi+3];
        double n = 1.0 / sqrt( w*w + x*x + y*y + z*z );
        w *= n; x *= n; y *= n; z *= n;
        double wx = scale * dqv[vi], wy = scale * dqv[vi+1], wz = scale * dqv[vi+2];
        double qw = w - 0.5 * ( x*wx + y*wy + z*wz );
        double qx = x + 0.5 * ( w*wx + y*wz - z*wy );
        double qy = y + 0.5 * ( w*wy + z*wx - x*wz );
        double qz = z + 0.5 * ( w*wz + x*wy - y*wx );
        double nn = 1.0 / sqrt( qw*qw + qx*qx + qy*qy + qz*qz );
        q[qi] = qw*nn; q[qi+1] = qx*nn; q[qi+2] = qy*nn; q[qi+3] = qz*nn;
        double R[3][3];
        R[0][0] = 2*( w*w + x*x - 0.5 ); R[0][1] = 2*( x*y - w*z );       R[0][2] = 2*( x*z + w*y );
        R[1][0] = 2*( x*y + w*z );       R[1][1] = 2*( w*w + y*y - 0.5 ); R[1][2] = 2*( y*z - w*x );
        R[2][0] = 2*( x*z - w*y );       R[2][1] = 2*( y*z + w*x );       R[2][2] = 2*( w*w + z*z - 0.5 );
        double vv[3] = { scale * dqv[vi+3], scale * dqv[vi+4], scale * dqv[vi+5] };
        for ( int r = 0 ; r < 3 ; r++ )
            q[ qi + 4 + r ] += R[r][0]*vv[0] + R[r][1]*vv[1] + R[r][2]*vv[2];
    }
}

// Damped Gauss-Newton: the raw step is always a descent direction of
// |delta|^2, but from a badly torn state a full step can overshoot the
// linearization (rotations especially), so the tangent step norm is capped
// and the step is backtracked until |delta| actually decreases. Returns
// nonzero if it stalls -- which happens at genuine local minima of |delta|^2
// when the target is outside the mechanism's reachable workspace, so callers
// should treat failure as "state unreachable", not as a numerical mishap.
static int projectClosure( const Model *mdl, double *q ) {
    if ( !mdl->nloop ) return 0;
    static double K[ MODEL_MAX_NC * MODEL_MAX_V ];
    double delta[ MODEL_MAX_NC ], qdd[ MODEL_MAX_V ];
    double zeroV[ MODEL_MAX_V ] = { 0 }, zeroF[ MODEL_MAX_BODIES * 6 ] = { 0 };
    int nc = mdl->nc, nv = mdl->nv;

    // Geometry-only call: K/delta don't depend on forces at all, so a zero
    // act (reusing the already-zeroed zeroV buffer) is exact, not a stand-in.
    fd_dense_full( mdl, q, zeroV, zeroV, zeroF, zeroV, qdd, K, delta );
    double n2 = 0;
    for ( int c = 0 ; c < nc ; c++ ) n2 += delta[c] * delta[c];

    for ( int iter = 0 ; iter < 200 ; iter++ ) {
        if ( sqrt( n2 ) < 1e-12 ) return 0;

        double G[ MODEL_MAX_NC * MODEL_MAX_NC ], rhs[ MODEL_MAX_NC ], lam[ MODEL_MAX_NC ];
        for ( int r = 0 ; r < nc ; r++ ) {
            rhs[r] = delta[r];
            for ( int c = 0 ; c < nc ; c++ ) {
                double s = 0;
                for ( int k = 0 ; k < nv ; k++ ) s += K[ r*nv + k ] * K[ c*nv + k ];
                G[ r*nc + c ] = s;
            }
        }
        linsolve( nc, G, rhs, lam );
        double dqv[ MODEL_MAX_V ], step2 = 0;
        for ( int k = 0 ; k < nv ; k++ ) {
            double s = 0;
            for ( int c = 0 ; c < nc ; c++ ) s += K[ c*nv + k ] * lam[c];
            dqv[k] = -s;
            step2 += s * s;
        }
        double scale = sqrt( step2 ) > 0.5 ? 0.5 / sqrt( step2 ) : 1.0;

        int accepted = 0;
        for ( int bt = 0 ; bt < 30 ; bt++ ) {
            double qtry[ MODEL_MAX_Q ];
            for ( int k = 0 ; k < mdl->nq ; k++ ) qtry[k] = q[k];
            tangentStep( mdl, qtry, dqv, scale );
            double dtry[ MODEL_MAX_NC ], n2try = 0;
            fd_dense_full( mdl, qtry, zeroV, zeroV, zeroF, zeroV, qdd, K, dtry );
            for ( int c = 0 ; c < nc ; c++ ) n2try += dtry[c] * dtry[c];
            if ( n2try < n2 ) {
                for ( int k = 0 ; k < mdl->nq ; k++ ) q[k] = qtry[k];
                for ( int c = 0 ; c < nc ; c++ ) delta[c] = dtry[c];
                n2 = n2try;
                accepted = 1;
                break;
            }
            scale *= 0.5;
        }
        if ( !accepted ) return 1;
    }
    return sqrt( n2 ) < 1e-12 ? 0 : 1;
}

// ---------------------------------------------------------------------------
// Random-walk state generation on the closure manifold, shared by the
// validation and Chrono harnesses (both need many varied assembled states).
// Fully random q tears the loops open by meters, and at such states the
// method-2 matrix A = K H^-1 K^T can be almost arbitrarily ill-conditioned
// -- the quartercar model hits cond(A) ~ 1e10 within 1000 draws, which turns
// roundoff into 1e-7 noise between implementations that actually agree.
// Near the manifold cond(A) stays that of the physical mechanism (~1e4
// there). Direct projection of wild states is also unreliable: Gauss-Newton
// stalls at local minima of the closure residual whenever the random target
// is outside the reachable workspace. The walk avoids both problems: seed
// once, then take a moderate random tangent step each trial and project
// back -- every projection starts near the manifold, so it is local.
// ---------------------------------------------------------------------------

// Seeds the walk: q = 0 with identity quaternions first (most models are
// authored to close there), then random draws. The q = 0 attempt can
// legitimately fail for a valid model -- fourbar's q = 0 is the
// all-links-collinear fold, a stationary point of the closure residual.
// Returns nonzero if no seed projects, i.e. the mechanism never assembles.
static int manifoldSeed( const Model *mdl, double *qman ) {
    for ( int a = 0 ; a < 100 ; a++ ) {
        if ( a == 0 ) {
            for ( int k = 0 ; k < mdl->nq ; k++ ) qman[k] = 0;
            for ( int i = 0 ; i < mdl->nbod ; i++ )
                if ( mdl->body[i].jtype == JOINT_FLOATING ) qman[ mdl->qidx[i] ] = 1.0;
        } else {
            for ( int k = 0 ; k < mdl->nq ; k++ ) qman[k] = 3.0 * frand();
        }
        if ( !projectClosure( mdl, qman ) ) return 0;
    }
    return 1;
}

// One walk step: random tangent step, projected back. A step can leave the
// reachable workspace (a genuine local minimum of the closure residual, not
// a solver failure); the walk then stays at the previous on-manifold state,
// so qman is always valid on return.
static void manifoldStep( const Model *mdl, double *qman ) {
    double dqv[ MODEL_MAX_V ], qprev[ MODEL_MAX_Q ];
    for ( int k = 0 ; k < mdl->nv ; k++ ) dqv[k] = 0.5 * frand();
    for ( int k = 0 ; k < mdl->nq ; k++ ) qprev[k] = qman[k];
    tangentStep( mdl, qman, dqv, 1.0 );
    if ( projectClosure( mdl, qman ) )
        for ( int k = 0 ; k < mdl->nq ; k++ ) qman[k] = qprev[k];
}

// Velocity projection needs no tangent map: it operates on qd directly,
// which is already indexed by velocity coordinate for every joint type.
static void projectVelocity( const Model *mdl, const double *q, double *qd ) {
    if ( !mdl->nloop ) return;
    static double K[ MODEL_MAX_NC * MODEL_MAX_V ];
    double qdd[ MODEL_MAX_V ];
    double zeroV[ MODEL_MAX_V ] = { 0 }, zeroF[ MODEL_MAX_BODIES * 6 ] = { 0 };
    int nc = mdl->nc, nv = mdl->nv;
    fd_dense_full( mdl, q, zeroV, zeroV, zeroF, zeroV, qdd, K, NULL );

    double G[ MODEL_MAX_NC * MODEL_MAX_NC ], rhs[ MODEL_MAX_NC ], lam[ MODEL_MAX_NC ];
    for ( int r = 0 ; r < nc ; r++ ) {
        rhs[r] = 0;
        for ( int k = 0 ; k < nv ; k++ ) rhs[r] += K[ r*nv + k ] * qd[k];
        for ( int c = 0 ; c < nc ; c++ ) {
            double s = 0;
            for ( int k = 0 ; k < nv ; k++ ) s += K[ r*nv + k ] * K[ c*nv + k ];
            G[ r*nc + c ] = s;
        }
    }
    linsolve( nc, G, rhs, lam );
    for ( int k = 0 ; k < nv ; k++ ) {
        double s = 0;
        for ( int c = 0 ; c < nc ; c++ ) s += K[ c*nv + k ] * lam[c];
        qd[k] -= s;
    }
}
