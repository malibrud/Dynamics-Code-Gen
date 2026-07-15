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

// Plucker motion transform from rotation E and translation r (RBDA eq 2.24
// with our generator's parametrization): X = [ E, 0; -E*rx, E ].
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

// Rotation by q about a canonical axis, in the generator's convention
// (XFORM_RZ: E = [ c -s 0; s c 0; 0 0 1 ], and the cyclic analogues).
static void rotAxis( char axis, double q, double E[3][3] ) {
    double c = cos( q ), s = sin( q );
    if ( axis == 'x' ) {
        E[0][0] = 1; E[0][1] = 0; E[0][2] =  0;
        E[1][0] = 0; E[1][1] = c; E[1][2] = -s;
        E[2][0] = 0; E[2][1] = s; E[2][2] =  c;
    } else if ( axis == 'y' ) {
        E[0][0] =  c; E[0][1] = 0; E[0][2] = s;
        E[1][0] =  0; E[1][1] = 1; E[1][2] = 0;
        E[2][0] = -s; E[2][1] = 0; E[2][2] = c;
    } else {
        E[0][0] = c; E[0][1] = -s; E[0][2] = 0;
        E[1][0] = s; E[1][1] =  c; E[1][2] = 0;
        E[2][0] = 0; E[2][1] =  0; E[2][2] = 1;
    }
}

// Joint transform XJ and motion subspace S (6 x ndof, stored as columns of a
// 6x6 array) for a joint at position coordinates qseg (this body's slice of
// the packed q vector). Revolute: rotation about the axis, S = (axis; 0).
// Prismatic: translation along the axis, S = (0; axis). Floating: quaternion
// orientation plus position, S = identity(6) -- a genuine 6-DOF joint with
// no gimbal lock.
static int jointXS( const ModelBody *b, const double *qseg, double XJ[6][6], double S[6][6] ) {
    double a3[3], Eident[3][3] = { {1,0,0},{0,1,0},{0,0,1} }, r0[3] = { 0, 0, 0 };
    for ( int r = 0 ; r < 6 ; r++ )
    for ( int c = 0 ; c < 6 ; c++ ) S[r][c] = 0;

    if ( b->jtype == JOINT_FLOATING ) {
        // Unit quaternion (w,x,y,z), normalized here; rotation matrix in
        // this codebase's E convention (transpose of spatial_v2's rq.m).
        double w = qseg[0], x = qseg[1], y = qseg[2], z = qseg[3];
        double n = 1.0 / sqrt( w*w + x*x + y*y + z*z );
        w *= n; x *= n; y *= n; z *= n;
        double E[3][3];
        E[0][0] = 2.0*( w*w + x*x - 0.5 ); E[0][1] = 2.0*( x*y - w*z );       E[0][2] = 2.0*( x*z + w*y );
        E[1][0] = 2.0*( x*y + w*z );       E[1][1] = 2.0*( w*w + y*y - 0.5 ); E[1][2] = 2.0*( y*z - w*x );
        E[2][0] = 2.0*( x*z - w*y );       E[2][1] = 2.0*( y*z + w*x );       E[2][2] = 2.0*( w*w + z*z - 0.5 );
        plux( E, qseg + 4, XJ );
        for ( int k = 0 ; k < 6 ; k++ ) S[k][k] = 1;
        return 6;
    }

    axisVec( b->axis, a3 );
    if ( b->jtype == JOINT_REVOLUTE ) {
        double E[3][3];
        rotAxis( b->axis, qseg[0], E );
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
// Dense forward dynamics: H*qdd = tau - C, with C from the RNEA (RBDA Table
// 5.1, aJ = 0) and H from the CRBA (RBDA Table 6.4). Joints may have 1 DOF
// (revolute/prismatic) or 6 (floating), so H is nv x nv with per-joint
// blocks; q is addressed through mdl->qidx and qd/tau/qdd through mdl->vidx.
// The large working arrays are static: fd_dense is not reentrant, but stays
// comfortably off the stack even at MODEL_MAX_V.
// ---------------------------------------------------------------------------
static void fd_dense( const Model *mdl, const double *q, const double *qd,
                      const double *tau, const double *fext, double *qdd_out ) {
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

    // Solve H qdd = tau - C.
    static double A[ MODEL_MAX_V * MODEL_MAX_V ];
    double rhs[ MODEL_MAX_V ];
    for ( int r = 0 ; r < nv ; r++ ) {
        rhs[r] = tau[r] - C[r];
        for ( int c = 0 ; c < nv ; c++ ) A[ r*nv + c ] = H[r][c];
    }
    linsolve( nv, A, rhs, qdd_out );
}

// ---------------------------------------------------------------------------
// Deterministic pseudo-random state generator, shared by the harnesses.
// ---------------------------------------------------------------------------
static unsigned long long rngState = 0x853c49e6748fea9bULL;
static double frand( void ) { // uniform in [-1, 1]
    rngState = rngState * 6364136223846793005ULL + 1442695040888963407ULL;
    return ( ( rngState >> 17 ) & 0xFFFFFF ) / (double)0x800000 - 1.0;
}
