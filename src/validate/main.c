// Dense-matrix validation reference for the generated sparse dynamics.
//
// This program implements the same forward-dynamics computation as the code
// generator (src/codegen -> build/dynamics.c) but with full 6x6 spatial
// matrices and 6-vectors, straight from the equations in RBDA:
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
// It runs both implementations on random states and compares qdd.
// Exit code 0 = PASS, 1 = FAIL.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NBOD 3

// The generated implementation under test (compiled from build/dynamics.c).
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *out_qdd );

// ---------------------------------------------------------------------------
// Model description. MUST MATCH src/codegen/main.cpp: a planar chain of NBOD
// links, every joint revolute about z, every tree transform a translation of
// 1 along x, parent(i) = i-1.
// ---------------------------------------------------------------------------
static const double masses[ NBOD ] = { 1.0, 1.0, 1.0 };
static const double cgs[ NBOD ][ 3 ] = {
    { 0.5, 0.0, 0.0 },
    { 1.0, 0.5, 0.0 },
    { 1.5, 1.0, 0.0 },
};
// Rotational inertia about the CG: xx, yy, zz, xy, yz, zx (actual matrix
// entries, off-diagonal minus signs already absorbed -- same convention as
// the generator's Inertia struct).
static const double icm[ NBOD ][ 6 ] = {
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
};
static const double linkLen = 1.0;
// -a_grav, the seed of the RNEA acceleration recursion (the generator's
// sys.grav): physical gravity is (0,0,-9.81), so this is +9.81 in linear z.
static const double gravNeg[ 6 ] = { 0, 0, 0, 0, 0, 9.81 };

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

// Rotation about z by q, in the generator's XFORM_RZ convention:
// E = [ c -s 0; s c 0; 0 0 1 ].
static void rotz( double q, double E[3][3] ) {
    double c = cos( q ), s = sin( q );
    E[0][0] = c; E[0][1] = -s; E[0][2] = 0;
    E[1][0] = s; E[1][1] =  c; E[1][2] = 0;
    E[2][0] = 0; E[2][1] =  0; E[2][2] = 1;
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
// 5.1, aJ = 0) and H from the CRBA (RBDA Table 6.4). parent(i) = i-1.
// ---------------------------------------------------------------------------
static void fd_dense( const double *q, const double *qd, const double *tau,
                      const double *fext, double *qdd_out ) {
    double Xup[ NBOD ][6][6]; // i_X_lambda(i)
    double Xa[ NBOD ][6][6];  // i_X_0
    double I6[ NBOD ][6][6];
    double v[ NBOD ][6], avp[ NBOD ][6], fvp[ NBOD ][6];

    for ( int i = 0 ; i < NBOD ; i++ ) {
        mcI( masses[i], cgs[i], icm[i], I6[i] );

        // Xup = XJ * XT, exactly as the generator composes
        // Xlink = xJoint o xTree.
        double E[3][3], Eident[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
        double rT[3] = { linkLen, 0, 0 }, r0[3] = { 0, 0, 0 };
        double XJ[6][6], XT[6][6];
        rotz( q[i], E );
        plux( E, r0, XJ );
        plux( Eident, rT, XT );
        m6_mul( XJ, XT, Xup[i] );

        if ( i == 0 ) {
            for ( int r = 0 ; r < 6 ; r++ )
            for ( int c = 0 ; c < 6 ; c++ ) Xa[i][r][c] = Xup[i][r][c];
        } else {
            m6_mul( Xup[i], Xa[i-1], Xa[i] );
        }

        // v_i = Xup * v_parent + vJ,  vJ = S * qd_i = (0,0,qd,0,0,0)
        double vJ[6] = { 0, 0, qd[i], 0, 0, 0 };
        if ( i == 0 ) {
            for ( int k = 0 ; k < 6 ; k++ ) v[i][k] = vJ[k];
        } else {
            m6_mulvec( Xup[i], v[i-1], v[i] );
            for ( int k = 0 ; k < 6 ; k++ ) v[i][k] += vJ[k];
        }

        // avp_i = Xup * avp_parent + crm(v_i)*vJ   (avp_parent(0) = -a_grav)
        double crmV[6][6], vCrossVJ[6];
        m6_mulvec( Xup[i], i == 0 ? gravNeg : avp[i-1], avp[i] );
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

        double XaT[36], b[6], fb[6];
        for ( int r = 0 ; r < 6 ; r++ )
        for ( int c = 0 ; c < 6 ; c++ ) XaT[ r*6 + c ] = Xa[i][c][r];
        for ( int k = 0 ; k < 6 ; k++ ) b[k] = fext[ i*6 + k ];
        linsolve( 6, XaT, b, fb );
        for ( int k = 0 ; k < 6 ; k++ ) fvp[i][k] -= fb[k];
    }

    // Inward pass: C_i = S^T fvp_i; fvp_parent += Xup^T fvp_i.
    double C[ NBOD ];
    for ( int i = NBOD - 1 ; i >= 0 ; i-- ) {
        C[i] = fvp[i][2];
        if ( i > 0 ) {
            double up[6];
            m6_tmulvec( Xup[i], fvp[i], up );
            for ( int k = 0 ; k < 6 ; k++ ) fvp[i-1][k] += up[k];
        }
    }

    // CRBA: composite inertias, then H.
    double IC[ NBOD ][6][6];
    for ( int i = 0 ; i < NBOD ; i++ )
        for ( int r = 0 ; r < 6 ; r++ )
        for ( int c = 0 ; c < 6 ; c++ ) IC[i][r][c] = I6[i][r][c];

    double H[ NBOD ][ NBOD ] = { 0 };
    for ( int i = NBOD - 1 ; i >= 0 ; i-- ) {
        if ( i > 0 ) {
            double folded[6][6];
            m6_congruence( Xup[i], IC[i], folded );
            for ( int r = 0 ; r < 6 ; r++ )
            for ( int c = 0 ; c < 6 ; c++ ) IC[i-1][r][c] += folded[r][c];
        }
        double F[6];
        for ( int k = 0 ; k < 6 ; k++ ) F[k] = IC[i][k][2]; // IC * S, S = e_z
        H[i][i] = F[2];
        for ( int j = i ; j > 0 ; ) {
            double Fup[6];
            m6_tmulvec( Xup[j], F, Fup );
            for ( int k = 0 ; k < 6 ; k++ ) F[k] = Fup[k];
            j = j - 1;
            H[i][j] = F[2]; // F . S_j
            H[j][i] = H[i][j];
        }
    }

    // Solve H qdd = tau - C.
    double A[ NBOD * NBOD ], rhs[ NBOD ];
    for ( int i = 0 ; i < NBOD ; i++ ) {
        rhs[i] = tau[i] - C[i];
        for ( int j = 0 ; j < NBOD ; j++ ) A[ i*NBOD + j ] = H[i][j];
    }
    linsolve( NBOD, A, rhs, qdd_out );
}

// ---------------------------------------------------------------------------
// Comparison harness.
// ---------------------------------------------------------------------------
static unsigned long long rngState = 0x853c49e6748fea9bULL;
static double frand( void ) { // uniform in [-1, 1], deterministic
    rngState = rngState * 6364136223846793005ULL + 1442695040888963407ULL;
    return ( ( rngState >> 17 ) & 0xFFFFFF ) / (double)0x800000 - 1.0;
}

int main( void ) {
    const int NTRIALS = 1000;
    const double TOL = 1e-9;
    double maxErr = 0;
    int worstTrial = -1;

    for ( int t = 0 ; t < NTRIALS ; t++ ) {
        double q[ NBOD ], qd[ NBOD ], tau[ NBOD ], fext[ NBOD * 6 ];

        // First trials pin down special cases; the rest are fully random.
        for ( int i = 0 ; i < NBOD ; i++ ) {
            q[i]   = ( t == 0 )          ? 0 : 3.0 * frand();
            qd[i]  = ( t <= 1 )          ? 0 : 2.0 * frand();
            tau[i] = ( t <= 2 )          ? 0 : 5.0 * frand();
        }
        for ( int k = 0 ; k < NBOD * 6 ; k++ )
            fext[k] = ( t <= 3 ) ? 0 : 5.0 * frand();

        double qddGen[ NBOD ], qddRef[ NBOD ];
        fwd_dyn_crba( q, qd, tau, fext, qddGen );
        fd_dense( q, qd, tau, fext, qddRef );

        for ( int i = 0 ; i < NBOD ; i++ ) {
            double err = fabs( qddGen[i] - qddRef[i] );
            if ( err > maxErr ) { maxErr = err; worstTrial = t; }
        }
    }

    printf( "validate: %d trials, max |qdd_generated - qdd_dense| = %.3e (worst at trial %d)\n",
            NTRIALS, maxErr, worstTrial );
    if ( maxErr < TOL ) {
        printf( "validate: PASS\n" );
        return 0;
    }
    printf( "validate: FAIL (tolerance %.1e)\n", TOL );
    return 1;
}
