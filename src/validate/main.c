// Validation harness: runs the generated sparse dynamics (build/dynamics.c)
// and the dense 6x6 reference (dense.h) on random states and compares qdd.
// Both take the mechanism from the same model file. Exit code 0 = PASS,
// 1 = FAIL. See dense.h for the reference implementation.

#include <stdio.h>
#include "dense.h"

// The generated implementation under test (compiled from build/dynamics.c).
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *out_qdd );
extern const int fwd_dyn_crba_nq;
extern const int fwd_dyn_crba_nv;

// Numerical rank of the nc x nv constraint matrix K by Gaussian elimination
// with full pivoting. The pipeline requires properly constrained loops
// (rank = nc); anything less means redundant constraints and a singular
// K H^-1 K^T, which the generated method-2 solve cannot handle.
static int rankK( int nc, int nv, const double *Kin ) {
    double K[ MODEL_MAX_NC ][ MODEL_MAX_V ];
    for ( int r = 0 ; r < nc ; r++ )
    for ( int c = 0 ; c < nv ; c++ ) K[r][c] = Kin[ r*nv + c ];

    double maxAbs = 0;
    for ( int r = 0 ; r < nc ; r++ )
    for ( int c = 0 ; c < nv ; c++ ) if ( fabs( K[r][c] ) > maxAbs ) maxAbs = fabs( K[r][c] );
    double tol = 1e-8 * ( maxAbs > 0 ? maxAbs : 1.0 );

    int rank = 0;
    for ( int step = 0 ; step < nc ; step++ ) {
        int pr = -1, pc = -1;
        double best = tol;
        for ( int r = step ; r < nc ; r++ )
        for ( int c = 0 ; c < nv ; c++ ) if ( fabs( K[r][c] ) > best ) { best = fabs( K[r][c] ); pr = r; pc = c; }
        if ( pr < 0 ) break;
        rank++;
        for ( int c = 0 ; c < nv ; c++ ) { double t = K[step][c]; K[step][c] = K[pr][c]; K[pr][c] = t; }
        for ( int r = step + 1 ; r < nc ; r++ ) {
            double f = K[r][pc] / K[step][pc];
            for ( int c = 0 ; c < nv ; c++ ) K[r][c] -= f * K[step][c];
        }
    }
    return rank;
}

int main( int argc, char **argv ) {
    const char *path = argc > 1 ? argv[ 1 ] : "../models/threelink.mech";
    Model mdl;
    char err[ 256 ];
    if ( modelLoad( path, &mdl, err, sizeof err ) ) {
        fprintf( stderr, "validate: %s\n", err );
        return 1;
    }
    const char *unsup = modelUnsupported( &mdl );
    if ( unsup ) {
        fprintf( stderr, "validate: %s\n", unsup );
        return 1;
    }
    if ( fwd_dyn_crba_nq != mdl.nq || fwd_dyn_crba_nv != mdl.nv ) {
        fprintf( stderr, "validate: dynamics.c was generated for nq=%d nv=%d but '%s' has nq=%d nv=%d -- regenerate\n",
                 fwd_dyn_crba_nq, fwd_dyn_crba_nv, path, mdl.nq, mdl.nv );
        return 1;
    }
    int nb = mdl.nbod, nq = mdl.nq, nv = mdl.nv;

    const int NTRIALS = 1000;
    const double TOL = 1e-9;
    double maxErr = 0;
    int worstTrial = -1;

    for ( int t = 0 ; t < NTRIALS ; t++ ) {
        double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ], tau[ MODEL_MAX_V ];
        double fext[ MODEL_MAX_BODIES * 6 ];

        // First trials pin down special cases; the rest are fully random.
        // Loop models always draw random q: a staged q = 0 can sit at a
        // variable-mobility singularity (RBDA sec 8.10, e.g. a four-bar with
        // all links collinear), where rank(K) legitimately drops and the
        // method-2 solve is undefined. Random configurations hit such
        // singular sets with probability zero, so the rank check below then
        // really does detect structural overconstraint only.
        for ( int k = 0 ; k < nq ; k++ ) q[k]   = ( t == 0 && !mdl.nloop ) ? 0 : 3.0 * frand();
        for ( int k = 0 ; k < nv ; k++ ) qd[k]  = ( t <= 1 ) ? 0 : 2.0 * frand();
        for ( int k = 0 ; k < nv ; k++ ) tau[k] = ( t <= 2 ) ? 0 : 5.0 * frand();
        for ( int k = 0 ; k < nb * 6 ; k++ )
            fext[k] = ( t <= 3 ) ? 0 : 5.0 * frand();

        // An all-zero quaternion is not a pose; the at-rest trial needs the
        // identity quaternion for every floating body. Random trials keep
        // whatever they drew -- both implementations normalize identically.
        if ( t == 0 ) {
            for ( int i = 0 ; i < nb ; i++ )
                if ( mdl.body[i].jtype == JOINT_FLOATING ) q[ mdl.qidx[i] ] = 1.0;
        }

        double qddGen[ MODEL_MAX_V ], qddRef[ MODEL_MAX_V ];
        static double K[ MODEL_MAX_NC * MODEL_MAX_V ];
        fwd_dyn_crba( q, qd, tau, fext, qddGen );
        fd_dense_full( &mdl, q, qd, tau, fext, qddRef, mdl.nloop ? K : NULL, NULL );

        if ( mdl.nloop ) {
            int r = rankK( mdl.nc, nv, K );
            if ( r < mdl.nc ) {
                fprintf( stderr, "validate: FAIL: rank(K) = %d < nc = %d at trial %d -- the loop joints\n"
                                 "impose redundant constraints (RBDA sec 8.10). Close the loop with a\n"
                                 "joint imposing fewer constraints (e.g. point-on-line for planar loops).\n",
                         r, mdl.nc, t );
                return 1;
            }
        }

        for ( int k = 0 ; k < nv ; k++ ) {
            // Scaled error: absolute for O(1) values, relative for large
            // ones, so ill-conditioned states (near-singular H) measure
            // agreement rather than conditioning. A real bug still shows
            // up as a large relative error.
            double e = fabs( qddGen[k] - qddRef[k] ) / fmax( 1.0, fabs( qddRef[k] ) );
            if ( e > maxErr ) { maxErr = e; worstTrial = t; }
        }
    }

    printf( "validate: %s: %d trials, max scaled |qdd_generated - qdd_dense| = %.3e (worst at trial %d)\n",
            path, NTRIALS, maxErr, worstTrial );
    if ( maxErr < TOL ) {
        printf( "validate: PASS\n" );
        return 0;
    }
    printf( "validate: FAIL (tolerance %.1e)\n", TOL );
    return 1;
}
