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
        for ( int k = 0 ; k < nq ; k++ ) q[k]   = ( t == 0 ) ? 0 : 3.0 * frand();
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
        fwd_dyn_crba( q, qd, tau, fext, qddGen );
        fd_dense( &mdl, q, qd, tau, fext, qddRef );

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
