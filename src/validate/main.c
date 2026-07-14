// Validation harness: runs the generated sparse dynamics (build/dynamics.c)
// and the dense 6x6 reference (dense.h) on random states and compares qdd.
// Both take the mechanism from the same model file. Exit code 0 = PASS,
// 1 = FAIL. See dense.h for the reference implementation.

#include <stdio.h>
#include "dense.h"

// The generated implementation under test (compiled from build/dynamics.c).
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *out_qdd );
extern const int fwd_dyn_crba_nbod;

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
    if ( fwd_dyn_crba_nbod != mdl.nbod ) {
        fprintf( stderr, "validate: dynamics.c was generated for %d bodies but '%s' has %d -- regenerate\n",
                 fwd_dyn_crba_nbod, path, mdl.nbod );
        return 1;
    }
    int nb = mdl.nbod;

    const int NTRIALS = 1000;
    const double TOL = 1e-9;
    double maxErr = 0;
    int worstTrial = -1;

    for ( int t = 0 ; t < NTRIALS ; t++ ) {
        double q[ MODEL_MAX_BODIES ], qd[ MODEL_MAX_BODIES ], tau[ MODEL_MAX_BODIES ];
        double fext[ MODEL_MAX_BODIES * 6 ];

        // First trials pin down special cases; the rest are fully random.
        for ( int i = 0 ; i < nb ; i++ ) {
            q[i]   = ( t == 0 ) ? 0 : 3.0 * frand();
            qd[i]  = ( t <= 1 ) ? 0 : 2.0 * frand();
            tau[i] = ( t <= 2 ) ? 0 : 5.0 * frand();
        }
        for ( int k = 0 ; k < nb * 6 ; k++ )
            fext[k] = ( t <= 3 ) ? 0 : 5.0 * frand();

        double qddGen[ MODEL_MAX_BODIES ], qddRef[ MODEL_MAX_BODIES ];
        fwd_dyn_crba( q, qd, tau, fext, qddGen );
        fd_dense( &mdl, q, qd, tau, fext, qddRef );

        for ( int i = 0 ; i < nb ; i++ ) {
            double e = fabs( qddGen[i] - qddRef[i] );
            if ( e > maxErr ) { maxErr = e; worstTrial = t; }
        }
    }

    printf( "validate: %s: %d trials, max |qdd_generated - qdd_dense| = %.3e (worst at trial %d)\n",
            path, NTRIALS, maxErr, worstTrial );
    if ( maxErr < TOL ) {
        printf( "validate: PASS\n" );
        return 0;
    }
    printf( "validate: FAIL (tolerance %.1e)\n", TOL );
    return 1;
}
