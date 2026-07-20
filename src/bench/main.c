// Benchmark: generated sparse dynamics (build/dynamics.c) vs the dense 6x6
// reference (src/validate/dense.h), same inputs, wall-clock timed with
// QueryPerformanceCounter. Build with /O2 -- timing unoptimized code says
// nothing useful (build.bat does this). Both take the mechanism from the
// same model file.
//
// Each implementation cycles through a pool of pre-generated random states so
// the compiler can't specialize on constant inputs, and a running checksum of
// the outputs is printed so the calls can't be dead-code eliminated. Timing
// runs in batches until a minimum elapsed time, so the per-call figure is
// stable regardless of how fast the machine is.

#include <stdio.h>
#include <windows.h>
#include "../validate/dense.h"

// The generated implementation (compiled from build/dynamics.c).
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *act, double *out_qdd );
extern const int fwd_dyn_crba_nq;
extern const int fwd_dyn_crba_nv;
extern const int fwd_dyn_crba_nact;

#define NSTATES 256 // power of two so the index wrap is a mask

static Model g_mdl;
static double qs[ NSTATES ][ MODEL_MAX_Q ];
static double qds[ NSTATES ][ MODEL_MAX_V ];
static double taus[ NSTATES ][ MODEL_MAX_V ];
static double fexts[ NSTATES ][ MODEL_MAX_BODIES * 6 ];
static double acts[ NSTATES ][ MODEL_MAX_ELEM ];

typedef void (*FdFn)( const double *q, const double *qd, const double *tau,
                      const double *fext, const double *act, double *qdd );

// The generated function's signature is non-const; the inputs are never
// written, so adapting via cast is safe here.
static void genAdapter( const double *q, const double *qd, const double *tau,
                        const double *fext, const double *act, double *qdd ) {
    fwd_dyn_crba( (double*)q, (double*)qd, (double*)tau, (double*)fext, (double*)act, qdd );
}

static void denseAdapter( const double *q, const double *qd, const double *tau,
                          const double *fext, const double *act, double *qdd ) {
    fd_dense( &g_mdl, q, qd, tau, fext, act, qdd );
}

// Runs fn in batches until minSeconds has elapsed; returns ns per call.
static double benchOne( FdFn fn, double minSeconds, double *checksumPerCall ) {
    LARGE_INTEGER freq, t0, t1;
    QueryPerformanceFrequency( &freq );

    double qdd[ MODEL_MAX_V ];
    double checksum = 0;
    long long calls = 0;
    double elapsed = 0;

    QueryPerformanceCounter( &t0 );
    do {
        for ( int b = 0 ; b < 10000 ; b++ ) {
            int s = (int)( calls & ( NSTATES - 1 ) );
            fn( qs[s], qds[s], taus[s], fexts[s], acts[s], qdd );
            checksum += qdd[ 0 ] + qdd[ g_mdl.nv - 1 ];
            calls++;
        }
        QueryPerformanceCounter( &t1 );
        elapsed = (double)( t1.QuadPart - t0.QuadPart ) / (double)freq.QuadPart;
    } while ( elapsed < minSeconds );

    *checksumPerCall = checksum / (double)calls;
    return elapsed * 1e9 / (double)calls;
}

int main( int argc, char **argv ) {
    const char *path = argc > 1 ? argv[ 1 ] : "../models/threelink.mech";
    char err[ 256 ];
    if ( modelLoad( path, &g_mdl, err, sizeof err ) ) {
        fprintf( stderr, "bench: %s\n", err );
        return 1;
    }
    const char *unsup = modelUnsupported( &g_mdl );
    if ( unsup ) {
        fprintf( stderr, "bench: %s\n", unsup );
        return 1;
    }
    if ( fwd_dyn_crba_nq != g_mdl.nq || fwd_dyn_crba_nv != g_mdl.nv || fwd_dyn_crba_nact != g_mdl.nact ) {
        fprintf( stderr, "bench: dynamics.c was generated for nq=%d nv=%d nact=%d but '%s' has nq=%d nv=%d nact=%d -- regenerate\n",
                 fwd_dyn_crba_nq, fwd_dyn_crba_nv, fwd_dyn_crba_nact, path, g_mdl.nq, g_mdl.nv, g_mdl.nact );
        return 1;
    }

    for ( int s = 0 ; s < NSTATES ; s++ ) {
        for ( int k = 0 ; k < g_mdl.nq ; k++ )     qs[s][k]   = 3.0 * frand();
        for ( int k = 0 ; k < g_mdl.nv ; k++ )     qds[s][k]  = 2.0 * frand();
        for ( int k = 0 ; k < g_mdl.nv ; k++ )     taus[s][k] = 5.0 * frand();
        for ( int k = 0 ; k < g_mdl.nbod * 6 ; k++ ) fexts[s][k] = 5.0 * frand();
        for ( int k = 0 ; k < g_mdl.nact ; k++ )     acts[s][k]  = 5.0 * frand();
    }

    // Warm-up (page in code, settle clocks), then measure.
    double cs;
    benchOne( genAdapter,   0.05, &cs );
    benchOne( denseAdapter, 0.05, &cs );

    double csGen, csDense;
    double nsGen   = benchOne( genAdapter,   1.0, &csGen );
    double nsDense = benchOne( denseAdapter, 1.0, &csDense );

    printf( "bench: %s (%d bodies, nq %d, nv %d)\n", path, g_mdl.nbod, g_mdl.nq, g_mdl.nv );
    printf( "bench: generated sparse  %8.1f ns/call  (checksum/call %+.6f)\n", nsGen, csGen );
    printf( "bench: dense reference   %8.1f ns/call  (checksum/call %+.6f)\n", nsDense, csDense );
    printf( "bench: speedup %.1fx\n", nsDense / nsGen );
    return 0;
}
