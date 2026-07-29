// Trajectory harness: steps the generated dynamics forward in time and checks
// the invariants that only a trajectory can check.
//
// src/validate proves the generated qdd matches the dense reference at a
// state. That says nothing about whether qd is really dq/dt, whether the
// quaternion update is right, or whether a force model is passive -- all of
// which are properties of a trajectory. This harness covers that gap.
//
// Stage 0 (this file, today) runs the integrator-only check: with a
// symplectic scheme the energy error is O(dt), so halving dt must halve the
// drift. That is the same discipline the Chrono cross-check uses to separate
// truncation error from a real discrepancy -- a bug in the integrator or in
// the qd/dq-dt convention shows up as a drift that does NOT scale.
//
// Usage: sim_<model>.exe <model>.mech

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "integrate.h"
#include "../contact/track_contact.h"

// The generated implementation under test.
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *act, double *out_qdd );
extern const int fwd_dyn_crba_nbod;
extern const int fwd_dyn_crba_nq;
extern const int fwd_dyn_crba_nv;
extern const int fwd_dyn_crba_nact;

// Seed q from the model's design pose. Only the track bands have a nonzero
// one; everything else is zeros with an identity quaternion per floating body
// (see modelLoad). This is the first consumer in the repo to actually start a
// band at its design pose -- the Chrono harness reads q0 to PLACE its bodies
// but then simulates from a scrambled state instead.
static void seedDesignPose( const Model *mdl, double *q, double *qd ) {
    memcpy( q, mdl->q0, mdl->nq * sizeof *q );
    memset( qd, 0, mdl->nv * sizeof *qd );
}

// Runs nsteps of semi-implicit Euler from the design pose with no applied
// load, and returns the largest |E(t) - E(0)| seen along the way.
//
// tau, fext and act are all zero, so the only forces are gravity and whatever
// the model's own force elements produce. Damping is what would break the
// check, so models carrying dampers are excluded by the caller.
static double energyDrift( const Model *mdl, double dt, int nsteps ) {
    static double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ], qdd[ MODEL_MAX_V ];
    static double tau[ MODEL_MAX_V ], fext[ MODEL_MAX_BODIES * 6 ], act[ MODEL_MAX_ELEM ];

    seedDesignPose( mdl, q, qd );
    memset( tau,  0, sizeof tau );
    memset( fext, 0, sizeof fext );
    memset( act,  0, sizeof act );

    // A mechanism released from rest at a symmetric pose can just sit there,
    // which would pass any energy check trivially. Give it a nudge so the
    // trajectory is actually exercised.
    for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] = 0.3;

    double e0 = simEnergy( mdl, q, qd ), worst = 0;
    for ( int s = 0 ; s < nsteps ; s++ ) {
        fwd_dyn_crba( q, qd, tau, fext, act, qdd );
        simStep( mdl, q, qd, qdd, dt );
        double e = fabs( simEnergy( mdl, q, qd ) - e0 );
        if ( e > worst ) worst = e;
    }
    return worst;
}

// ---------------------------------------------------------------------------
// Track models: settle the band onto its rollers, then check the trajectory
// invariants that only contact-in-the-loop can exercise.
//
// WHY GRAVITY DOES NOT NEED TURNING OFF. Ground contact is not implemented, so
// nothing holds the vehicle up and the whole assembly is in free fall. That is
// harmless here and it is not a fudge: gravity is uniform, so it contributes
// the same acceleration to every body and cancels out of every RELATIVE
// position and velocity the contact law depends on. The band relaxes onto its
// wheels exactly as it would in orbit; only the absolute altitude runs away.
// It also leaves the momentum check exact rather than approximate, since
// dP/dt = M*g holds identically when all other forces are internal.
//
// WHY THE GAINS RAMP. At the design pose the band wraps only the two end
// pulleys, so the five road wheels per side poke through the lower run by
// 0.17-0.25 m. At full stiffness that is order 1e6 N per side against a
// 1.1e5 N vehicle weight -- the band would be fired off the hull. Ramping the
// stiffness in (and damping it heavily at first) lets the band relax through
// the same force law that will run afterwards, with no separate "assembly"
// code path to get wrong.
// ---------------------------------------------------------------------------

// Worst penetration of any pin into any roller, and how many pairs are loaded.
static double bandWorstPenetration( const Model *mdl, const TrackContactParams *prm,
                                    const double *q, int *nLoaded ) {
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];
    static double qd0[ MODEL_MAX_V ];
    memset( qd0, 0, mdl->nv * sizeof *qd0 );
    denseKinematics( mdl, q, qd0, Xa, vb );

    double worst = 0;
    *nLoaded = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[t];
        for ( int r = 0 ; r < tr->nroller ; r++ ) {
            const ModelTrackRoller *ro = &tr->roller[r];
            int s = ro->body;
            double Es[3][3], c[3];
            for ( int a = 0 ; a < 3 ; a++ ) for ( int e = 0 ; e < 3 ; e++ ) Es[a][e] = Xa[s][a][e];
            xformOrigin( Xa[s], c );
            for ( int i = 0 ; i < tr->shoes ; i++ ) {
                int j = tr->shoe0 + i;
                double p[3], rel[3], d[3];
                xformOrigin( Xa[j], p );
                for ( int k = 0 ; k < 3 ; k++ ) rel[k] = p[k] - c[k];
                for ( int k = 0 ; k < 3 ; k++ )
                    d[k] = Es[k][0]*rel[0] + Es[k][1]*rel[1] + Es[k][2]*rel[2];
                double delta = ro->radius - sqrt( d[0]*d[0] + d[2]*d[2] );
                if ( delta > worst ) worst = delta;
                if ( delta > 0 ) (*nLoaded)++;
            }
        }
    }
    (void)prm;
    return worst;
}

static int runTrackSim( const Model *mdl, const char *path ) {
    static double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ], qdd[ MODEL_MAX_V ];
    static double tau[ MODEL_MAX_V ], fext[ MODEL_MAX_BODIES * 6 ], act[ MODEL_MAX_ELEM ];
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];

    char err[ 256 ];
    if ( trackContactCheckModel( mdl, err, sizeof err ) ) {
        printf( "sim: FAIL: %s\n", err );
        return 1;
    }

    TrackContactParams prm;
    trackContactDefaults( mdl, &prm );

    memcpy( q, mdl->q0, mdl->nq * sizeof *q );
    memset( qd,  0, mdl->nv * sizeof *qd );
    memset( tau, 0, mdl->nv * sizeof *tau );
    memset( act, 0, mdl->nact * sizeof *act );

    int nLoaded0;
    double pen0 = bandWorstPenetration( mdl, &prm, q, &nLoaded0 );

    const double DT = 1e-4, T_RAMP = 0.2, T_SETTLE = 0.6;
    int nsteps = (int)( T_SETTLE / DT );
    double mTotal = 0;
    for ( int i = 0 ; i < mdl->nbod ; i++ ) mTotal += mdl->body[i].mass;

    double ePrev = 0, eFirst = 0, eLast = 0, keScale = 0;
    int    energyRises = 0;
    double worstRise = 0;

    // Momentum is only conserved once the artificial viscosity is gone: that
    // term multiplies the whole of qd, rigid-body modes included, so while it
    // is active it is an external force and dP/dt = M*g does not hold. Latch a
    // reference state at the moment it switches off and measure from there.
    double Pref[3] = { 0, 0, 0 }, tRef = 0;
    int    haveRef = 0;

    for ( int s = 0 ; s < nsteps ; s++ ) {
        double t = s * DT;
        double x = t < T_RAMP ? t / T_RAMP : 1.0;
        double sm = x*x*( 3.0 - 2.0*x );
        prm.gain = sm;
        prm.zeta = 1.0 + ( 0.2 - 1.0 ) * sm;      // over-damped early, 0.2 running

        denseKinematics( mdl, q, qd, Xa, vb );
        memset( fext, 0, mdl->nbod * 6 * sizeof *fext );
        trackContactFext( mdl, &prm, Xa, vb, fext );
        fwd_dyn_crba( q, qd, tau, fext, act, qdd );
        simStep( mdl, q, qd, qdd, DT );

        // Artificial viscosity, faded out with the ramp: it drains the
        // relaxation transient without biasing the settled state.
        double beta = 20.0 * ( 1.0 - sm );
        if ( beta > 0 ) { double f = 1.0 - beta*DT; for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] *= f; }

        // Once the gains are at full strength and the artificial viscosity is
        // gone, every remaining force is either conservative or dissipative,
        // so the total energy must not increase. This is the trajectory-level
        // statement of the passivity the contact tests check pointwise, and it
        // is what a sign error or an energy-injecting damper would break.
        if ( t > T_RAMP * 1.5 ) {
            // Full budget: kinetic + gravity + the model's own springs +
            // the contact springs. The contact term is only well defined at
            // constant gain, which is exactly why this window starts after
            // the ramp -- while the gain moves, the potential is explicitly
            // time-dependent and energy is not expected to be conserved.
            static double Xe[ MODEL_MAX_BODIES ][6][6], ve[ MODEL_MAX_BODIES ][6];
            denseKinematics( mdl, q, qd, Xe, ve );
            double ke;
            double e = simEnergyK( mdl, q, qd, &ke ) + trackContactPotential( mdl, &prm, Xe );
            if ( ke > keScale ) keScale = ke;
            if ( ePrev != 0 && e > ePrev ) {
                double rise = e - ePrev;
                if ( rise > worstRise ) worstRise = rise;
                energyRises++;
            }
            if ( !haveRef ) eFirst = e;
            eLast = e;
            ePrev = e;

            if ( !haveRef ) {
                simMomentum( mdl, q, qd, Pref );
                tRef = t + DT;          // state is post-step
                haveRef = 1;
            }
        }
    }

    int nLoaded1;
    double pen1 = bandWorstPenetration( mdl, &prm, q, &nLoaded1 );

    printf( "sim: %s: settled %d steps (%.2f s, dt=%.0e): worst pin penetration %.4f -> %.4f m, "
            "loaded pairs %d -> %d\n",
            path, nsteps, T_SETTLE, DT, pen0, pen1, nLoaded0, nLoaded1 );

    // Momentum: with every contact force internal and gravity uniform,
    // dP/dt = M*g exactly, so over the post-ramp window dP must be M*g*dt.
    // This is the check that a dropped reaction cannot survive -- and note it
    // is blind to the contact model's details, which is the point.
    double P[3], worstMom = 0;
    double span = nsteps * DT - tRef;
    simMomentum( mdl, q, qd, P );
    // Normalize against the momentum actually present, not against the
    // component's own target: gravity is z-only, so dividing by |M*g_k*span|
    // would make the x and y components an absolute test with a denominator
    // of 1 kg.m/s -- which on an 11 t vehicle is no test at all.
    double Pmag = 0;
    for ( int k = 0 ; k < 3 ; k++ ) Pmag += fabs( P[k] );
    if ( Pmag < 1.0 ) Pmag = 1.0;
    for ( int k = 0 ; k < 3 ; k++ ) {
        double want = mTotal * mdl->gravity[k] * span;
        double got  = P[k] - Pref[k];
        double e = fabs( got - want ) / Pmag;
        if ( e > worstMom ) worstMom = e;
    }
    printf( "sim: %s: linear momentum over the %.2f s post-ramp window matches M*g*dt "
            "to %.3e relative\n", path, span, worstMom );

    if ( pen1 > pen0 ) {
        printf( "sim: FAIL: settling made the interpenetration worse (%.4f -> %.4f m)\n", pen0, pen1 );
        return 1;
    }
    // Energy: everything left after the ramp is either conservative or
    // dissipative, so the budget must shrink. Two statements, because a
    // strictly monotone assertion would be wrong: semi-implicit Euler's energy
    // error OSCILLATES within an O(dt) band rather than decreasing, so a
    // first-order integrator on a stiff penalty law ripples upward now and
    // then no matter how correct the forces are.
    //   - net change over the window must be a loss;
    //   - no single-step rise may reach a truncation-scale fraction of the
    //     kinetic energy. A genuine energy-injecting force (a sign error, or
    //     the gated potential without its radial gradient, which is how this
    //     check earned its keep) shows up orders of magnitude above this.
    double eNet   = eLast - eFirst;
    double eBound = 1e-4 * fmax( keScale, 1.0 );
    printf( "sim: %s: energy over the post-ramp window: net %+.4g J, worst single-step rise "
            "%.4g J against a %.4g J bound (%d rises, KE scale %.4g J)\n",
            path, eNet, worstRise, eBound, energyRises, keScale );
    if ( eNet > eBound ) {
        printf( "sim: FAIL: total energy grew by %.4g J after the ramp -- "
                "the contact law is injecting energy\n", eNet );
        return 1;
    }
    if ( worstRise > eBound ) {
        printf( "sim: FAIL: a single step added %.4g J (bound %.4g J) -- "
                "that is too large to be truncation\n", worstRise, eBound );
        return 1;
    }
    // Contact is a stiff penalty law and the integrator is first order, so the
    // momentum residual is truncation, not a modelling error; a dropped
    // reaction would show up as a residual of order one, not of order dt.
    if ( worstMom > 1e-2 ) {
        printf( "sim: FAIL: linear momentum does not follow M*g*t (relative error %.3e) -- "
                "a contact reaction is being dropped\n", worstMom );
        return 1;
    }
    printf( "sim: PASS\n" );
    return 0;
}

int main( int argc, char **argv ) {
    const char *path = argc > 1 ? argv[ 1 ] : "../models/fourbar.mech";
    Model mdl;
    char err[ 256 ];
    if ( modelLoad( path, &mdl, err, sizeof err ) ) {
        fprintf( stderr, "sim: %s\n", err );
        return 1;
    }
    const char *unsup = modelUnsupported( &mdl );
    if ( unsup ) {
        fprintf( stderr, "sim: %s\n", unsup );
        return 1;
    }
    if ( fwd_dyn_crba_nbod != mdl.nbod || fwd_dyn_crba_nq != mdl.nq
      || fwd_dyn_crba_nv != mdl.nv || fwd_dyn_crba_nact != mdl.nact ) {
        fprintf( stderr, "sim: dynamics was generated for nbod=%d nq=%d nv=%d nact=%d but '%s' has\n"
                         "nbod=%d nq=%d nv=%d nact=%d -- regenerate\n",
                 fwd_dyn_crba_nbod, fwd_dyn_crba_nq, fwd_dyn_crba_nv, fwd_dyn_crba_nact,
                 path, mdl.nbod, mdl.nq, mdl.nv, mdl.nact );
        return 1;
    }

    // Track models get the contact-in-the-loop run instead of the bare
    // integrator check: they have both damping and closure loops, so the
    // energy-scaling test below does not apply to them, and settling the band
    // is the thing actually worth asserting.
    if ( mdl.ntrack ) return runTrackSim( &mdl, path );

    // Two classes of model cannot gate an energy check, for reasons that are
    // properties of the model rather than of the integrator:
    //
    //  - Damping dissipates, so energy is not conserved at all.
    //  - Loop models close their constraints with Baumgarte stabilization,
    //    which is a feedback term on constraint error: it injects and removes
    //    energy by construction and is not symplectic. Its drift is real and
    //    is NOT first order in dt, so gating on it would be gating on the
    //    stabilization scheme -- the same trap the validator avoids by
    //    sampling near the closure manifold. Loop models are exercised by the
    //    other checks instead.
    //
    // Enough tree models remain to make this a real gate, including satellite,
    // whose floating bodies are the only thing that exercises integrateQ's
    // quaternion path -- the one piece of this file that is not q += qd*dt.
    for ( int e = 0 ; e < mdl.nelem ; e++ ) {
        if ( mdl.elem[e].c != 0.0 ) {
            printf( "sim: %s: model has damping (element '%s', c = %g) -- energy check skipped\n",
                    path, mdl.elem[e].name, mdl.elem[e].c );
            return 0;
        }
    }
    if ( mdl.nloop ) {
        printf( "sim: %s: %d closure loop(s) -- Baumgarte stabilization is not symplectic,\n"
                "sim: energy check skipped (see the note in src/sim/main.c)\n", path, mdl.nloop );
        return 0;
    }

    // Same simulated duration at both step sizes, so the comparison is drift
    // per unit time and not drift per step.
    const double T = 0.5, DT = 1e-4;
    double dHi = energyDrift( &mdl, DT,       (int)( T / DT ) );
    double dLo = energyDrift( &mdl, DT / 2.0, (int)( T / ( DT / 2.0 ) ) );
    double ratio = dLo > 0 ? dHi / dLo : 0;

    printf( "sim: %s: energy drift over %.2f s: %.4e at dt=%.1e, %.4e at dt=%.1e, ratio %.3f\n",
            path, T, dHi, DT, dLo, DT / 2.0, ratio );

    // First order means halving dt halves the drift, i.e. ratio -> 2. The band
    // is wide because the drift of a symplectic scheme oscillates rather than
    // accumulating, so the max over a window is not a clean power law; what
    // this rejects is a drift that does not shrink at all (ratio ~ 1, a
    // convention or integrator bug) or one that blows up.
    if ( dHi < 1e-12 ) {
        printf( "sim: %s: drift below the noise floor -- nothing to scale\n", path );
        printf( "sim: PASS\n" );
        return 0;
    }
    if ( ratio > 1.5 && ratio < 4.0 ) {
        printf( "sim: PASS\n" );
        return 0;
    }
    printf( "sim: FAIL: energy drift is not first order in dt (ratio %.3f, expected ~2)\n", ratio );
    return 1;
}
