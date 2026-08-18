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
#include "../contact/ground_contact.h"

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
// Track models: settle the band onto its rollers and the vehicle onto the
// ground, then check the trajectory invariants that only contact-in-the-loop
// can exercise.
//
// WHAT GROUND CONTACT CHANGED HERE. Every track contact is internal, so before
// ground existed the whole assembly was in free fall and dP/dt = M*g held
// identically -- which made the momentum check exact and completely blind to
// the contact model's details. Ground reactions are EXTERNAL, so that identity
// is gone and the check had to be rebuilt around the impulse the ground
// actually applied: dP = integral( M*g + F_ground ) dt. It keeps the property
// that earned it its keep, since a dropped reaction still shows up as a
// residual of order one rather than of order dt.
//
// In exchange, ground contact unlocks a stronger check that an internal-only
// band could not offer at all: STATIC EQUILIBRIUM. Once settled, the total
// ground normal load must equal the vehicle's weight. That is "it stands up",
// asserted numerically, and it is independent of every parameter in the force
// law.
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
        double runZ = tr->padZ + tr->thickness;
        double fade = prm->faceFade * tr->linkLen;
        for ( int r = 0 ; r < tr->nroller ; r++ ) {
            const ModelTrackRoller *ro = &tr->roller[r];
            int s = ro->body;
            double Es[3][3], c[3];
            for ( int a = 0 ; a < 3 ; a++ ) for ( int e = 0 ; e < 3 ; e++ ) Es[a][e] = Xa[s][a][e];
            xformOrigin( Xa[s], c );

            // The metric has to match the primitive. A road wheel rolls on the
            // shoe's inner running face, so its penetration is measured to that
            // face -- and its PINS are not constrained at all, sitting a chord's
            // worth OUTSIDE the pitch circle by construction. Measuring them
            // against the circle would report a large penetration for a wheel
            // resting exactly where it should.
            if ( !ro->isEnd ) {
                for ( int i = 0 ; i < tr->shoes ; i++ ) {
                    int j = tr->shoe0 + i;
                    double Ej[3][3], rw[3], rel[3], cc[3];
                    for ( int a = 0 ; a < 3 ; a++ ) for ( int e = 0 ; e < 3 ; e++ ) Ej[a][e] = Xa[j][a][e];
                    xformOrigin( Xa[j], rw );
                    for ( int k = 0 ; k < 3 ; k++ ) rel[k] = c[k] - rw[k];
                    for ( int k = 0 ; k < 3 ; k++ )
                        cc[k] = Ej[k][0]*rel[0] + Ej[k][1]*rel[1] + Ej[k][2]*rel[2];
                    double dd = cc[2] - runZ;
                    if ( dd <= 0 ) continue;
                    if ( trkFaceW( cc[0], tr->linkLen, fade, NULL ) <= 0 ) continue;
                    double delta = ro->radius - dd;
                    if ( delta > worst ) worst = delta;
                    if ( delta > 0 ) (*nLoaded)++;
                }
                continue;
            }

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

// Worst penetration of any shoe pad into the terrain, and how many pad corners
// are loaded. Mirrors bandWorstPenetration; the loaded count is what the
// realized ground damping ratio is computed from (see groundContactDefaults).
static double groundWorstPenetration( const Model *mdl, const double *q, int *nLoaded ) {
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];
    static double qd0[ MODEL_MAX_V ];

    double worst = 0;
    *nLoaded = 0;
    if ( !mdl->ground.has ) return 0;

    memset( qd0, 0, mdl->nv * sizeof *qd0 );
    denseKinematics( mdl, q, qd0, Xa, vb );

    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[t];
        if ( !tr->hasPad ) continue;
        for ( int i = 0 ; i < tr->shoes ; i++ ) {
            int j = tr->shoe0 + i;
            double r[3];
            xformOrigin( Xa[j], r );
            for ( int c = 0 ; c < 4 ; c++ ) {
                double tl[3], p[3], h, hx, hy;
                gndPadCorner( tr, c, tl );
                for ( int k = 0 ; k < 3 ; k++ )
                    p[k] = r[k] + Xa[j][0][k]*tl[0] + Xa[j][1][k]*tl[1] + Xa[j][2][k]*tl[2];
                if ( !terrainAt( &mdl->ground.terrain, p[0], p[1], &h, &hx, &hy ) ) continue;
                double g = h - p[2];
                if ( g > worst ) worst = g;
                if ( g > 0 ) (*nLoaded)++;
            }
        }
    }
    return worst;
}

// Optional dt override from the command line. The momentum balance below is
// asserted to be truncation rather than a modelling error, and the only honest
// way to say that is to re-run at a smaller step and watch it shrink -- the
// same discipline the loop-free energy check uses. Zero means "derive it".
static double gDtOverride = 0;

static int runTrackSim( const Model *mdl, const char *path ) {
    static double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ], qdd[ MODEL_MAX_V ];
    static double tau[ MODEL_MAX_V ], fext[ MODEL_MAX_BODIES * 6 ], act[ MODEL_MAX_ELEM ];
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];

    char err[ 256 ];
    if ( trackContactCheckModel( mdl, err, sizeof err ) ) {
        printf( "sim: FAIL: %s\n", err );
        return 1;
    }
    if ( groundContactCheckModel( mdl, err, sizeof err ) ) {
        printf( "sim: FAIL: %s\n", err );
        return 1;
    }

    TrackContactParams prm;
    trackContactDefaults( mdl, &prm );
    GroundContactParams gprm;
    groundContactDefaults( mdl, &gprm );

    memcpy( q, mdl->q0, mdl->nq * sizeof *q );
    memset( qd,  0, mdl->nv * sizeof *qd );
    memset( tau, 0, mdl->nv * sizeof *tau );
    memset( act, 0, mdl->nact * sizeof *act );

    int nLoaded0, nGnd0;
    double pen0  = bandWorstPenetration( mdl, &prm, q, &nLoaded0 );
    double gpen0 = groundWorstPenetration( mdl, q, &nGnd0 );

    // The design pose must not start meaningfully inside the ground. Asserted
    // rather than quietly corrected: a model whose ride height and ground plane
    // disagree is an authoring error, and a harness that silently lifted the
    // vehicle would hide it. Note that dropping a model ONTO the ground would
    // be just as wrong -- m113's design pose already is its ride height, and
    // lowering it until the lowest pad touched would bury the road wheels,
    // which carry no contact points, half a metre under the terrain.
    //
    // The threshold is the static deflection the model is going to settle to
    // anyway, W/(n*kn). A model authored exactly at its ride height (as
    // track_demo is) starts at a penetration of zero plus rounding, so a bare
    // "> 0" would be a rounding test rather than a modelling one.
    double sinkStatic = 0;
    if ( mdl->ground.has && mdl->ground.kn > 0 ) {
        double W = 0, gMag = 0;
        for ( int i = 0 ; i < mdl->nbod ; i++ ) W += mdl->body[i].mass;
        for ( int k = 0 ; k < 3 ; k++ ) gMag += mdl->gravity[k] * mdl->gravity[k];
        sinkStatic = W * sqrt( gMag ) / ( groundLoadedEstimate( mdl ) * mdl->ground.kn );
    }
    if ( gpen0 > sinkStatic ) {
        printf( "sim: FAIL: the design pose starts %.4g m inside the ground, past the %.4g m "
                "it would settle to -- move the model or the 'ground' plane\n", gpen0, sinkStatic );
        return 1;
    }

    // Step size. Ground friction, not stiffness, is what limits an explicit
    // integrator here (see groundContactStabilityStep), and a light-shoe model
    // like track_demo lands well under the 1e-4 the band alone was happy with.
    // Derive it rather than hard-coding one that happens to suit m113.
    // 4 s, not the 0.6 s the band alone needed. Standing the vehicle up excites
    // a slow heave/pitch mode on the torsion bars -- about 1 s per cycle on
    // m113 -- and the settled load is not meaningful until that has rung down.
    // T_RAMP IS SHORT ON PURPOSE, and it used to be 0.2 s. The ramp existed to
    // relax a band jammed 0.22 m into its own rollers by the two-pulley belt
    // idealization. Now that the belt is built around the full roller hull the
    // design pose starts within 0.009 m of equilibrium, and a long ramp does
    // the opposite of what it was for: while contact is still at low gain the
    // whole vehicle free-falls, and 0.2 s of that is 1/2*9.81*0.2^2 = 0.196 m
    // -- almost exactly the depth the road wheels were then found buried in the
    // band. The gains still ease in, over a window short enough that the drop
    // is sub-millimetre, and the artificial viscosity still drains what is left.
    const double T_RAMP = 0.02, T_SETTLE = 4.0;
    //
    // BOTH engines get a say, and which one binds depends on the model. The
    // ground's friction sets it for a light-shoe band like track_demo; the
    // band's own face friction sets it for m113, where the partition of unity
    // puts a whole road wheel's load on one shoe. Taking only the ground's
    // estimate chose 6.2e-5 s for m113 against the band's 4.1e-5 limit, and the
    // run answered with 10469 single-step energy rises of ~0.05 J while its NET
    // energy still fell -- a marginally unstable damper, not truncation.
    double DT = 1e-4;
    double dtGnd = groundContactStabilityStep( mdl, &gprm );
    double dtBand = trackContactStabilityStep( mdl, &prm );
    double dtStab = dtGnd;
    if ( dtBand > 0 && ( dtStab == 0 || dtBand < dtStab ) ) dtStab = dtBand;
    if ( dtStab > 0 && 0.5 * dtStab < DT ) DT = 0.5 * dtStab;
    if ( gDtOverride > 0 ) DT = gDtOverride;
    int nsteps = (int)( T_SETTLE / DT );
    double mTotal = 0;
    for ( int i = 0 ; i < mdl->nbod ; i++ ) mTotal += mdl->body[i].mass;
    if ( mdl->ground.has )
        printf( "sim: %s: ground: %d pads, kn %.3g N/m/pt, mu %.3g, mEff %.3g kg, "
                "stability dt %.2e s (ground %.2e, band %.2e) -> using dt %.2e s\n",
                path, mdl->ground.npoint, mdl->ground.kn, mdl->ground.mu, gprm.mEff,
                dtStab, dtGnd, dtBand, DT );

    double ePrev = 0, eFirst = 0, eLast = 0, keScale = 0;
    int    energyRises = 0;
    double worstRise = 0;

    // Momentum is only conserved once the artificial viscosity is gone: that
    // term multiplies the whole of qd, rigid-body modes included, so while it
    // is active it is an external force and dP/dt = M*g does not hold. Latch a
    // reference state at the moment it switches off and measure from there.
    double Pref[3] = { 0, 0, 0 }, tRef = 0;
    int    haveRef = 0;

    // Ground impulse, and its value latched at the same instant as Pref. The
    // ground wrench comes straight out of groundContactFext rather than being
    // recomputed, so the balance is measured against the very numbers the
    // dynamics was handed.
    double gImp[3] = { 0, 0, 0 }, gImpRef[3] = { 0, 0, 0 };
    double gAvg[3] = { 0, 0, 0 };
    int    nAvg = 0, avgFrom = nsteps - nsteps / 10;

    for ( int s = 0 ; s < nsteps ; s++ ) {
        double t = s * DT;
        double x = t < T_RAMP ? t / T_RAMP : 1.0;
        double sm = x*x*( 3.0 - 2.0*x );
        prm.gain = sm;
        prm.zeta = 1.0 + ( 0.2 - 1.0 ) * sm;      // over-damped early, 0.2 running
        gprm.gain = sm;
        gprm.zeta = 1.0 + ( 0.3 - 1.0 ) * sm;

        double fGnd[6];
        denseKinematics( mdl, q, qd, Xa, vb );
        memset( fext, 0, mdl->nbod * 6 * sizeof *fext );
        trackContactFext( mdl, &prm, Xa, vb, fext );
        groundContactFext( mdl, &gprm, Xa, vb, fext, fGnd );
        fwd_dyn_crba( q, qd, tau, fext, act, qdd );
        simStep( mdl, q, qd, qdd, DT );

        // Accumulated before the latch below, so gImp at latch time covers
        // exactly the steps whose forces acted before tRef.
        for ( int k = 0 ; k < 3 ; k++ ) gImp[k] += fGnd[ 3 + k ] * DT;
        if ( s >= avgFrom ) { for ( int k = 0 ; k < 3 ; k++ ) gAvg[k] += fGnd[ 3 + k ]; nAvg++; }

        // Artificial viscosity, faded out with the ramp: it drains the
        // relaxation transient without biasing the settled state.
        double beta = 20.0 * ( 1.0 - sm );
        if ( beta > 0 ) { double f = 1.0 - beta*DT; for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] *= f; }

        // Momentum reference, latched the instant the artificial viscosity is
        // gone -- that term multiplies the whole of qd, rigid-body modes
        // included, so while it is active it is an external force and the
        // balance below does not hold. Kept OUT of the energy window: the
        // momentum check costs nothing per step and wants the longest span it
        // can get, since it is the integrated impulse that makes a dropped
        // reaction visible.
        if ( !haveRef && t > T_RAMP * 1.5 ) {
            simMomentum( mdl, q, qd, Pref );
            for ( int k = 0 ; k < 3 ; k++ ) gImpRef[k] = gImp[k];
            tRef = t + DT;              // state is post-step
            haveRef = 1;
        }

        // Once the gains are at full strength and the artificial viscosity is
        // gone, every remaining force is either conservative or dissipative,
        // so the total energy must not increase. This is the trajectory-level
        // statement of the passivity the contact tests check pointwise, and it
        // is what a sign error or an energy-injecting damper would break.
        //
        // Restricted to the settled tail rather than the whole post-ramp
        // window, and both reasons point the same way: the budget costs a
        // second denseKinematics plus simEnergyK every step, which on 153
        // bodies dominates the run; and energy injection is easiest to see once
        // the relaxation transient is over, so the quiet quarter is the
        // stronger test as well as the cheaper one.
        if ( s >= avgFrom - nsteps / 4 && t > T_RAMP * 1.5 ) {
            // Full budget: kinetic + gravity + the model's own springs + both
            // sets of contact springs. Each term is only well defined at
            // constant gain, which is why this window is well clear of the
            // ramp -- while a gain moves, that potential is explicitly
            // time-dependent and energy is not expected to be conserved.
            static double Xe[ MODEL_MAX_BODIES ][6][6], ve[ MODEL_MAX_BODIES ][6];
            denseKinematics( mdl, q, qd, Xe, ve );
            double ke;
            double e = simEnergyK( mdl, q, qd, &ke ) + trackContactPotential( mdl, &prm, Xe )
                     + groundContactPotential( mdl, &gprm, Xe );
            if ( ke > keScale ) keScale = ke;
            if ( ePrev != 0 && e > ePrev ) {
                double rise = e - ePrev;
                if ( rise > worstRise ) worstRise = rise;
                energyRises++;
            }
            if ( eFirst == 0 ) eFirst = e;
            eLast = e;
            ePrev = e;
        }
    }

    int nLoaded1, nGnd1;
    double pen1  = bandWorstPenetration( mdl, &prm, q, &nLoaded1 );
    double gpen1 = groundWorstPenetration( mdl, q, &nGnd1 );

    printf( "sim: %s: settled %d steps (%.2f s, dt=%.0e): worst pin penetration %.4f -> %.4f m, "
            "loaded pairs %d -> %d\n",
            path, nsteps, T_SETTLE, DT, pen0, pen1, nLoaded0, nLoaded1 );

    // Static equilibrium: the check ground contact makes possible and an
    // internal-only band could not offer. Averaged over the last tenth of the
    // run because a penalty contact settles into an oscillation about the
    // static load rather than onto it.
    int gndOk = 1;
    if ( mdl->ground.has && nAvg > 0 ) {
        for ( int k = 0 ; k < 3 ; k++ ) gAvg[k] /= nAvg;
        // At rest the ground carries the whole weight vector, sum F = -M*g, so
        // comparing components is exact on a slope too and not just on the
        // flat. The lateral components are printed but NOT asserted: they are
        // the right check for slope holding, which is the next stage's job, and
        // on flat ground a band with nothing but a penalty guide horn drifts
        // slowly sideways -- track_demo shows ~10% of weight in y, which is the
        // demo's placeholder geometry rather than the ground law.
        double W = -mTotal * mdl->gravity[2];
        double relW = W != 0 ? fabs( gAvg[2] - W ) / fabs( W ) : fabs( gAvg[2] );
        // zeta_eff = zeta * sqrt( n_loaded / npoint ); see groundContactDefaults.
        double zEff = mdl->ground.npoint > 0
                    ? 0.3 * sqrt( (double)nGnd1 / mdl->ground.npoint ) : 0;
        printf( "sim: %s: ground: penetration %.4f -> %.4f m, loaded pads %d -> %d of %d "
                "(realized zeta %.3f)\n",
                path, gpen0, gpen1, nGnd0, nGnd1, mdl->ground.npoint, zEff );
        printf( "sim: %s: settled ground load %.4g N vs vehicle weight %.4g N (%.2f%% off); "
                "lateral (%.3g, %.3g) N\n",
                path, gAvg[2], W, relW * 100.0, gAvg[0], gAvg[1] );
        // 5%: the vehicle is still ringing down a 1.7 m/s landing at the end of
        // the window, and this is a settling assertion, not a precision one.
        // A vehicle that never stood up, or one held up by a sign error, misses
        // by order 100%, not by 5.
        if ( relW > 0.05 ) {
            printf( "sim: FAIL: the ground is not carrying the vehicle's weight "
                    "(%.4g N against %.4g N) -- it has not settled or is not standing up\n",
                    gAvg[2], W );
            gndOk = 0;
        }
        if ( gpen1 <= 0 ) {
            printf( "sim: FAIL: no pad is touching the ground after settling\n" );
            gndOk = 0;
        }
    }
    if ( !gndOk ) return 1;

    // Momentum. Before ground contact every force was internal and gravity
    // uniform, so dP/dt = M*g held identically. Ground reactions are external,
    // so the target now carries their impulse as well -- measured from the
    // wrench groundContactFext actually applied, not recomputed. What the check
    // still catches is unchanged: a dropped reaction leaves a residual of order
    // one, not of order dt.
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
    // Denominator. Pmag alone was right when every force was internal and the
    // whole vehicle was in free fall, because then P really was the largest
    // thing around. Standing on the ground it is not: the vehicle ends up
    // nearly at rest, so P collapses while the impulses that have to cancel to
    // produce that rest are each ~4e5 N.s. Dividing by Pmag then measures how
    // still the vehicle happens to have ended up, not whether a force went
    // missing -- it moved by a factor of 10 across three step sizes while the
    // absolute residual barely changed. So the scale is the larger of the
    // momentum present and the impulse actually applied. One scalar for all
    // three components, not a per-component target: gravity is z-only, so a
    // per-component denominator would leave x and y divided by ~zero.
    double impScale = 0;
    for ( int k = 0 ; k < 3 ; k++ )
        impScale += fabs( mTotal * mdl->gravity[k] * span ) + fabs( gImp[k] - gImpRef[k] );
    double scale = Pmag > impScale ? Pmag : impScale;
    double worstAbs = 0;
    for ( int k = 0 ; k < 3 ; k++ ) {
        double want = mTotal * mdl->gravity[k] * span + ( gImp[k] - gImpRef[k] );
        double got  = P[k] - Pref[k];
        double a = fabs( got - want );
        if ( a > worstAbs ) worstAbs = a;
        double e = a / scale;
        if ( e > worstMom ) worstMom = e;
    }
    printf( "sim: %s: linear momentum over the %.2f s post-ramp window matches "
            "M*g*dt + ground impulse to %.3e relative (%.4g N.s against a %.4g N.s scale)\n",
            path, span, worstMom, worstAbs, scale );

    // Band interpenetration. "It must go down" was the right assertion when the
    // only thing a settle did was relax a band jammed into its own geometry.
    // With the vehicle standing on the ground it is wrong: the load path runs
    // road wheel -> pin -> pad -> terrain, so a band that starts TANGENT to its
    // rollers (track_demo, by construction) has to compress into them to carry
    // anything at all, and 0 -> 0.0015 m is the mechanism working rather than
    // failing. What is still worth asserting is both halves of the original
    // intent, stated separately:
    //
    //   - the band must never be crushed onto its axles, loaded or not;
    //   - a band that DID start jammed must have relaxed.
    //
    // The crush threshold is a quarter of the smallest roller radius: a pin
    // that far inside the pitch circle is not in contact with a wheel in any
    // meaningful sense. m113 settles to 0.0348 m against 0.0613, having started
    // at 0.2237.
    double rMin = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ )
        for ( int r = 0 ; r < mdl->track[t].nroller ; r++ ) {
            double rr = mdl->track[t].roller[r].radius;
            if ( rMin == 0 || rr < rMin ) rMin = rr;
        }
    double crush = 0.25 * rMin;
    if ( pen1 > crush ) {
        printf( "sim: FAIL: the band is crushed into its rollers (%.4f m against a %.4f m "
                "limit, a quarter of the smallest roller radius)\n", pen1, crush );
        return 1;
    }
    if ( pen0 > crush && pen1 >= pen0 ) {
        printf( "sim: FAIL: the band started jammed and settling did not relax it "
                "(%.4f -> %.4f m)\n", pen0, pen1 );
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
    // momentum residual is truncation, not a modelling error. That is asserted
    // rather than assumed: on m113, halving dt twice from 6.22e-5 takes the
    // absolute residual 33.4 -> 4.91 -> 1.64 N.s against a scale that stays put
    // at 8.25e5, so it is converging at least first order. A dropped reaction
    // does not converge at all -- it is order one, and here that means order
    // 1e5. Hence the tolerance can stay where it was even though the quantity
    // being normalized changed underneath it; m113 passes with ~250x margin.
    // Re-run with an explicit dt argument to repeat the scaling check.
    if ( worstMom > 1e-2 ) {
        printf( "sim: FAIL: linear momentum does not follow M*g*t + ground impulse "
                "(relative error %.3e) -- a contact reaction is being dropped\n", worstMom );
        return 1;
    }
    printf( "sim: PASS\n" );
    return 0;
}

int main( int argc, char **argv ) {
    const char *path = argc > 1 ? argv[ 1 ] : "../models/fourbar.mech";
    if ( argc > 2 ) gDtOverride = atof( argv[ 2 ] );   // track models only
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
