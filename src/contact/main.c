// Invariant and unit tests for the track contact engine. No trajectory and no
// dynamics: this harness only ever calls trackContactFext / trkPairForce, so
// it runs in milliseconds and can be applied to every model in build.bat --
// models without a track band trivially produce no contact forces, which is
// itself worth asserting.
//
// The checks, in the order they were written (which is also their order of
// value):
//
//   1. Newton's third law. Every contact here is INTERNAL, so the sum of all
//      contact wrenches must vanish. Each pair contributes [p x F; F] and its
//      exact negation, computed from the same p and the same F, so the
//      cancellation is exact PAIRWISE -- but not bitwise in the total, because
//      the two halves land in different accumulators (a shoe's slot sums a
//      handful of rollers, a roller's slot sums 63 shoes) and floating-point
//      addition is not associative. What survives is roundoff proportional to
//      the magnitudes involved, so the check is relative. It still catches
//      every sign error, every wrong body index, every dropped reaction:
//      those produce a residual of order the force itself, not of order eps.
//   2. Determinism. The engine must be a pure function of (q, qd). Everything
//      else in the project's validation strategy depends on that property.
//   3. Smoothness. A force discontinuity as pins enter and leave the wrap
//      makes an explicit integrator chatter.
//   4. The tooth law's fixed points, amplitude and periodicity.
//   5. Design-pose regression: the band must not be jammed into the pulleys
//      that define its own geometry.
//
// Usage: contact_<model>.exe <model>.mech

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../contact/track_contact.h"
#include "../contact/ground_contact.h"

static unsigned long rngState = 123456789UL;
static double frand2( void ) {          // uniform in [-1, 1)
    rngState = rngState * 1103515245UL + 12345UL;
    return ( (double)( ( rngState >> 16 ) & 0x7FFF ) / 16384.0 ) - 1.0;
}

static void randomState( const Model *mdl, double *q, double *qd ) {
    for ( int k = 0 ; k < mdl->nq ; k++ ) q[k]  = 3.0 * frand2();
    for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] = 2.0 * frand2();
    // A floating body's quaternion slice must be a pose, not four random
    // numbers that happen to be near zero; jointXS normalizes, but an
    // all-but-zero quaternion would amplify noise for no reason.
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        if ( mdl->body[i].jtype != JOINT_FLOATING ) continue;
        double *s = q + mdl->qidx[i];
        s[0] += 2.0;
    }
}

// ---------------------------------------------------------------------------
// 1 + 2: internal-force balance and determinism, over random states.
// ---------------------------------------------------------------------------
static int checkBalanceAndPurity( const Model *mdl, const TrackContactParams *prm ) {
    static double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ];
    static double fext[ MODEL_MAX_BODIES * 6 ], fext2[ MODEL_MAX_BODIES * 6 ];
    const int NTRIALS = 1000;
    const double TOL = 1e-12;   // relative; see the header note on why not 0
    int touched = 0;
    double worstRel = 0;

    for ( int t = 0 ; t < NTRIALS ; t++ ) {
        randomState( mdl, q, qd );
        trackContactFromState( mdl, prm, q, qd, fext );

        // Residual against the magnitude that produced it: sum |fext| over
        // bodies is the right scale, because a set of wrenches that cancels
        // perfectly still has each half of every pair sitting in the array at
        // full size.
        double sum[6] = { 0, 0, 0, 0, 0, 0 }, scale[6] = { 0, 0, 0, 0, 0, 0 };
        for ( int i = 0 ; i < mdl->nbod ; i++ )
            for ( int k = 0 ; k < 6 ; k++ ) {
                sum[k]   += fext[ i*6 + k ];
                scale[k] += fabs( fext[ i*6 + k ] );
            }

        for ( int i = 0 ; i < mdl->nbod * 6 ; i++ ) if ( fext[i] != 0.0 ) touched = 1;

        for ( int k = 0 ; k < 6 ; k++ ) {
            double rel = scale[k] > 0 ? fabs( sum[k] ) / scale[k] : fabs( sum[k] );
            if ( rel > worstRel ) worstRel = rel;
            if ( rel > TOL ) {
                printf( "contact: FAIL: contact wrenches do not cancel at trial %d: "
                        "sum[%d] = %.6g against a magnitude of %.6g (relative %.3e, tolerance %.1e)\n",
                        t, k, sum[k], scale[k], rel, TOL );
                return 1;
            }
        }

        trackContactFromState( mdl, prm, q, qd, fext2 );
        for ( int i = 0 ; i < mdl->nbod * 6 ; i++ ) {
            if ( fext[i] != fext2[i] ) {
                printf( "contact: FAIL: engine is not a pure function of (q,qd): "
                        "fext[%d] = %.17g then %.17g at trial %d\n", i, fext[i], fext2[i], t );
                return 1;
            }
        }
    }

    printf( "contact: %d random states: internal wrenches cancel to %.3e relative, "
            "output deterministic%s\n",
            NTRIALS, worstRel, touched ? "" : " (no contact pairs in range -- model has no band)" );
    return 0;
}

// ---------------------------------------------------------------------------
// 3: smoothness of the force law, swept directly through trkPairForce.
// ---------------------------------------------------------------------------
static int checkSmoothness( const Model *mdl, const TrackContactParams *prm ) {
    if ( mdl->ntrack == 0 ) return 0;
    const ModelTrack *tr = &mdl->track[0];

    // Prefer the drive sprocket so the tooth term is in the sweep.
    const ModelTrackRoller *ro = &tr->roller[0];
    for ( int r = 0 ; r < tr->nroller ; r++ ) if ( tr->roller[r].isDrive ) ro = &tr->roller[r];
    double R = ro->radius;

    double u0[3] = { 0, 0, 0 };   // static sweep: no damping terms, so the
                                  // only slope discontinuity that can appear
                                  // is a real one in the elastic law.

    // --- radial sweep -----------------------------------------------------
    // From well outside the gate to well inside contact.
    const int NS = 20000;
    double lo = R - 0.06, hi = R + 0.05, h = ( hi - lo ) / NS;
    double prev = 0, prevOk = 0, maxJump = 0;
    for ( int i = 0 ; i <= NS ; i++ ) {
        double rho = lo + i*h;
        double d[3] = { rho, 0, 0 };
        TrkPairForce pf;
        trkPairForce( prm, tr, ro, d, u0, &pf );
        double f = sqrt( pf.f[0]*pf.f[0] + pf.f[1]*pf.f[1] + pf.f[2]*pf.f[2] );
        if ( prevOk ) { double jump = fabs( f - prev ); if ( jump > maxJump ) maxJump = jump; }
        prev = f; prevOk = 1;
    }
    // The steepest legitimate slope is kn (radial) plus fSkip*teeth/R (tooth),
    // so a step of h may legitimately change |f| by that much. Allow 4x, which
    // still rejects any genuine jump: at the gate a discontinuity would be a
    // finite fraction of fSkip = 2e4 N, thousands of times this bound.
    double slopeMax = prm->gain * ( tr->kn + prm->fSkip * ro->teeth / R );
    double bound = 4.0 * slopeMax * h;
    if ( maxJump > bound ) {
        printf( "contact: FAIL: radial sweep of |f| jumps by %.4g over a step of %.2e "
                "(bound %.4g) -- the engagement gate is discontinuous\n", maxJump, h, bound );
        return 1;
    }

    // --- gate C1 ----------------------------------------------------------
    // Strictly OUTSIDE contact (delta < 0) only the gated tooth and lateral
    // terms act, and there the law must be C1 -- that is what the smoothstep
    // is for. Inside contact it is only C0: a linear penalty spring turns on
    // with slope kn at first touch, which is inherent to the law (Hunt-Crossley
    // would smooth it) and is deliberately not asserted here.
    {
        double a = R + 0.0005, b = R + 0.034;     // delta in (-0.034, -0.0005)
        int    NG = 4000;
        double hg = ( b - a ) / NG, maxCurv = 0;
        for ( int i = 1 ; i < NG ; i++ ) {
            double f[3];
            for ( int k = -1 ; k <= 1 ; k++ ) {
                double d[3] = { a + ( i + k )*hg, 0, 0 };
                TrkPairForce pf;
                trkPairForce( prm, tr, ro, d, u0, &pf );
                f[k+1] = sqrt( pf.f[0]*pf.f[0] + pf.f[1]*pf.f[1] + pf.f[2]*pf.f[2] );
            }
            double curv = fabs( f[2] - 2*f[1] + f[0] ) / ( hg*hg );
            if ( curv > maxCurv ) maxCurv = curv;
        }
        // A C1 function has bounded second difference; a kink would send this
        // to infinity as hg shrinks. The smoothstep's own curvature scales as
        // fSkip/gateFade^2, so bound generously against that.
        double curvBound = 50.0 * prm->gain * prm->fSkip / ( prm->gateFade * prm->gateFade );
        if ( maxCurv > curvBound ) {
            printf( "contact: FAIL: gate is not C1: max |f''| = %.4g exceeds %.4g\n",
                    maxCurv, curvBound );
            return 1;
        }
    }

    // --- circumferential sweep, two full tooth pitches through +/-pi -------
    if ( ro->isDrive && ro->teeth > 0 ) {
        const int NC = 20000;
        double span = 4.0 * M_PI / ro->teeth;
        double th0 = M_PI - span/2, hc = span / NC;
        double prevC = 0, maxJumpC = 0;
        for ( int i = 0 ; i <= NC ; i++ ) {
            double th = th0 + i*hc;
            double d[3] = { R*cos(th), 0, R*sin(th) };   // on the circle: delta = 0
            TrkPairForce pf;
            trkPairForce( prm, tr, ro, d, u0, &pf );
            double f = sqrt( pf.f[0]*pf.f[0] + pf.f[1]*pf.f[1] + pf.f[2]*pf.f[2] );
            if ( i ) { double jump = fabs( f - prevC ); if ( jump > maxJumpC ) maxJumpC = jump; }
            prevC = f;
        }
        // d|f|/dth peaks at fSkip*teeth.
        double boundC = 4.0 * prm->gain * prm->fSkip * ro->teeth * hc;
        if ( maxJumpC > boundC ) {
            printf( "contact: FAIL: circumferential sweep jumps by %.4g crossing th = +/-pi "
                    "(bound %.4g) -- the tooth law is not periodic there\n", maxJumpC, boundC );
            return 1;
        }
        printf( "contact: sweeps smooth: radial max step %.3e N, circumferential %.3e N "
                "(crossing th = +/-pi)\n", maxJump, maxJumpC );
    } else {
        printf( "contact: radial sweep smooth: max step %.3e N (no drive sprocket to sweep)\n", maxJump );
    }
    return 0;
}

// Evaluates the tooth term at an angle, on the pitch circle: delta = 0 there,
// so the radial term is off and the tooth force is all that acts in the plane.
static double toothFT( const TrackContactParams *prm, const ModelTrack *tr,
                       const ModelTrackRoller *ro, double th ) {
    double u0[3] = { 0, 0, 0 };
    double d[3] = { ro->radius*cos(th), 0, ro->radius*sin(th) };
    TrkPairForce pf;
    trkPairForce( prm, tr, ro, d, u0, &pf );
    return pf.fT;
}

// ---------------------------------------------------------------------------
// 4: the tooth law's fixed points, slopes, amplitude and periodicity.
// ---------------------------------------------------------------------------
static int checkToothLaw( const Model *mdl, const TrackContactParams *prm ) {
    if ( mdl->ntrack == 0 ) return 0;
    const ModelTrack *tr = &mdl->track[0];
    const ModelTrackRoller *ro = NULL;
    for ( int r = 0 ; r < tr->nroller ; r++ ) if ( tr->roller[r].isDrive ) ro = &tr->roller[r];
    if ( !ro || ro->teeth <= 0 ) return 0;

    double R = ro->radius, u0[3] = { 0, 0, 0 };
    int    N = ro->teeth;
    double pitch = 2.0 * M_PI / N;

    double eps = 1e-6;
    struct { const char *what; double th; int wantSign; } pts[] = {
        { "pocket centre", prm->toothPhase,               -1 },   // restoring: stable
        { "tooth tip",     prm->toothPhase + pitch/2,     +1 },   // repelling: unstable
    };

    for ( int i = 0 ; i < 2 ; i++ ) {
        double f0 = toothFT( prm, tr, ro, pts[i].th );
        if ( fabs( f0 ) > 1e-9 * prm->fSkip ) {
            printf( "contact: FAIL: tooth force at the %s is %.4g, expected 0\n", pts[i].what, f0 );
            return 1;
        }
        double slope = ( toothFT( prm, tr, ro, pts[i].th + eps ) - toothFT( prm, tr, ro, pts[i].th - eps ) ) / ( 2*eps );
        if ( slope * pts[i].wantSign <= 0 ) {
            printf( "contact: FAIL: tooth force slope at the %s is %.4g, expected sign %+d "
                    "(a pocket must restore, a tip must repel)\n", pts[i].what, slope, pts[i].wantSign );
            return 1;
        }
    }

    // Quarter pitch is the peak: |fT| = fSkip exactly.
    double peak = fabs( toothFT( prm, tr, ro, prm->toothPhase + pitch/4 ) );
    if ( fabs( peak - prm->gain * prm->fSkip ) > 1e-9 * prm->fSkip ) {
        printf( "contact: FAIL: peak tooth force is %.6g, expected fSkip = %.6g\n",
                peak, prm->gain * prm->fSkip );
        return 1;
    }

    // Periodic at the pocket pitch, to machine precision, at an arbitrary
    // phase (not a fixed point, where equality would be trivial).
    double thA = prm->toothPhase + 0.3*pitch;
    double dPer = fabs( toothFT( prm, tr, ro, thA ) - toothFT( prm, tr, ro, thA + pitch ) );
    if ( dPer > 1e-9 * prm->fSkip ) {
        printf( "contact: FAIL: tooth force is not periodic at 2*pi/%d: differs by %.4g\n", N, dPer );
        return 1;
    }

    // Continuous across the atan2 branch cut. This is the property the
    // sinusoid buys for free and an explicit nearest-pocket law would have to
    // get right by hand.
    double dCut = fabs( toothFT( prm, tr, ro, M_PI - 1e-9 ) - toothFT( prm, tr, ro, -M_PI + 1e-9 ) );
    if ( dCut > 1e-6 * prm->fSkip ) {
        printf( "contact: FAIL: tooth force jumps by %.4g across th = +/-pi\n", dCut );
        return 1;
    }
    
    printf( "contact: tooth law: %d pockets at %.6f rad, zero and restoring at the pocket, "
            "zero and repelling at the tip, peak %.0f N, periodic and branch-cut clean\n",
            N, pitch, peak );
    return 0;
}

// ---------------------------------------------------------------------------
// 5: design-pose regression on the path-defining pulleys.
//
// The belt path is built tangent to circles of the declared radii, so the
// design pose ought to be the zero-force pose for the two 'end' rollers. It
// very nearly is -- but NOT because delta == 0. The loader takes each shoe's
// direction from the chord between pins while taking its length from the ARC
// spacing ds, and on a wrap arc the true chord 2R*sin(ds/2R) is shorter than
// ds, so every link overshoots and the reconstructed pins sit a few mm OUTSIDE
// the pitch circle. The radial force is zero because the one-sided clamp
// fires. Assert exactly that, and separately that the pins have not drifted
// far -- a regression that threw the band off the pulleys entirely would
// otherwise still "pass".
// ---------------------------------------------------------------------------
static int checkDesignPose( const Model *mdl, const TrackContactParams *prm ) {
    if ( mdl->ntrack == 0 ) return 0;
    static double Xa[ MODEL_MAX_BODIES ][6][6];
    static double vb[ MODEL_MAX_BODIES ][6];
    static double qd0[ MODEL_MAX_V ];
    memset( qd0, 0, sizeof qd0 );
    denseKinematics( mdl, mdl->q0, qd0, Xa, vb );

    double worstDelta = -1e30, worstClear = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[t];
        for ( int r = 0 ; r < tr->nroller ; r++ ) {
            const ModelTrackRoller *ro = &tr->roller[r];
            if ( !ro->isEnd ) continue;             // road wheels are known to interpenetrate
            int s = ro->body;
            double minClear = 1e30;                 // closest approach of any pin to THIS pulley
            double Es[3][3], c[3];
            for ( int a = 0 ; a < 3 ; a++ ) for ( int b = 0 ; b < 3 ; b++ ) Es[a][b] = Xa[s][a][b];
            xformOrigin( Xa[s], c );

            for ( int i = 0 ; i < tr->shoes ; i++ ) {
                int j = tr->shoe0 + i;
                double p[3], rel[3], d[3];
                xformOrigin( Xa[j], p );
                for ( int k = 0 ; k < 3 ; k++ ) rel[k] = p[k] - c[k];
                for ( int k = 0 ; k < 3 ; k++ )
                    d[k] = Es[k][0]*rel[0] + Es[k][1]*rel[1] + Es[k][2]*rel[2];

                double rho = sqrt( d[0]*d[0] + d[2]*d[2] );
                double delta = ro->radius - rho;
                if ( delta > worstDelta ) worstDelta = delta;
                if ( -delta < minClear ) minClear = -delta;

                double u0[3] = { 0, 0, 0 };
                TrkPairForce pf;
                trkPairForce( prm, tr, ro, d, u0, &pf );
                if ( pf.fN != 0.0 ) {
                    printf( "contact: FAIL: design pose has a nonzero radial force: track '%s' "
                            "roller body %d vs shoe %d, delta = %.6g, fN = %.6g\n",
                            tr->name, s, j, delta, pf.fN );
                    return 1;
                }
            }

            // The belt was built tangent to THIS circle, so some pin must
            // still be hugging it. Measuring the CLOSEST approach is the
            // meaningful test; the largest standoff is not, because a pin one
            // link past a tangent point is legitimately s^2/(2R) away -- 0.046
            // m on the m113 sprocket -- purely because a straight run leaves a
            // circle quadratically.
            //
            // The bound has to be geometric rather than a fixed number of
            // millimetres, because how far a chain of chords can sit off a
            // circle depends entirely on how long a chord is relative to the
            // radius. One chord's sagitta, linkLen^2/(8R), is the natural
            // scale: 11.5 mm on the m113 sprocket (0.150 m links on a 0.245 m
            // radius) but 108 mm on track_demo, whose 0.464 m links on a 0.25
            // m radius are placeholder geometry the model file admits to.
            double sagitta = tr->linkLen * tr->linkLen / ( 8.0 * ro->radius );
            if ( minClear > worstClear ) worstClear = minClear;
            if ( minClear > sagitta ) {
                printf( "contact: FAIL: no design-pose pin comes within one chord sagitta (%.4g m) "
                        "of roller body %d -- closest is %.4g m, so the belt geometry has regressed\n",
                        sagitta, s, minClear );
                return 1;
            }
        }
    }

    if ( worstDelta > 0 ) {
        printf( "contact: FAIL: design pose penetrates a path-defining pulley by %.6g m\n", worstDelta );
        return 1;
    }
    printf( "contact: design pose: no penetration of any path-defining pulley (worst delta %.3g m), "
            "every band still hugs its pulleys (worst closest approach %.3g m), all radial forces zero\n",
            worstDelta, worstClear );
    return 0;
}

// ---------------------------------------------------------------------------
// 6: the drive-torque path, i.e. that the [p x F; F] convention really does
// deliver the tooth force as torque about the sprocket axle.
//
// Two independent routes to the same number. The assembled route reads the
// sprocket's fext 6-vector, whose moment is about the BASE ORIGIN, and shifts
// it to the axle: M_c = M_0 - c x F, then projects on the axle direction. The
// direct route never leaves roller coordinates and just sums (d x f)_y over
// the pairs. If the moment reference point were wrong -- the single easiest
// mistake to make with this convention -- the two would disagree while the
// wrench-balance check in 1 kept passing, because a wrong reference point
// still cancels between the shoe and the roller.
// ---------------------------------------------------------------------------
static int checkDriveTorque( const Model *mdl, const TrackContactParams *prm ) {
    if ( mdl->ntrack == 0 ) return 0;
    static double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ];
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];
    static double fext[ MODEL_MAX_BODIES * 6 ];

    double worstRel = 0, reported = 0;
    for ( int t = 0 ; t < 200 ; t++ ) {
        randomState( mdl, q, qd );
        denseKinematics( mdl, q, qd, Xa, vb );
        memset( fext, 0, mdl->nbod * 6 * sizeof *fext );
        trackContactFext( mdl, prm, Xa, vb, fext );

        for ( int b = 0 ; b < mdl->ntrack ; b++ ) {
            const ModelTrack *tr = &mdl->track[b];
            for ( int r = 0 ; r < tr->nroller ; r++ ) {
                const ModelTrackRoller *ro = &tr->roller[r];
                if ( !ro->isDrive ) continue;
                int s = ro->body;

                double Es[3][3], c[3];
                for ( int a = 0 ; a < 3 ; a++ ) for ( int e = 0 ; e < 3 ; e++ ) Es[a][e] = Xa[s][a][e];
                xformOrigin( Xa[s], c );

                // Assembled: shift the base-origin moment to the axle, then
                // project onto the axle direction (roller local y in world).
                const double *fx = fext + s*6;
                double Mc[3];
                Mc[0] = fx[0] - ( c[1]*fx[5] - c[2]*fx[4] );
                Mc[1] = fx[1] - ( c[2]*fx[3] - c[0]*fx[5] );
                Mc[2] = fx[2] - ( c[0]*fx[4] - c[1]*fx[3] );
                double tqAssembled = 0;
                for ( int k = 0 ; k < 3 ; k++ ) tqAssembled += Mc[k] * Es[1][k];

                // Direct: sum the per-pair moment in roller coordinates. The
                // force on the ROLLER is -f, so its moment about the axle is
                // -(d x f)_y = -(d[2]*f[0] - d[0]*f[2]).
                double ws[3], vs[3];
                for ( int k = 0 ; k < 3 ; k++ ) {
                    ws[k] = Es[0][k]*vb[s][0] + Es[1][k]*vb[s][1] + Es[2][k]*vb[s][2];
                    vs[k] = Es[0][k]*vb[s][3] + Es[1][k]*vb[s][4] + Es[2][k]*vb[s][5];
                }
                double tqDirect = 0, scale = 0;
                for ( int i = 0 ; i < tr->shoes ; i++ ) {
                    int j = tr->shoe0 + i;
                    double p[3], rel[3], d[3], uw[3], u[3];
                    xformOrigin( Xa[j], p );
                    for ( int k = 0 ; k < 3 ; k++ ) rel[k] = p[k] - c[k];
                    for ( int k = 0 ; k < 3 ; k++ )
                        d[k] = Es[k][0]*rel[0] + Es[k][1]*rel[1] + Es[k][2]*rel[2];
                    for ( int k = 0 ; k < 3 ; k++ )
                        uw[k] = Xa[j][0][k]*vb[j][3] + Xa[j][1][k]*vb[j][4] + Xa[j][2][k]*vb[j][5];
                    uw[0] -= vs[0] + ws[1]*rel[2] - ws[2]*rel[1];
                    uw[1] -= vs[1] + ws[2]*rel[0] - ws[0]*rel[2];
                    uw[2] -= vs[2] + ws[0]*rel[1] - ws[1]*rel[0];
                    for ( int k = 0 ; k < 3 ; k++ )
                        u[k] = Es[k][0]*uw[0] + Es[k][1]*uw[1] + Es[k][2]*uw[2];

                    TrkPairForce pf;
                    if ( !trkPairForce( prm, tr, ro, d, u, &pf ) ) continue;
                    double mY = d[2]*pf.f[0] - d[0]*pf.f[2];
                    tqDirect -= mY;
                    scale    += fabs( mY );
                }

                double rel = scale > 0 ? fabs( tqAssembled - tqDirect ) / scale : 0;
                if ( rel > worstRel ) worstRel = rel;
                if ( scale > reported ) reported = scale;
                if ( rel > 1e-9 ) {
                    printf( "contact: FAIL: sprocket axle torque disagrees: %.6g assembled from fext "
                            "vs %.6g summed directly (relative %.3e) -- the moment reference point "
                            "in the fext convention is wrong\n", tqAssembled, tqDirect, rel );
                    return 1;
                }
            }
        }
    }
    printf( "contact: drive torque: axle moment from the assembled fext matches the direct sum "
            "to %.3e relative (torques up to %.3g N.m)\n", worstRel, reported );
    return 0;
}

// ---------------------------------------------------------------------------
// 7: passivity. The velocity-dependent part of the law must never do work on
// the system.
//
// A contact pair delivers power F.u into the bodies, with F the force on the
// shoe and u the shoe's velocity relative to the roller at the contact point.
// Splitting the law by evaluating it at the true u and again at u = 0 isolates
// everything velocity-dependent -- the radial and tooth dampers, the lateral
// damper, and friction, which vanishes at u = 0 by construction -- and that
// part must dissipate at every state, with no tolerance.
//
// This is where the one-sided clamp on f_n earns its keep: without it, the
// radial damper pulls a separating shoe back and injects energy, which shows
// up here as a positive power.
// ---------------------------------------------------------------------------
static int checkPassivity( const Model *mdl, const TrackContactParams *prm ) {
    if ( mdl->ntrack == 0 ) return 0;
    const ModelTrack *tr = &mdl->track[0];
    double worst = 0;

    for ( int t = 0 ; t < 200000 ; t++ ) {
        const ModelTrackRoller *ro = &tr->roller[ ( t / 977 ) % tr->nroller ];
        // Sweep the pin over a band straddling the surface, and the relative
        // velocity over a range that reaches well past vEps in every
        // direction, so separation, approach, spin and lateral slip all occur.
        double d[3], u[3];
        d[0] = ro->radius * ( 1.0 + 0.15*frand2() );
        d[1] = 0.05 * frand2();
        d[2] = ro->radius * ( 0.15*frand2() );
        for ( int k = 0 ; k < 3 ; k++ ) u[k] = 5.0 * frand2();

        TrkPairForce pf, p0;
        double uz[3] = { 0, 0, 0 };
        trkPairForce( prm, tr, ro, d, u,  &pf );
        trkPairForce( prm, tr, ro, d, uz, &p0 );

        double power = 0, mag = 0;
        for ( int k = 0 ; k < 3 ; k++ ) {
            power += ( pf.f[k] - p0.f[k] ) * u[k];
            mag   += fabs( pf.f[k] - p0.f[k] ) * fabs( u[k] );
        }
        if ( power > worst ) worst = power;
        if ( power > 1e-12 * ( mag > 1 ? mag : 1 ) ) {
            printf( "contact: FAIL: velocity-dependent force does %.6g W of work at "
                    "d = (%.4f, %.4f, %.4f), u = (%.3f, %.3f, %.3f) -- the law is not passive\n",
                    power, d[0], d[1], d[2], u[0], u[1], u[2] );
            return 1;
        }
    }
    printf( "contact: passivity: velocity-dependent force dissipates at every sampled state "
            "(worst power %.3e W)\n", worst );
    return 0;
}

// ===========================================================================
// 8, 9, 10: the ROAD-WHEEL FACE law.
//
// Every check above drives trkPairForce, which since the face primitive landed
// is the SPROCKET AND IDLER path only. A road wheel never reaches it, so all
// of them stayed clean while the face law was free to do as it liked --
// checkPassivity in particular printed "worst power 0.000e+00 W" for a model
// whose face friction was injecting 5e6 J over a four-second settle, because
// the states it sampled went through a function the wheels no longer use.
// These three drive trkFaceForce directly, the same way.
// ===========================================================================

// Any non-'end' roller will do; returns 0 if the model has no road wheels.
static const ModelTrack *faceFixture( const Model *mdl, const ModelTrackRoller **roOut ) {
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];
        for ( int r = 0 ; r < tr->nroller ; r++ )
            if ( !tr->roller[ r ].isEnd ) { *roOut = &tr->roller[ r ]; return tr; }
    }
    return 0;
}

// An axle position straddling the face in all three directions, and velocities
// reaching well past vEps, so approach, separation, spin-up and lateral slip
// all occur.
static void randomFaceState( const ModelTrack *tr, const ModelTrackRoller *ro,
                             double c[3], double uc[3], double us[3] ) {
    double runZ = tr->padZ + tr->thickness;
    c[0] = tr->linkLen * ( 0.5 + 1.5*frand2() );        // reaches past both ends
    c[1] = 0.05 * frand2();
    c[2] = runZ + ro->radius * ( 1.0 + 0.20*frand2() );
    for ( int k = 0 ; k < 3 ; k++ ) { uc[k] = 5.0 * frand2(); us[k] = 5.0 * frand2(); }
}

// ---------------------------------------------------------------------------
// 8: the partition of unity, as an identity rather than through a model. Shoes
// tile the band end to end, so one world point is at local x on one shoe,
// x - linkLen on the next and x + linkLen on the previous. Those weights must
// sum to exactly 1 and their slopes to exactly 0, or support is dropped or
// doubled at every pin -- which an explicit integrator answers with chatter.
// dw is checked against central differences too, since it is dw that carries
// the tangential half of the gradient.
// ---------------------------------------------------------------------------
static int checkFacePartition( const Model *mdl, const TrackContactParams *prm ) {
    const ModelTrackRoller *ro;
    const ModelTrack *tr = faceFixture( mdl, &ro );
    if ( !tr ) return 0;

    double L = tr->linkLen, fade = prm->faceFade * L;
    double worstSum = 0, worstSlope = 0, worstD = 0;

    for ( int i = 0 ; i <= 20000 ; i++ ) {
        double x = L * ( (double)i / 20000.0 );         // one full link
        double sum = 0, slope = 0;
        for ( int k = -3 ; k <= 3 ; k++ ) {             // every shoe that can reach
            double dw, w = trkFaceW( x + k*L, L, fade, &dw );
            sum += w; slope += dw;
        }
        if ( fabs( sum - 1.0 ) > worstSum )   worstSum   = fabs( sum - 1.0 );
        if ( fabs( slope )     > worstSlope ) worstSlope = fabs( slope );

        double h = 1e-6 * L, dw;
        trkFaceW( x, L, fade, &dw );
        double fd = ( trkFaceW( x + h, L, fade, 0 ) - trkFaceW( x - h, L, fade, 0 ) ) / ( 2*h );
        double err = fabs( dw - fd ) * fade;             // scaled by the 1/fade slope
        if ( err > worstD ) worstD = err;
    }
    if ( worstSum > 1e-12 || worstSlope > 1e-9 ) {
        printf( "contact: FAIL: face weights do not partition unity across a joint "
                "(worst |sum - 1| %.3e, worst |sum of slopes| %.3e)\n", worstSum, worstSlope );
        return 1;
    }
    if ( worstD > 1e-5 ) {
        printf( "contact: FAIL: face weight slope disagrees with central differences by "
                "%.3e relative\n", worstD );
        return 1;
    }
    printf( "contact: face: weights partition unity across every joint (worst |sum - 1| %.3e, "
            "slopes %.3e), dw matches central differences to %.3e\n",
            worstSum, worstSlope, worstD );
    return 0;
}

// ---------------------------------------------------------------------------
// 9: force = grad U for the face law, by finite difference against an
// INDEPENDENTLY WRITTEN potential -- the weight is restated here rather than
// reused, since it is the weight that carries the trap. U = w(x) * U0 is a
// gated potential, so the gradient has a tangential term U0*w'(x) beside the
// normal one; that is the same shape as the sprocket tooth law's missing
// g'(delta)*U0 term, and it is worth testing at a state rather than waiting
// for a trajectory to notice. trkFaceForce reports the force ON THE SHOE,
// which is +grad U, not -grad U.
// ---------------------------------------------------------------------------
static double faceU( const ModelTrack *tr, const ModelTrackRoller *ro,
                     const TrackContactParams *prm, const double c[3] ) {
    double runZ = tr->padZ + tr->thickness;
    double d = c[2] - runZ;
    if ( d <= 0 ) return 0;
    double delta = ro->radius - d;
    if ( delta <= 0 ) return 0;

    double fade = prm->faceFade * tr->linkLen;
    double a = ( c[0] + fade ) / ( 2*fade );
    double b = ( c[0] - tr->linkLen + fade ) / ( 2*fade );
    double sa = a <= 0 ? 0 : a >= 1 ? 1 : a*a*( 3 - 2*a );
    double sb = b <= 0 ? 0 : b >= 1 ? 1 : b*b*( 3 - 2*b );
    double w  = sa * ( 1 - sb );

    double kn   = prm->gain * tr->kn;
    double kLat = prm->gain * prm->kLatFrac * tr->kn;
    return w * ( 0.5*kn*delta*delta + 0.5*kLat*c[1]*c[1] );
}

static int checkFaceGradient( const Model *mdl, const TrackContactParams *prm ) {
    const ModelTrackRoller *ro;
    const ModelTrack *tr = faceFixture( mdl, &ro );
    if ( !tr ) return 0;

    // Elastic only: dampers and friction are not gradients of anything.
    TrackContactParams sp = *prm;
    sp.enableFriction = 0;
    sp.zeta           = 0;

    double worst = 0, scale = 0;
    int tested = 0;

    for ( int t = 0 ; t < 200000 ; t++ ) {
        double c[3], ua[3], ub[3], zero[3] = { 0, 0, 0 };
        randomFaceState( tr, ro, c, ua, ub );

        TrkFaceForce ff;
        if ( !trkFaceForce( &sp, tr, ro, c, zero, zero, &ff ) ) continue;
        // Skip the gate edges, where the one-sided cut makes a central
        // difference meaningless -- smoothness there is check 8's business.
        if ( ff.delta < 1e-4 || ff.w < 1e-3 || ff.w > 1 - 1e-3 ) continue;
        tested++;

        double h = 1e-7;
        for ( int k = 0 ; k < 3 ; k++ ) {
            double cp[3], cm[3];
            for ( int j = 0 ; j < 3 ; j++ ) { cp[j] = c[j]; cm[j] = c[j]; }
            cp[k] += h; cm[k] -= h;
            double fd  = ( faceU( tr, ro, &sp, cp ) - faceU( tr, ro, &sp, cm ) ) / ( 2*h );
            double err = fabs( ff.fCons[k] - fd );
            double s   = fabs( fd ) > 1 ? fabs( fd ) : 1.0;
            if ( fabs( fd ) > scale ) scale = fabs( fd );
            if ( err / s > worst ) worst = err / s;
            if ( err / s > 1e-5 ) {
                printf( "contact: FAIL: face elastic force disagrees with grad U in axis %d: "
                        "%.6g vs %.6g at c = (%.5f, %.5f, %.5f)\n",
                        k, ff.fCons[k], fd, c[0], c[1], c[2] );
                return 1;
            }
        }
    }
    printf( "contact: face: elastic force matches grad U to %.3e relative over %d loaded "
            "states (force scale %.3e N)\n", worst, tested, scale );
    return 0;
}

// ---------------------------------------------------------------------------
// 10: passivity of the face law, and its friction cone.
//
// Two application points means two velocities: the conservative half acts at
// the AXLE and is conjugate to uc, friction acts at the FOOT and is conjugate
// to us. And unlike the pin law, both of those are the ROLLER against the
// shoe while the reported force is still the force ON the shoe -- so the power
// the pair delivers is -f.u, not +f.u. Getting that backwards is precisely the
// bug this check exists to catch, and it is invisible to every other check
// here: the wrench balance still cancels, the gradient is still right, and the
// design pose is still clean, because at rest there is no slip to be wrong
// about.
// ---------------------------------------------------------------------------
static int checkFacePassivity( const Model *mdl, const TrackContactParams *prm ) {
    const ModelTrackRoller *ro;
    const ModelTrack *tr = faceFixture( mdl, &ro );
    if ( !tr ) return 0;

    double worst = 0, worstCone = 0;
    int tested = 0;

    for ( int t = 0 ; t < 200000 ; t++ ) {
        double c[3], uc[3], us[3], zero[3] = { 0, 0, 0 };
        randomFaceState( tr, ro, c, uc, us );

        TrkFaceForce ff, f0;
        if ( !trkFaceForce( prm, tr, ro, c, uc, us, &ff ) ) continue;
        trkFaceForce( prm, tr, ro, c, zero, zero, &f0 );
        tested++;

        double power = 0, mag = 0;
        for ( int k = 0 ; k < 3 ; k++ ) {
            power -= ( ff.fCons[k] - f0.fCons[k] ) * uc[k];
            power -= ( ff.fFric[k] - f0.fFric[k] ) * us[k];
            mag   += fabs( ff.fCons[k] - f0.fCons[k] ) * fabs( uc[k] )
                   + fabs( ff.fFric[k] - f0.fFric[k] ) * fabs( us[k] );
        }
        if ( power > worst ) worst = power;
        if ( power > 1e-12 * ( mag > 1 ? mag : 1 ) ) {
            printf( "contact: FAIL: the face law does %.6g W of work at c = (%.4f, %.4f, %.4f), "
                    "uc = (%.3f, %.3f, %.3f), us = (%.3f, %.3f, %.3f) -- it is not passive\n",
                    power, c[0], c[1], c[2], uc[0], uc[1], uc[2], us[0], us[1], us[2] );
            return 1;
        }

        // The cone, on the same states. Friction must also stay in the face
        // plane, or it would modulate the load it is itself proportional to.
        if ( ff.fFric[2] != 0 ) {
            printf( "contact: FAIL: face friction has a normal component of %.6g N\n", ff.fFric[2] );
            return 1;
        }
        double ftMag = sqrt( ff.fFric[0]*ff.fFric[0] + ff.fFric[1]*ff.fFric[1] );
        double lim   = tr->mu * ff.fN;
        double fS    = ff.fN > 1 ? ff.fN : 1.0;
        if ( ( ftMag - lim ) / fS > worstCone ) worstCone = ( ftMag - lim ) / fS;
        if ( ftMag - lim > 1e-9 * fS ) {
            printf( "contact: FAIL: face friction %.6g N exceeds the cone mu*fN = %.6g N\n",
                    ftMag, lim );
            return 1;
        }
    }
    printf( "contact: face: velocity-dependent force dissipates at every sampled state (worst "
            "power %.3e W over %d loaded states); friction inside the cone (worst excess %.3e)\n",
            worst, tested, worstCone );
    return 0;
}

// ===========================================================================
// Ground contact. The same discipline as above, but check 1 had to change
// SHAPE rather than just target: every track contact is internal, so its
// wrenches sum to zero, whereas the ground is not a body and its reactions
// have nowhere to go. There is no cancellation left to assert. What replaces
// it is a second, independently written sum of the same wrenches, plus a
// statement about WHERE in fext they are allowed to land -- because with no
// reaction body, a wrench written to the wrong body index would otherwise
// leave no trace at all in any total.
// ===========================================================================

// A random contact point straddling the terrain, and a random material
// velocity. depth < 0 is clear of the ground, > 0 is penetrating.
static void randomGroundPoint( const ModelGround *gnd, double p[3], double u[3] ) {
    double h, hx, hy;
    p[0] = 20.0 * frand2();
    p[1] = 20.0 * frand2();
    terrainAt( &gnd->terrain, p[0], p[1], &h, &hx, &hy );
    p[2] = h - 0.05 * frand2();          // +/- 50 mm about the surface
    for ( int k = 0 ; k < 3 ; k++ ) u[k] = 5.0 * frand2();
}

// ---------------------------------------------------------------------------
// G1 + G2: total wrench by two routes, where the wrenches land, and purity.
// ---------------------------------------------------------------------------
static int checkGroundBalance( const Model *mdl, const GroundContactParams *prm ) {
    static double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ];
    static double fext[ MODEL_MAX_BODIES * 6 ], fext2[ MODEL_MAX_BODIES * 6 ];
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];
    const int NTRIALS = 500;
    const double TOL = 1e-12;
    int touched = 0;
    double worstRel = 0;

    // Which bodies are allowed to receive a ground wrench at all.
    static char mayReceive[ MODEL_MAX_BODIES ];
    memset( mayReceive, 0, sizeof mayReceive );
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[t];
        if ( !tr->hasPad ) continue;
        for ( int i = 0 ; i < tr->shoes ; i++ ) mayReceive[ tr->shoe0 + i ] = 1;
    }

    for ( int t = 0 ; t < NTRIALS ; t++ ) {
        randomState( mdl, q, qd );

        double fTot[6];
        groundContactFromState( mdl, prm, q, qd, fext, fTot );

        // Route 2: re-derive every pad point from the kinematics and sum
        // [p x F; F] directly, sharing nothing with the fext bookkeeping but
        // gndPointForce itself. This is what catches a wrong moment reference
        // point -- an error a total that cancels can never see, which is
        // exactly why the track harness needed its separate drive-torque check.
        double ref[6] = { 0, 0, 0, 0, 0, 0 };
        denseKinematics( mdl, q, qd, Xa, vb );
        for ( int tr_i = 0 ; tr_i < mdl->ntrack ; tr_i++ ) {
            const ModelTrack *tr = &mdl->track[ tr_i ];
            if ( !tr->hasPad ) continue;
            for ( int i = 0 ; i < tr->shoes ; i++ ) {
                int j = tr->shoe0 + i;
                double r[3];
                xformOrigin( Xa[j], r );
                for ( int c = 0 ; c < 4 ; c++ ) {
                    double tl[3], p[3], ub[3], u[3];
                    gndPadCorner( tr, c, tl );
                    for ( int k = 0 ; k < 3 ; k++ )
                        p[k] = r[k] + Xa[j][0][k]*tl[0] + Xa[j][1][k]*tl[1] + Xa[j][2][k]*tl[2];
                    ub[0] = vb[j][3] + vb[j][1]*tl[2] - vb[j][2]*tl[1];
                    ub[1] = vb[j][4] + vb[j][2]*tl[0] - vb[j][0]*tl[2];
                    ub[2] = vb[j][5] + vb[j][0]*tl[1] - vb[j][1]*tl[0];
                    for ( int k = 0 ; k < 3 ; k++ )
                        u[k] = Xa[j][0][k]*ub[0] + Xa[j][1][k]*ub[1] + Xa[j][2][k]*ub[2];
                    GndPointForce pf;
                    if ( !gndPointForce( prm, &mdl->ground, p, u, &pf ) ) continue;
                    ref[0] += p[1]*pf.f[2] - p[2]*pf.f[1];
                    ref[1] += p[2]*pf.f[0] - p[0]*pf.f[2];
                    ref[2] += p[0]*pf.f[1] - p[1]*pf.f[0];
                    for ( int k = 0 ; k < 3 ; k++ ) ref[3+k] += pf.f[k];
                }
            }
        }

        double sum[6] = { 0, 0, 0, 0, 0, 0 }, scale[6] = { 0, 0, 0, 0, 0, 0 };
        for ( int i = 0 ; i < mdl->nbod ; i++ )
            for ( int k = 0 ; k < 6 ; k++ ) {
                sum[k]   += fext[ i*6 + k ];
                scale[k] += fabs( fext[ i*6 + k ] );
            }
        for ( int i = 0 ; i < mdl->nbod * 6 ; i++ ) if ( fext[i] != 0.0 ) touched = 1;

        for ( int k = 0 ; k < 6 ; k++ ) {
            double sc  = scale[k] > 1 ? scale[k] : 1.0;
            double rel = fabs( sum[k] - ref[k] ) / sc;
            if ( rel > worstRel ) worstRel = rel;
            if ( rel > TOL ) {
                printf( "contact: FAIL: ground wrench disagrees with the direct sum at trial %d: "
                        "component %d is %.6g via fext against %.6g direct (relative %.3e)\n",
                        t, k, sum[k], ref[k], rel );
                return 1;
            }
            // fTotal is the engine's own report of what it applied, and
            // src/sim/main.c's momentum balance is measured against it.
            if ( fabs( sum[k] - fTot[k] ) > TOL * sc ) {
                printf( "contact: FAIL: the reported total ground wrench does not match what "
                        "landed in fext at trial %d: component %d, %.6g reported, %.6g applied\n",
                        t, k, fTot[k], sum[k] );
                return 1;
            }
        }

        // Ground wrenches must land only on bodies that carry pads. With no
        // reaction body there is nothing to cancel against, so a wrench written
        // to the wrong index is invisible in every total above.
        for ( int i = 0 ; i < mdl->nbod ; i++ ) {
            if ( mayReceive[i] ) continue;
            for ( int k = 0 ; k < 6 ; k++ ) if ( fext[ i*6 + k ] != 0.0 ) {
                printf( "contact: FAIL: body %d carries no ground pad but received a wrench "
                        "(component %d = %.6g) at trial %d\n", i, k, fext[ i*6 + k ], t );
                return 1;
            }
        }

        groundContactFromState( mdl, prm, q, qd, fext2, NULL );
        for ( int i = 0 ; i < mdl->nbod * 6 ; i++ ) if ( fext[i] != fext2[i] ) {
            printf( "contact: FAIL: ground engine is not a pure function of (q,qd): "
                    "fext[%d] = %.17g then %.17g at trial %d\n", i, fext[i], fext2[i], t );
            return 1;
        }
    }

    printf( "contact: ground: %d random states: wrench matches an independent sum to %.3e "
            "relative, lands only on padded bodies, output deterministic%s\n",
            NTRIALS, worstRel, touched ? "" : " (no pads in range -- model has no ground)" );
    return 0;
}

// ---------------------------------------------------------------------------
// G3 + G4: smoothness, vertically through the surface and horizontally across
// it. The vertical sweep is the counterpart of the track harness's radial one.
//
// The HORIZONTAL sweep is the new one, and it is aimed at terrain that does not
// exist yet: dragging a point along the surface must not make the force jump.
// On a plane that is trivially true, which is the point of writing it now --
// when a sampled grid arrives, bilinear interpolation gives C0 heights and
// therefore a normal direction that snaps at every cell edge, and this check is
// what will reject it and force a bicubic.
// ---------------------------------------------------------------------------
static int checkGroundSmoothness( const Model *mdl, const GroundContactParams *prm ) {
    if ( !mdl->ground.has ) return 0;
    const int N = 20000;
    double u[3] = { 0.3, -0.2, 0.4 };

    // Vertical: 100 mm of travel straddling the surface.
    double h, hx, hy, worstV = 0, prev[3] = { 0, 0, 0 };
    terrainAt( &mdl->ground.terrain, 1.0, 0.5, &h, &hx, &hy );
    for ( int i = 0 ; i <= N ; i++ ) {
        double p[3] = { 1.0, 0.5, h + 0.05 - 0.10 * i / N };
        GndPointForce pf;
        gndPointForce( prm, &mdl->ground, p, u, &pf );
        if ( i ) {
            double d = 0;
            for ( int k = 0 ; k < 3 ; k++ ) d += ( pf.f[k] - prev[k] ) * ( pf.f[k] - prev[k] );
            d = sqrt( d );
            if ( d > worstV ) worstV = d;
        }
        for ( int k = 0 ; k < 3 ; k++ ) prev[k] = pf.f[k];
    }

    // Horizontal: 20 m of travel at a fixed depth, so only the terrain varies.
    double worstH = 0;
    for ( int i = 0 ; i <= N ; i++ ) {
        double x = -10.0 + 20.0 * i / N, hh, hxx, hyy;
        terrainAt( &mdl->ground.terrain, x, 0.5, &hh, &hxx, &hyy );
        double p[3] = { x, 0.5, hh - 0.01 };       // 10 mm down, always in contact
        GndPointForce pf;
        gndPointForce( prm, &mdl->ground, p, u, &pf );
        if ( i ) {
            double d = 0;
            for ( int k = 0 ; k < 3 ; k++ ) d += ( pf.f[k] - prev[k] ) * ( pf.f[k] - prev[k] );
            d = sqrt( d );
            if ( d > worstH ) worstH = d;
        }
        for ( int k = 0 ; k < 3 ; k++ ) prev[k] = pf.f[k];
    }

    // The vertical bound is the elastic slope: crossing 100 mm in N steps at
    // stiffness kn cannot legitimately step by more than kn*100mm/N, and the
    // damper adds its own ramp. Allow 4x that; a genuine discontinuity is
    // orders of magnitude over, not a factor of a few.
    double vBound = 4.0 * ( mdl->ground.kn * 0.10 / N
                          + 2.0 * prm->zeta * sqrt( mdl->ground.kn * prm->mEff ) * 0.10 / N );
    if ( worstV > vBound && vBound > 0 ) {
        printf( "contact: FAIL: the vertical force sweep steps by %.4g N against a %.4g N bound "
                "-- the ground law is discontinuous\n", worstV, vBound );
        return 1;
    }
    printf( "contact: ground: sweeps smooth: vertical max step %.3e N (bound %.3e), "
            "horizontal %.3e N over 20 m\n", worstV, vBound, worstH );
    return 0;
}

// ---------------------------------------------------------------------------
// G5: the elastic force is exactly -grad U, checked end to end.
//
// This is the check the whole vertical-gap formulation exists to pass, and it
// is aimed squarely at the failure that cost the sprocket tooth law two rounds
// of debugging: a force that is only PART of the gradient of its own potential
// leaks energy around any closed path, and no pointwise test of either the
// force or the potential alone can see it.
//
// Method: translate the entire model by +/-eps and central-difference
// groundContactPotential against the total force from groundContactFext.
// Translating every FLOATING body's position slice translates the whole model,
// since every other body's pose comes from tree offsets beneath one. Evaluated
// at qd = 0, which kills the damper and friction and leaves exactly the elastic
// part -- the only part that is supposed to be a gradient.
// ---------------------------------------------------------------------------
static int checkGroundGradient( const Model *mdl, const GroundContactParams *prm ) {
    if ( !mdl->ground.has ) return 0;
    static double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ], qp[ MODEL_MAX_Q ];
    static double fext[ MODEL_MAX_BODIES * 6 ];
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];
    const double EPS = 1e-6;
    double worst = 0, worstScale = 0;
    int tested = 0;

    for ( int t = 0 ; t < 200 ; t++ ) {
        randomState( mdl, q, qd );
        memset( qd, 0, mdl->nv * sizeof *qd );        // elastic part only

        double fTot[6];
        groundContactFromState( mdl, prm, q, qd, fext, fTot );
        double mag = fabs( fTot[3] ) + fabs( fTot[4] ) + fabs( fTot[5] );
        if ( mag < 1.0 ) continue;                     // nothing in contact
        tested++;

        for ( int k = 0 ; k < 3 ; k++ ) {
            double U[2];
            for ( int s = 0 ; s < 2 ; s++ ) {
                memcpy( qp, q, mdl->nq * sizeof *qp );
                double d = s ? EPS : -EPS;
                for ( int i = 0 ; i < mdl->nbod ; i++ )
                    if ( mdl->body[i].jtype == JOINT_FLOATING )
                        qp[ mdl->qidx[i] + 4 + k ] += d;
                denseKinematics( mdl, qp, qd, Xa, vb );
                U[s] = groundContactPotential( mdl, prm, Xa );
            }
            double dUdx = ( U[1] - U[0] ) / ( 2.0 * EPS );
            double e    = fabs( -dUdx - fTot[ 3 + k ] );
            if ( e / mag > worst ) { worst = e / mag; worstScale = mag; }
        }
    }

    if ( !tested ) {
        printf( "contact: ground: gradient check found no loaded states to sample\n" );
        return 0;
    }
    // Central differences at eps 1e-6 on a stiff quadratic: the floor is
    // roundoff in U divided by eps, a few parts in 1e-7 of the force scale.
    if ( worst > 1e-5 ) {
        printf( "contact: FAIL: the elastic ground force is not -grad U (worst relative "
                "mismatch %.3e over %d loaded states) -- the law is not conservative\n",
                worst, tested );
        return 1;
    }
    printf( "contact: ground: elastic force matches -grad U to %.3e relative over %d loaded "
            "states (force scale %.4g N)\n", worst, tested, worstScale );
    return 0;
}

// ---------------------------------------------------------------------------
// G6: the terrain's analytic gradient against central differences of its own
// height. Cheap, and it catches a bad dh/dx in every terrain kind that will
// ever be added -- including the ones this file cannot otherwise reach, since
// the gradient is what the conservativeness argument in terrain.h rests on.
// ---------------------------------------------------------------------------
static int checkTerrainGradient( const Model *mdl ) {
    if ( !mdl->ground.has ) return 0;
    const ModelTerrain *tn = &mdl->ground.terrain;
    const double EPS = 1e-5;
    double worst = 0;

    for ( int t = 0 ; t < 5000 ; t++ ) {
        double x = 50.0 * frand2(), y = 50.0 * frand2();
        double h, hx, hy, a, b, gx, gy;
        if ( !terrainAt( tn, x, y, &h, &hx, &hy ) ) return 0;
        terrainAt( tn, x + EPS, y, &a, &gx, &gy );
        terrainAt( tn, x - EPS, y, &b, &gx, &gy );
        double fdx = ( a - b ) / ( 2.0 * EPS );
        terrainAt( tn, x, y + EPS, &a, &gx, &gy );
        terrainAt( tn, x, y - EPS, &b, &gx, &gy );
        double fdy = ( a - b ) / ( 2.0 * EPS );

        double sc = 1.0 + fabs( hx ) + fabs( hy );
        double e  = ( fabs( fdx - hx ) + fabs( fdy - hy ) ) / sc;
        if ( e > worst ) worst = e;
        if ( e > 1e-6 ) {
            printf( "contact: FAIL: terrain gradient disagrees with finite differences at "
                    "(%.4f, %.4f): analytic (%.6g, %.6g), differenced (%.6g, %.6g)\n",
                    x, y, hx, hy, fdx, fdy );
            return 1;
        }
    }
    printf( "contact: ground: terrain gradient matches central differences to %.3e relative\n", worst );
    return 0;
}

// ---------------------------------------------------------------------------
// G7 + G8: passivity, and the friction cone.
//
// Passivity is the same statement as the track harness's, and here it is
// provable rather than empirical -- power = -(kn*g + c*ramp*gd)*gd, whose
// second term is a square. The check is still worth running: it is what would
// catch the one-sided clamp being dropped, or gd losing its sign.
//
// The cone is separate because it is about the friction term specifically:
// |f_t| <= mu*|f_n|, and f_t must be perpendicular to the surface normal. A
// tangential force with a normal component would quietly add to or subtract
// from the load it is itself proportional to.
// ---------------------------------------------------------------------------
static int checkGroundPassivityAndCone( const Model *mdl, const GroundContactParams *prm ) {
    if ( !mdl->ground.has ) return 0;
    double worstPower = 0, worstCone = 0, worstPerp = 0;

    for ( int t = 0 ; t < 200000 ; t++ ) {
        double p[3], u[3], uz[3] = { 0, 0, 0 };
        randomGroundPoint( &mdl->ground, p, u );

        GndPointForce pf, p0;
        if ( !gndPointForce( prm, &mdl->ground, p, u, &pf ) ) continue;
        gndPointForce( prm, &mdl->ground, p, uz, &p0 );

        double power = 0, mag = 0;
        for ( int k = 0 ; k < 3 ; k++ ) {
            power += ( pf.f[k] - p0.f[k] ) * u[k];
            mag   += fabs( pf.f[k] - p0.f[k] ) * fabs( u[k] );
        }
        if ( power > worstPower ) worstPower = power;
        if ( power > 1e-12 * ( mag > 1 ? mag : 1 ) ) {
            printf( "contact: FAIL: the ground law does %.6g W of work at p = (%.4f, %.4f, %.4f), "
                    "u = (%.3f, %.3f, %.3f) -- it is not passive\n",
                    power, p[0], p[1], p[2], u[0], u[1], u[2] );
            return 1;
        }

        // Split the total back into normal and tangential about the unit
        // normal the law itself reported.
        double fn = 0;
        for ( int k = 0 ; k < 3 ; k++ ) fn += pf.f[k] * pf.n[k];
        double ft[3], ftMag = 0;
        for ( int k = 0 ; k < 3 ; k++ ) { ft[k] = pf.f[k] - fn * pf.n[k]; ftMag += ft[k]*ft[k]; }
        ftMag = sqrt( ftMag );

        double perp = 0;
        for ( int k = 0 ; k < 3 ; k++ ) perp += ft[k] * pf.n[k];
        double fScale = pf.fN > 1 ? pf.fN : 1.0;
        if ( fabs( perp ) / fScale > worstPerp ) worstPerp = fabs( perp ) / fScale;

        // mu*fN is the cone; the regularization only ever shrinks the friction
        // below it, never above, so this is a hard bound with no tolerance
        // beyond roundoff.
        double lim = mdl->ground.mu * pf.fN;
        double over = ftMag - lim;
        if ( over / fScale > worstCone ) worstCone = over / fScale;
        if ( over > 1e-9 * fScale ) {
            printf( "contact: FAIL: friction %.6g N exceeds the cone mu*fN = %.6g N at "
                    "p = (%.4f, %.4f, %.4f)\n", ftMag, lim, p[0], p[1], p[2] );
            return 1;
        }
    }
    printf( "contact: ground: passivity holds everywhere (worst power %.3e W); friction stays "
            "inside the cone (worst excess %.3e relative) and normal to the surface (%.3e)\n",
            worstPower, worstCone, worstPerp );
    return 0;
}

int main( int argc, char **argv ) {
    const char *path = argc > 1 ? argv[ 1 ] : "../models/m113.mech";
    Model mdl;
    char err[ 256 ];
    if ( modelLoad( path, &mdl, err, sizeof err ) ) {
        fprintf( stderr, "contact: %s\n", err );
        return 1;
    }
    if ( mdl.ntrack && trackContactCheckModel( &mdl, err, sizeof err ) ) {
        fprintf( stderr, "contact: FAIL: %s\n", err );
        return 1;
    }
    if ( groundContactCheckModel( &mdl, err, sizeof err ) ) {
        fprintf( stderr, "contact: FAIL: %s\n", err );
        return 1;
    }

    TrackContactParams prm;
    trackContactDefaults( &mdl, &prm );
    GroundContactParams gprm;
    groundContactDefaults( &mdl, &gprm );

    printf( "contact: %s: %d track band(s), %s\n", path, mdl.ntrack,
            mdl.ground.has ? "ground declared" : "no ground" );
    if ( checkBalanceAndPurity( &mdl, &prm ) ) return 1;
    if ( checkSmoothness( &mdl, &prm ) )       return 1;
    if ( checkToothLaw( &mdl, &prm ) )         return 1;
    if ( checkDesignPose( &mdl, &prm ) )       return 1;
    if ( checkDriveTorque( &mdl, &prm ) )      return 1;
    if ( checkPassivity( &mdl, &prm ) )        return 1;

    // The road-wheel face law, which none of the above reaches: every one of
    // them drives trkPairForce, and road wheels stopped going through it when
    // the face primitive landed. Inert on a band with no non-'end' rollers.
    if ( checkFacePartition( &mdl, &prm ) )    return 1;
    if ( checkFaceGradient( &mdl, &prm ) )     return 1;
    if ( checkFacePassivity( &mdl, &prm ) )    return 1;

    // Ground. Every one of these is inert on a model with no 'ground' block,
    // except the balance check, which asserts the engine produces nothing at
    // all -- the same thing the track checks assert for a model with no band.
    if ( checkGroundBalance( &mdl, &gprm ) )          return 1;
    if ( checkGroundSmoothness( &mdl, &gprm ) )       return 1;
    if ( checkGroundGradient( &mdl, &gprm ) )         return 1;
    if ( checkTerrainGradient( &mdl ) )               return 1;
    if ( checkGroundPassivityAndCone( &mdl, &gprm ) ) return 1;

    printf( "contact: PASS\n" );
    return 0;
}
