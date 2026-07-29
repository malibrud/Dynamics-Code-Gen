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

    TrackContactParams prm;
    trackContactDefaults( &mdl, &prm );

    printf( "contact: %s: %d track band(s)\n", path, mdl.ntrack );
    if ( checkBalanceAndPurity( &mdl, &prm ) ) return 1;
    if ( checkSmoothness( &mdl, &prm ) )       return 1;
    if ( checkToothLaw( &mdl, &prm ) )         return 1;
    if ( checkDesignPose( &mdl, &prm ) )       return 1;
    if ( checkDriveTorque( &mdl, &prm ) )      return 1;
    if ( checkPassivity( &mdl, &prm ) )        return 1;

    printf( "contact: PASS\n" );
    return 0;
}
