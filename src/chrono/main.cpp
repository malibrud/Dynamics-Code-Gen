// Chrono cross-check: builds the mechanism described by a .mech model file
// (src/model/model.h) as a Project Chrono multibody system and compares
// Chrono's constraint-based accelerations against this project's dense RBDA
// reference (src/validate/dense.h) on random states.
//
// Unlike src/validate (which shares the spatial-vector conventions with the
// generator by construction), Chrono is a completely independent dynamics
// engine -- maximal-coordinate bodies, quaternion states, joint constraints
// solved as a KKT system -- so agreement here is an external check of the
// .mech conventions themselves, not just of the code generation.
//
// Convention notes (matching src/model/model.h, encoded in fk() below):
//   - q is the physical joint coordinate: a revolute q means a right-handed
//     rotation of +q about the +axis, and qd = dq/dt. No sign games.
//     This checker is what established that: before the joint transforms
//     were fixed to E = R(q)^T, the placement and the rate were
//     opposite-handed (q ran one way, qd the other), and the velocity-only
//     trials below came back as exactly -qdd. src/validate could not see it,
//     because the dense reference and the generator shared the convention.
//   - The floating-joint quaternion is the ordinary body attitude
//     (body-to-base) -- the same convention Chrono itself uses.
//   - fext[i] is a spatial force in base coordinates: [ moment about the
//     base origin ; force ].
//   - For a floating body, qd/qdd slots are the body's spatial velocity /
//     its time derivative in the body's own frame: [ w ; v_origin ].
//
// Usage:
//   chrono_check <model.mech>                 cross-check vs dense reference
//   chrono_check <model.mech> --sim <T> <dt>  time-step the Chrono model and
//                                             print body positions
//
// Build: see build_chrono.bat (needs the Chrono SDK and Eigen).

#include <cstdio>
#include <cmath>
#include <vector>
#include <memory>

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChLinkLock.h"
#include "chrono/solver/ChDirectSolverLS.h"

#include <Eigen/Dense>

#include "../validate/dense.h" // fd_dense + frand (includes ../model/model.h)

// The generated implementation under test, compiled from build/dynamics_<model>.c
// (C, hence extern "C"). Linking it here means Chrono validates the artifact
// this project actually ships, rather than reaching it transitively through
// the dense reference.
extern "C" {
    void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *out_qdd );
    extern const int fwd_dyn_crba_nq;
    extern const int fwd_dyn_crba_nv;
}

using namespace chrono;

// ---------------------------------------------------------------------------
// Forward kinematics of the .mech model: world pose and velocity of every
// body frame, in the conventions summarized at the top of this file. This is
// what places/launches the Chrono bodies for a given (q, qd).
// ---------------------------------------------------------------------------
struct BodyKin {
    ChMatrix33<> R;    // world orientation of the body (mech) frame
    ChVector3d r;      // world position of the body frame origin
    ChVector3d w;      // world angular velocity
    ChVector3d vO;     // world velocity of the body-fixed point at the origin
};

static ChVector3d v3( const double *p ) { return ChVector3d( p[0], p[1], p[2] ); }

static ChVector3d axis3( char axis ) {
    if ( axis == 'x' ) return ChVector3d( 1, 0, 0 );
    if ( axis == 'y' ) return ChVector3d( 0, 1, 0 );
    return ChVector3d( 0, 0, 1 );
}

static void fk( const Model *mdl, const double *q, const double *qd, std::vector<BodyKin> &kin ) {
    kin.resize( mdl->nbod );
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];
        const double *qs = q + mdl->qidx[ i ];
        const double *vs = qd + mdl->vidx[ i ];

        // Pre-joint frame P: parent frame translated by the tree offset.
        BodyKin P;
        if ( b->parent == -1 ) {
            P.R = ChMatrix33<>( 1 );
            P.r = v3( b->tree );
            P.w = P.vO = ChVector3d( 0, 0, 0 );
        } else {
            const BodyKin &p = kin[ b->parent ];
            ChVector3d treeW = p.R * v3( b->tree );
            P.R = p.R;
            P.r = p.r + treeW;
            P.w = p.w;
            P.vO = p.vO + Vcross( p.w, treeW );
        }

        BodyKin &k = kin[ i ];
        if ( b->jtype == JOINT_FLOATING ) {
            // .mech quaternion is the ordinary body attitude, like Chrono's;
            // normalized here because the model accepts any nonzero one.
            double w = qs[0], x = qs[1], y = qs[2], z = qs[3];
            double n = 1.0 / sqrt( w*w + x*x + y*y + z*z );
            ChQuaterniond att( w*n, x*n, y*n, z*n );
            k.R = ChMatrix33<>( att );
            k.r = P.r + P.R * ChVector3d( qs[4], qs[5], qs[6] );
            // qd = spatial velocity in the body's own frame [ w ; v_origin ].
            k.w  = k.R * ChVector3d( vs[0], vs[1], vs[2] );
            k.vO = k.R * ChVector3d( vs[3], vs[4], vs[5] );
        } else if ( b->jtype == JOINT_REVOLUTE ) {
            // Right-handed rotation of +q about +axis, spinning at +qd.
            ChVector3d a = axis3( b->axis );
            double half = 0.5 * qs[0];
            ChQuaterniond rot( cos( half ), sin( half ) * a.x(), sin( half ) * a.y(), sin( half ) * a.z() );
            k.R = P.R * ChMatrix33<>( rot );
            k.r = P.r;
            k.w = P.w + ( P.R * a ) * vs[0];
            k.vO = P.vO;
        } else { // JOINT_PRISMATIC
            ChVector3d aW = P.R * axis3( b->axis );
            k.R = P.R;
            k.r = P.r + aW * qs[0];
            k.w = P.w;
            k.vO = P.vO + Vcross( P.w, aW * qs[0] ) + aW * vs[0];
        }
    }
}

// ---------------------------------------------------------------------------
// Advances the model's position coordinates by dt given qd -- the reference
// integration rule for .mech coordinates:
//   - single-DOF joints: q += qd*dt, because qd is exactly dq/dt.
//   - floating bodies: qd carries the body's spatial velocity in its OWN
//     frame while q carries a base-frame quaternion + position, so the rates
//     are Qdot = 0.5 * Q x (0, w_body) and rdot = R(Q) * v_body.
// ---------------------------------------------------------------------------
static void integrateQ( const Model *mdl, const double *q, const double *qd, double dt, double *qOut ) {
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];
        const double *qs = q + mdl->qidx[ i ];
        const double *vs = qd + mdl->vidx[ i ];
        double *out = qOut + mdl->qidx[ i ];

        if ( b->jtype != JOINT_FLOATING ) {
            out[0] = qs[0] + vs[0] * dt;
            continue;
        }
        ChQuaterniond Q( qs[0], qs[1], qs[2], qs[3] );
        Q.Normalize();
        ChQuaterniond Qd = Q * ChQuaterniond( 0, vs[0], vs[1], vs[2] ) * 0.5;
        ChQuaterniond Qn = Q + Qd * dt;
        Qn.Normalize();
        ChVector3d rdot = ChMatrix33<>( Q ) * ChVector3d( vs[3], vs[4], vs[5] );
        out[0] = Qn.e0(); out[1] = Qn.e1(); out[2] = Qn.e2(); out[3] = Qn.e3();
        out[4] = qs[4] + rdot.x() * dt;
        out[5] = qs[5] + rdot.y() * dt;
        out[6] = qs[6] + rdot.z() * dt;
    }
}

// Asserts the invariant that qd really is the time derivative of q: the body
// velocities fk() derives from qd must equal d/dt of the body poses fk()
// derives from q, with q advanced by integrateQ.
//
// This is the property that was silently violated while the joint transforms
// used E = R(q): q ran one way and qd the other, so an integrator stepping
// q += qd*dt drove every mechanism backwards. Neither src/validate (dense and
// generated share the convention) nor the dynamics cross-check below (which
// compares accelerations at a fixed state) can see it, so it is checked here
// directly -- no dynamics and no Chrono solver involved, just kinematics.
static int checkQdConsistency( const Model *mdl, const char *path ) {
    const double DT = 1e-6;
    const double TOL = 1e-6;
    double maxErr = 0;
    int worstTrial = -1;
    std::vector<BodyKin> k0, kp, km;

    for ( int t = 0 ; t < 50 ; t++ ) {
        double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ], qp[ MODEL_MAX_Q ], qm[ MODEL_MAX_Q ];
        for ( int k = 0 ; k < mdl->nq ; k++ ) q[k]  = 3.0 * frand();
        for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] = 2.0 * frand();

        integrateQ( mdl, q, qd, DT, qp );
        integrateQ( mdl, q, qd, -DT, qm );
        fk( mdl, q,  qd, k0 );
        fk( mdl, qp, qd, kp );
        fk( mdl, qm, qd, km );

        for ( int i = 0 ; i < mdl->nbod ; i++ ) {
            // Linear: central difference of the body origin vs fk's vO.
            ChVector3d vFD = ( kp[i].r - km[i].r ) / ( 2.0 * DT );
            double eLin = ( vFD - k0[i].vO ).Length() / fmax( 1.0, k0[i].vO.Length() );

            // Angular: dR/dt = skew(w) * R, so skew(w) = dR/dt * R^T.
            Eigen::Matrix3d dR = ( Eigen::Matrix3d( kp[i].R ) - Eigen::Matrix3d( km[i].R ) ) / ( 2.0 * DT );
            Eigen::Matrix3d wHat = dR * Eigen::Matrix3d( k0[i].R ).transpose();
            ChVector3d wFD( wHat(2,1), wHat(0,2), wHat(1,0) );
            double eAng = ( wFD - k0[i].w ).Length() / fmax( 1.0, k0[i].w.Length() );

            double e = fmax( eLin, eAng );
            if ( e > maxErr ) { maxErr = e; worstTrial = t; }
        }
    }

    printf( "chrono_check: %s: qd-vs-dq/dt consistency, max scaled error = %.3e (worst at trial %d)\n",
            path, maxErr, worstTrial );
    if ( maxErr < TOL ) return 0;
    printf( "chrono_check: FAIL: qd is not the time derivative of q (tolerance %.1e)\n", TOL );
    return 1;
}

// ---------------------------------------------------------------------------
// Chrono model: one ChBody per .mech body (Chrono body frame = CG frame,
// axes parallel to the mech body frame), one ChLinkLock joint per tree edge,
// one force accumulator per body for tau / fext. Built once per model; each
// trial repositions the bodies and refills the accumulators.
// ---------------------------------------------------------------------------
struct ChronoBody {
    std::shared_ptr<ChBody> body;
    unsigned int accum;    // force/torque accumulator index
    ChVector3d cg;         // mech body origin is at local -cg in the ChBody frame
};

struct ChronoModel {
    ChSystemNSC sys;
    std::shared_ptr<ChBody> ground;
    std::vector<ChronoBody> bod;
};

// Rotation mapping local z onto the joint axis, so ChLinkLock joints (which
// rotate about / slide along marker z) line up with .mech joint axes.
static ChQuaterniond zToAxis( char axis ) {
    if ( axis == 'x' ) return ChQuaterniond(  0.5, 0.5, 0.5, 0.5 );
    if ( axis == 'y' ) return ChQuaterniond( -0.5, 0.5, 0.5, 0.5 );
    return ChQuaterniond( 1, 0, 0, 0 );
}

static void buildChronoModel( const Model *mdl, ChronoModel &cm ) {
    cm.sys.SetGravitationalAcceleration( v3( mdl->gravity ) );
    // Bilateral-only system: a direct sparse solver gives machine-precision
    // constraint solves, which the acceleration comparison needs.
    cm.sys.SetSolver( chrono_types::make_shared<ChSolverSparseQR>() );

    cm.ground = chrono_types::make_shared<ChBody>();
    cm.ground->SetFixed( true );
    cm.sys.AddBody( cm.ground );

    // Bodies, placed at the zero configuration of the model.
    static double q0[ MODEL_MAX_Q ], qd0[ MODEL_MAX_V ];
    for ( int k = 0 ; k < mdl->nq ; k++ ) q0[k] = 0;
    for ( int k = 0 ; k < mdl->nv ; k++ ) qd0[k] = 0;
    for ( int i = 0 ; i < mdl->nbod ; i++ )
        if ( mdl->body[i].jtype == JOINT_FLOATING ) q0[ mdl->qidx[i] ] = 1.0;
    std::vector<BodyKin> kin;
    fk( mdl, q0, qd0, kin );

    cm.bod.resize( mdl->nbod );
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];
        ChronoBody &cb = cm.bod[ i ];
        cb.cg = v3( b->cg );

        cb.body = chrono_types::make_shared<ChBody>();
        cb.body->SetMass( b->mass );
        // .mech inertia entries are the actual tensor entries about the CG
        // (xx yy zz xy yz zx); Chrono wants (Ixx,Iyy,Izz) and (Ixy,Ixz,Iyz)
        // as tensor entries too.
        cb.body->SetInertiaXX( ChVector3d( b->inertia[0], b->inertia[1], b->inertia[2] ) );
        cb.body->SetInertiaXY( ChVector3d( b->inertia[3], b->inertia[5], b->inertia[4] ) );
        cb.body->SetPos( kin[i].r + kin[i].R * cb.cg );
        cb.body->SetRot( kin[i].R.GetQuaternion() );
        cb.accum = cb.body->AddAccumulator();
        cm.sys.AddBody( cb.body );

        if ( b->jtype == JOINT_FLOATING ) continue; // free body, no joint

        std::shared_ptr<ChBody> parent = b->parent == -1 ? cm.ground : cm.bod[ b->parent ].body;
        ChVector3d parCg = b->parent == -1 ? ChVector3d( 0, 0, 0 ) : cm.bod[ b->parent ].cg;
        ChQuaterniond A = zToAxis( b->axis );
        // Joint frame: at the child body origin (local -cg in the child, at
        // tree - cg in the parent), z axis along the joint axis. Marker
        // frames are local, so they stay valid when trials reposition bodies.
        ChFramed frameChild ( -cb.cg, A );
        ChFramed frameParent( v3( b->tree ) - parCg, A );

        std::shared_ptr<ChLinkLock> link;
        if ( b->jtype == JOINT_REVOLUTE )
            link = chrono_types::make_shared<ChLinkLockRevolute>();
        else
            link = chrono_types::make_shared<ChLinkLockPrismatic>();
        link->Initialize( cb.body, parent, true, frameChild, frameParent );
        cm.sys.AddLink( link );
    }
}

// Repositions all Chrono bodies to the state (q, qd) and loads the
// accumulators with tau and fext.
static void setChronoState( const Model *mdl, ChronoModel &cm,
                            const double *q, const double *qd,
                            const double *tau, const double *fext,
                            std::vector<BodyKin> &kin ) {
    fk( mdl, q, qd, kin );

    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];
        ChronoBody &cb = cm.bod[ i ];
        const BodyKin &k = kin[ i ];
        ChVector3d cgW = k.R * cb.cg;

        cb.body->SetPos( k.r + cgW );
        cb.body->SetRot( k.R.GetQuaternion() );
        cb.body->SetPosDt( k.vO + Vcross( k.w, cgW ) );  // CG velocity
        cb.body->SetAngVelParent( k.w );

        cb.body->EmptyAccumulator( cb.accum );

        // Joint generalized force tau: an equal-and-opposite wrench pair
        // between child and parent, conjugate to this joint's coordinate.
        const double *ts = tau + mdl->vidx[ i ];
        std::shared_ptr<ChBody> parent = b->parent == -1 ? cm.ground : cm.bod[ b->parent ].body;
        if ( b->jtype == JOINT_FLOATING ) {
            // tau = [ moment about body origin ; force ], in body coords.
            cb.body->AccumulateTorque( cb.accum, ChVector3d( ts[0], ts[1], ts[2] ), true );
            cb.body->AccumulateForce ( cb.accum, ChVector3d( ts[3], ts[4], ts[5] ), -cb.cg, true );
        } else if ( b->jtype == JOINT_REVOLUTE ) {
            // tau is conjugate to the motion subspace S = +a: a couple about
            // the +axis.
            ChVector3d T = axis3( b->axis ) * ts[0]; // same components in both frames
            cb.body->AccumulateTorque( cb.accum, T, true );
            if ( b->parent != -1 )
                cm.bod[ b->parent ].body->AccumulateTorque( cm.bod[ b->parent ].accum, -T, true );
        } else { // JOINT_PRISMATIC
            const BodyKin &Pk = b->parent == -1 ? BodyKin{ ChMatrix33<>( 1 ), ChVector3d(0,0,0), ChVector3d(0,0,0), ChVector3d(0,0,0) } : kin[ b->parent ];
            ChVector3d F = ( Pk.R * axis3( b->axis ) ) * ts[0];
            cb.body->AccumulateForce( cb.accum, F, k.r, false ); // at the joint, world coords
            if ( b->parent != -1 )
                cm.bod[ b->parent ].body->AccumulateForce( cm.bod[ b->parent ].accum, -F, k.r, false );
        }

        // External spatial force, base coords: [ n0 = moment about the base
        // origin ; f ]. Equivalent to force f acting at the world origin
        // plus a couple n0.
        if ( b->hasFext ) {
            ChVector3d n0( fext[ i*6 + 0 ], fext[ i*6 + 1 ], fext[ i*6 + 2 ] );
            ChVector3d f ( fext[ i*6 + 3 ], fext[ i*6 + 4 ], fext[ i*6 + 5 ] );
            cb.body->AccumulateForce ( cb.accum, f, ChVector3d( 0, 0, 0 ), false );
            cb.body->AccumulateTorque( cb.accum, n0, false );
        }
    }
}

// Solves for the exact accelerations of the Chrono system at its current
// state and scatters them into the bodies.
//
// Chrono's own DoAssembly() obtains accelerations by finite-differencing a
// small integration step, which caps its accuracy around 1e-5 and (with the
// linearized stepper) loses the centrifugal Cq_dot*v terms outright. For a
// tight cross-check we instead assemble the acceleration-level KKT system
//
//     [ M   Cq^T ] [ a      ]   [ F  ]
//     [ Cq  0    ] [ lambda ] = [ Qc ] ,   Qc = -d/dt(Cq) * v
//
// from Chrono's exported mass matrix, constraint Jacobian, and applied+
// gyroscopic force vector, and solve it densely with Eigen. Qc is obtained
// by central-differencing Cq*v along the state flow x(t) (the only quantity
// Chrono does not expose analytically); with a 1e-6 step that is accurate to
// ~1e-10, far below the comparison tolerance.
static bool solveAccelerations( ChSystemNSC &sys ) {
    sys.Setup();
    sys.Update();
    // The Get*Matrix functions read from the system descriptor, which is
    // only populated by the solver path; inject variables/constraints and
    // load the constraint Jacobians explicitly.
    sys.DescriptorPrepareInject();
    sys.LoadConstraintJacobians();

    int n  = (int)sys.GetNumCoordsVelLevel();
    int m  = (int)sys.GetNumConstraints();
    int np = (int)sys.GetNumCoordsPosLevel();

    ChState x( np, &sys );
    ChStateDelta v( n, &sys );
    double T;
    sys.StateGather( x, v, T );

    Eigen::VectorXd Qc = Eigen::VectorXd::Zero( m );
    Eigen::MatrixXd Cqd;
    if ( m > 0 ) {
        const double DEL = 1e-6;
        ChSparseMatrix Cqp( m, n ), Cqm( m, n );
        ChState xp( np, &sys ), xm( np, &sys );
        ChStateDelta dx( n, &sys );
        dx = v * DEL;
        sys.StateIncrementX( xp, x, dx );
        dx = v * -DEL;
        sys.StateIncrementX( xm, x, dx );
        sys.StateScatter( xp, v, T, UpdateFlags::UPDATE_ALL_NO_VISUAL );
        sys.LoadConstraintJacobians();
        sys.GetConstraintJacobianMatrix( Cqp );
        sys.StateScatter( xm, v, T, UpdateFlags::UPDATE_ALL_NO_VISUAL );
        sys.LoadConstraintJacobians();
        sys.GetConstraintJacobianMatrix( Cqm );
        Qc = -( Cqp * v - Cqm * v ) / ( 2.0 * DEL );
        // Restore the true state before reading F, M, Cq.
        sys.StateScatter( x, v, T, UpdateFlags::UPDATE_ALL_NO_VISUAL );
        sys.LoadConstraintJacobians();
    }

    ChVectorDynamic<> F( n );
    F.setZero();
    sys.LoadResidual_F( F, 1.0 );

    ChSparseMatrix M( n, n );
    sys.GetMassMatrix( M );

    Eigen::MatrixXd K = Eigen::MatrixXd::Zero( n + m, n + m );
    K.topLeftCorner( n, n ) = Eigen::MatrixXd( M );
    if ( m > 0 ) {
        ChSparseMatrix Cq( m, n );
        sys.GetConstraintJacobianMatrix( Cq );
        Cqd = Eigen::MatrixXd( Cq );
        K.bottomLeftCorner( m, n ) = Cqd;
        K.topRightCorner( n, m ) = Cqd.transpose();
    }
    Eigen::VectorXd rhs( n + m );
    rhs.head( n ) = F;
    rhs.tail( m ) = Qc;

    Eigen::FullPivLU<Eigen::MatrixXd> lu( K );
    if ( !lu.isInvertible() ) return false;
    Eigen::VectorXd sol = lu.solve( rhs );

    ChStateDelta a( n, &sys );
    a = sol.head( n );
    sys.StateScatterAcceleration( a );
    return true;
}

// Extracts the .mech generalized accelerations from the Chrono bodies after
// an acceleration-level assembly, by projecting relative body accelerations
// back onto each joint's motion subspace.
static void extractQdd( const Model *mdl, ChronoModel &cm,
                        const double *qd, const std::vector<BodyKin> &kin, double *qdd ) {
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];
        ChronoBody &cb = cm.bod[ i ];
        int vi = mdl->vidx[ i ];

        // World kinematics of this body's mech-frame origin point.
        ChVector3d aI = cb.body->PointAccelerationLocalToParent( -cb.cg );
        ChVector3d wI = cb.body->GetAngVelParent();
        ChVector3d alI = cb.body->GetAngAccParent();
        ChMatrix33<> RI = cb.body->GetRotMat();

        // Same for the parent's origin point (zero for the fixed base).
        ChVector3d aP( 0, 0, 0 ), wP( 0, 0, 0 ), alP( 0, 0, 0 ), rP( 0, 0, 0 );
        ChMatrix33<> RP( 1 );
        if ( b->parent != -1 ) {
            ChronoBody &pb = cm.bod[ b->parent ];
            aP  = pb.body->PointAccelerationLocalToParent( -pb.cg );
            wP  = pb.body->GetAngVelParent();
            alP = pb.body->GetAngAccParent();
            rP  = pb.body->TransformPointLocalToParent( -pb.cg );
            RP  = pb.body->GetRotMat();
        }

        if ( b->jtype == JOINT_FLOATING ) {
            // qdd = d/dt of the body-frame spatial velocity components:
            //   angular: R^T * alpha
            //   linear:  R^T * a_origin - w_body x v_origin_body
            ChVector3d wB  = RI.transpose() * wI;
            ChVector3d vOB( qd[ vi + 3 ], qd[ vi + 4 ], qd[ vi + 5 ] );
            ChVector3d dw = RI.transpose() * alI;
            ChVector3d dv = RI.transpose() * aI - Vcross( wB, vOB );
            qdd[ vi + 0 ] = dw.x(); qdd[ vi + 1 ] = dw.y(); qdd[ vi + 2 ] = dw.z();
            qdd[ vi + 3 ] = dv.x(); qdd[ vi + 4 ] = dv.y(); qdd[ vi + 5 ] = dv.z();
        } else if ( b->jtype == JOINT_REVOLUTE ) {
            // w_rel = R_i a * qd  =>  qdd = a . R_i^T (al_i - al_p
            // - w_i x w_rel)   (qd/qdd live in the S = +a basis).
            ChVector3d wRel = wI - wP;
            ChVector3d v = RI.transpose() * ( alI - alP - Vcross( wI, wRel ) );
            qdd[ vi ] = Vdot( axis3( b->axis ), v );
        } else { // JOINT_PRISMATIC
            // r_i = r_p + u, u = R_p (tree + a q):
            //   a_i = a_p + al_p x u + w_p x (w_p x u) + 2 w_p x sdot + a_w qdd
            ChVector3d rI = cb.body->TransformPointLocalToParent( -cb.cg );
            ChVector3d u = rI - rP;
            ChVector3d aW = RP * axis3( b->axis );
            ChVector3d sdot = aW * qd[ vi ];
            ChVector3d resid = aI - aP - Vcross( alP, u ) - Vcross( wP, Vcross( wP, u ) )
                             - Vcross( wP, sdot ) * 2.0;
            qdd[ vi ] = Vdot( aW, resid );
        }
    }
}

// ---------------------------------------------------------------------------
// Cross-check harness, mirroring src/validate/main.c: staged-then-random
// states, scaled acceleration error against both the dense RBDA reference
// and the generated code.
// ---------------------------------------------------------------------------
static int crossCheck( const Model *mdl, const char *path ) {
    // Coordinate conventions first: if qd is not dq/dt, the acceleration
    // comparison below is meaningless no matter how well it agrees.
    if ( checkQdConsistency( mdl, path ) ) return 1;

    ChronoModel cm;
    buildChronoModel( mdl, cm );

    int nb = mdl->nbod, nq = mdl->nq, nv = mdl->nv;
    const int NTRIALS = 200;
    const double TOL = 1e-6;
    double maxErr = 0, maxErrGen = 0;
    int worstTrial = -1;
    std::vector<BodyKin> kin;

    for ( int t = 0 ; t < NTRIALS ; t++ ) {
        double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ], tau[ MODEL_MAX_V ];
        double fext[ MODEL_MAX_BODIES * 6 ];

        // Staging: t0 at rest, then add random q, qd, tau, fext one at a time
        // so a convention mismatch points at the term that introduced it.
        for ( int k = 0 ; k < nq ; k++ ) q[k]   = ( t == 0 ) ? 0 : 3.0 * frand();
        for ( int k = 0 ; k < nv ; k++ ) qd[k]  = ( t <= 1 ) ? 0 : 2.0 * frand();
        for ( int k = 0 ; k < nv ; k++ ) tau[k] = ( t <= 2 ) ? 0 : 5.0 * frand();
        for ( int k = 0 ; k < nb * 6 ; k++ )
            fext[k] = ( t <= 3 ) ? 0 : 5.0 * frand();
        if ( t == 0 )
            for ( int i = 0 ; i < nb ; i++ )
                if ( mdl->body[i].jtype == JOINT_FLOATING ) q[ mdl->qidx[i] ] = 1.0;

        setChronoState( mdl, cm, q, qd, tau, fext, kin );
        if ( !solveAccelerations( cm.sys ) ) {
            fprintf( stderr, "chrono_check: singular KKT system at trial %d\n", t );
            return 1;
        }

        double qddCh[ MODEL_MAX_V ], qddRef[ MODEL_MAX_V ], qddGen[ MODEL_MAX_V ];
        extractQdd( mdl, cm, qd, kin, qddCh );
        fd_dense( mdl, q, qd, tau, fext, qddRef );
        fwd_dyn_crba( q, qd, tau, fext, qddGen );

        double trialErr = 0;
        for ( int k = 0 ; k < nv ; k++ ) {
            double e = fabs( qddCh[k] - qddRef[k] ) / fmax( 1.0, fabs( qddRef[k] ) );
            if ( e > trialErr ) trialErr = e;
            if ( e > maxErr ) { maxErr = e; worstTrial = t; }

            // Chrono against the generated code directly, so the shipped
            // artifact is validated externally rather than transitively.
            double eg = fabs( qddCh[k] - qddGen[k] ) / fmax( 1.0, fabs( qddGen[k] ) );
            if ( eg > maxErrGen ) { maxErrGen = eg; }
        }
        if ( t <= 5 && trialErr > TOL ) {
            printf( "  trial %d (stage %s): err %.3e\n", t,
                    t == 0 ? "rest" : t == 1 ? "q" : t == 2 ? "q,qd" : t == 3 ? "q,qd,tau" : "all", trialErr );
            for ( int k = 0 ; k < nv ; k++ )
                printf( "    qdd[%d]: chrono % .12f  dense % .12f\n", k, qddCh[k], qddRef[k] );
        }
    }

    printf( "chrono_check: %s: %d trials, max scaled |qdd_chrono - qdd_dense| = %.3e (worst at trial %d)\n",
            path, NTRIALS, maxErr, worstTrial );
    printf( "chrono_check: %s: %d trials, max scaled |qdd_chrono - qdd_generated| = %.3e\n",
            path, NTRIALS, maxErrGen );
    if ( maxErr < TOL && maxErrGen < TOL ) {
        printf( "chrono_check: PASS\n" );
        return 0;
    }
    printf( "chrono_check: FAIL (tolerance %.1e)\n", TOL );
    return 1;
}

// ---------------------------------------------------------------------------
// Simulation-loop comparison: integrate the GENERATED dynamics forward with
// the integration rule documented in src/model/model.h, time-step the Chrono
// model alongside it from the same initial state, and compare where the
// bodies physically end up.
//
// This tests something no single-shot comparison can: that q, qd and the
// integration rule fit together over time. The acceleration cross-check above
// evaluates f(q, qd) at a fixed state and never forms dq/dt, so it is blind
// to the q<->qd relationship -- which is exactly how the old opposite-handed
// convention survived. That bug would fail this test on the first few steps
// (q running backwards against its own qd) while passing the 1e-9
// acceleration comparison.
//
// The two sides use different first-order integrators (semi-implicit Euler
// here, Chrono's default EULER_IMPLICIT_LINEARIZED there), so they are not
// expected to agree exactly -- they agree to their common O(dt) truncation
// error. runSim is therefore run twice, at dt and dt/2: if the residual is
// integrator truncation it roughly halves, whereas a real discrepancy in the
// dynamics would not shrink with the step.
// ---------------------------------------------------------------------------

// A deterministic initial state that actually moves. Several models are
// static under gravity from rest (threelink's revolute axes are all parallel
// to gravity, so no gravity torque reaches any joint), hence the nonzero
// velocity -- which also puts the Coriolis and centrifugal terms into the
// trajectory rather than just the gravity terms.
static void initialState( const Model *mdl, double *q, double *qd ) {
    for ( int k = 0 ; k < mdl->nq ; k++ ) q[k]  = 0.4 * sin( 1.0 + 2.0 * k );
    for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] = 0.8 * cos( 2.0 + 3.0 * k );

    // A floating body's q slots are a quaternion plus a position, not angles;
    // give it a definite tilt and offset instead of the sin() fill above.
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        if ( mdl->body[i].jtype != JOINT_FLOATING ) continue;
        double *qs = q + mdl->qidx[ i ];
        qs[0] = 0.9; qs[1] = 0.2; qs[2] = -0.3; qs[3] = 0.1; // normalized on use
        qs[4] = 0.0; qs[5] = 0.0; qs[6] = 0.0;
    }
}

// One semi-implicit (symplectic) Euler step of the .mech state, driven by the
// generated dynamics: qdd = f(q, qd); qd += qdd*dt; then advance q with the
// new qd via the documented rule. Free motion -- no joint torques, no
// external forces, gravity only.
static void stepGenerated( const Model *mdl, double *q, double *qd, double dt ) {
    double tau[ MODEL_MAX_V ] = { 0 }, fext[ MODEL_MAX_BODIES * 6 ] = { 0 };
    double qdd[ MODEL_MAX_V ];
    double qNext[ MODEL_MAX_Q ];

    fwd_dyn_crba( q, qd, tau, fext, qdd );
    for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] += qdd[k] * dt;
    integrateQ( mdl, q, qd, dt, qNext );
    for ( int k = 0 ; k < mdl->nq ; k++ ) q[k] = qNext[k];
}

// Runs both simulations to time T at step dt; returns the largest distance
// between any body's generated-side and Chrono-side position over the run.
static double runSim( const Model *mdl, double T, double dt, int printRows ) {
    ChronoModel cm;
    buildChronoModel( mdl, cm );

    double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ];
    initialState( mdl, q, qd );

    // Place Chrono at the same initial state. tau and fext are zero, so the
    // force accumulators stay empty for the whole run and gravity comes from
    // the Chrono system itself.
    double tau0[ MODEL_MAX_V ] = { 0 }, fext0[ MODEL_MAX_BODIES * 6 ] = { 0 };
    std::vector<BodyKin> kin;
    setChronoState( mdl, cm, q, qd, tau0, fext0, kin );

    int nsteps = (int)( T / dt + 0.5 );
    int rowEvery = printRows > 0 ? nsteps / printRows : 0;
    double maxDiv = 0;

    for ( int n = 1 ; n <= nsteps ; n++ ) {
        stepGenerated( mdl, q, qd, dt );
        cm.sys.DoStepDynamics( dt );

        fk( mdl, q, qd, kin );
        double stepDiv = 0;
        for ( int i = 0 ; i < mdl->nbod ; i++ ) {
            ChVector3d pCh = cm.bod[i].body->TransformPointLocalToParent( -cm.bod[i].cg );
            double d = ( kin[i].r - pCh ).Length();
            if ( d > stepDiv ) stepDiv = d;
        }
        if ( stepDiv > maxDiv ) maxDiv = stepDiv;

        if ( rowEvery && n % rowEvery == 0 ) {
            // Track the last body: for a serial chain body 0's origin sits on
            // its own joint and never translates, so the tip is what shows
            // the motion. Also print the worst gap across all bodies.
            int tip = mdl->nbod - 1;
            ChVector3d pCh = cm.bod[tip].body->TransformPointLocalToParent( -cm.bod[tip].cg );
            printf( "    t=%6.3f  gen b%d(% .4f % .4f % .4f)  chrono b%d(% .4f % .4f % .4f)  max|dp|=%.2e\n",
                    n * dt, tip, kin[tip].r.x(), kin[tip].r.y(), kin[tip].r.z(),
                    tip, pCh.x(), pCh.y(), pCh.z(), stepDiv );
        }
    }
    return maxDiv;
}

static int simCompare( const Model *mdl, const char *path, double T, double dt ) {
    printf( "chrono_check: %s: simulation loop, generated dynamics vs Chrono, T=%g dt=%g\n",
            path, T, dt );
    double divCoarse = runSim( mdl, T, dt, 8 );
    double divFine   = runSim( mdl, T, dt / 2, 0 );

    // Both integrators are first order, so halving the step should shrink the
    // gap: a ratio of about 2 (or better) means the two dynamics agree and
    // only the integrators differ, while a ratio near 1 at a gap well above
    // roundoff would mean a real discrepancy that refining the step cannot
    // remove. Some models never get that far: `gantry` is a Cartesian PPP
    // robot whose dynamics are linear, so its two trajectories agree to
    // roundoff and there is no truncation to converge.
    const double TOL = 1e-2;
    const double ROUNDOFF = 1e-11;
    double ratio = divFine > 0 ? divCoarse / divFine : 0;

    printf( "chrono_check: max |position gap| = %.3e at dt, %.3e at dt/2", divCoarse, divFine );
    if ( divCoarse <= ROUNDOFF )
        printf( "  (already at roundoff)\n" );
    else if ( ratio >= 1.5 )
        printf( "  (ratio %.2f -- shrinks with the step, so the gap is integrator truncation)\n", ratio );
    else
        printf( "  (ratio %.2f -- does NOT shrink with the step, so this is a real"
                " discrepancy in the dynamics, not truncation)\n", ratio );

    if ( divCoarse < TOL ) {
        printf( "chrono_check: SIM PASS (tolerance %.1e)\n", TOL );
        return 0;
    }
    printf( "chrono_check: SIM FAIL (tolerance %.1e)\n", TOL );
    return 1;
}

int main( int argc, char **argv ) {
    if ( argc < 2 ) {
        fprintf( stderr, "usage: chrono_check <model.mech> [--sim <T> <dt>]\n" );
        return 1;
    }
    Model mdl;
    char err[ 256 ];
    if ( modelLoad( argv[ 1 ], &mdl, err, sizeof err ) ) {
        fprintf( stderr, "chrono_check: %s\n", err );
        return 1;
    }
    const char *unsup = modelUnsupported( &mdl );
    if ( unsup ) {
        fprintf( stderr, "chrono_check: %s\n", unsup );
        return 1;
    }
    // This binary links one model's generated dynamics; catch a mismatch
    // rather than letting it read past the end of q.
    if ( fwd_dyn_crba_nq != mdl.nq || fwd_dyn_crba_nv != mdl.nv ) {
        fprintf( stderr, "chrono_check: linked dynamics.c is for nq=%d nv=%d but '%s' has nq=%d nv=%d -- regenerate\n",
                 fwd_dyn_crba_nq, fwd_dyn_crba_nv, argv[ 1 ], mdl.nq, mdl.nv );
        return 1;
    }

    if ( argc >= 5 && !strcmp( argv[ 2 ], "--sim" ) )
        return simCompare( &mdl, argv[ 1 ], atof( argv[ 3 ] ), atof( argv[ 4 ] ) );
    if ( crossCheck( &mdl, argv[ 1 ] ) ) return 1;
    return simCompare( &mdl, argv[ 1 ], 0.25, 2e-5 );
}
