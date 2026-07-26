// Mechanism definition file: parser and in-memory model.
//
// Shared by the code generator (src/codegen), the dense validator
// (src/validate), and the benchmark (src/bench), so there is exactly one
// definition of the mechanism and one parser.
//
// File format: line-based, '#' starts a comment, blank lines ignored.
//
//   gravity <gx> <gy> <gz>       physical gravity vector (e.g. 0 0 -9.81)
//
//   body                         starts body i (bodies are numbered in file order)
//   parent <p>                   parent body index; -1 = fixed base. Must be < i.
//   joint <type> [axis]          type: revolute | prismatic | cylindrical |
//                                spherical | floating
//                                axis: x | y | z (required for revolute,
//                                prismatic, cylindrical)
//
// Joint coordinate conventions
// ----------------------------
// A joint coordinate q is the ordinary physical coordinate, and qd is its
// honest time derivative:
//
//   revolute:  q is a right-handed rotation of the child frame by +q about
//              the +axis, qd = dq/dt, and tau is a couple about the +axis.
//   prismatic: q is the displacement along the +axis, qd = dq/dt, and tau is
//              a force along the +axis.
//
// So a single-DOF joint integrates with plain q += qd*dt. This is worth
// stating because it is NOT what falls out of Featherstone's formulation for
// free: the Plucker E block is a coordinate transform (child-from-parent), so
// a child physically rotated by R has E = R^T. spatial_v2 folds that
// transpose into rotx/roty/rotz; this project instead writes the ordinary R
// and transposes it explicitly where the transform is built (see
// Xform::setJointRot in src/codegen/xform.h and jointXS in
// src/validate/dense.h). Using R directly as E -- which this project did
// until the Chrono cross-check caught it -- silently makes q and qd
// opposite-handed: forward dynamics stays self-consistent, but any
// integrator stepping q += qd*dt drives the mechanism backwards.
// src/chrono/main.cpp asserts qd == dq/dt on every model to keep this fixed.
//
// A 'floating' joint gives the body all six DOF relative to the fixed base:
// a genuine free-flying body (RBDA sec 9.4), with the orientation carried as
// a unit quaternion so there is no gimbal lock anywhere. It must be attached
// to the base (parent -1). Because a quaternion has four components, the
// position and velocity coordinate counts differ (nq != nv):
//
//   q  slots (7): [ qw qx qy qz ] unit quaternion (normalized on input, so
//                 any nonzero quaternion is accepted) giving the body's
//                 ordinary attitude, i.e. it rotates BODY coordinates into
//                 BASE coordinates; then [ x y z ], the position of the body
//                 frame in base coordinates
//   qd slots (6): the body's spatial velocity IN ITS OWN FRAME,
//                 [ wx wy wz vx vy vz ]
//
// out_qdd gets the corresponding 6 spatial-acceleration slots. Because q and
// qd are in different frames here, a floating body does NOT integrate with
// q += qd*dt; the rates are
//
//   Qdot = 0.5 * Q x (0, w)      (quaternion product, w = qd[0..2], body frame)
//   rdot = R(Q) * v              (v = qd[3..5], body frame)
//
// which is the integrator's job, not this function's. integrateQ() in
// src/chrono/main.cpp is a worked reference implementation. Use
// modelNq()/modelNv() (m->nq, m->nv) and the per-body offsets m->qidx[i],
// m->vidx[i] to address the packed q/qd vectors.
//   tree <x> <y> <z>             translation from the parent's frame to this
//                                joint's frame (no fixed rotation yet)
//   mass <m>
//   cg <x> <y> <z>               CG location in this body's frame
//   inertia <xx> <yy> <zz> <xy> <yz> <zx>
//                                rotational inertia about the CG; actual matrix
//                                entries (off-diagonal minus signs absorbed)
//   fext <0|1>                   1 = generate support for an external spatial
//                                force on this body
//
// Kinematic loops (RBDA ch. 8)
// ----------------------------
// After the bodies, a mechanism may declare loop-closing joints. The bodies
// always form a spanning tree; a loop joint adds CONSTRAINTS between two
// existing bodies, never coordinates, so nq/nv are untouched (method 2 of
// RBDA sec 8.5: solve for the constraint forces lambda first).
//
//   stab <Tstab>                 global: Baumgarte stabilization time constant
//                                in seconds (RBDA sec 8.3, alpha = beta =
//                                1/Tstab). Default 0.1.
//
//   loop                         starts a loop joint
//   predecessor <p>              body index, or -1 for the fixed base
//   successor <s>                body index, or -1 for the fixed base
//   joint <type> [axis]          type: revolute | prismatic | cylindrical |
//                                spherical | point-on-line | distance
//                                axis: x | y | z (all types except spherical
//                                and distance)
//   ptree <x> <y> <z>            translation, predecessor body frame -> loop
//                                joint frame (RBDA's XP)
//   stree <x> <y> <z>            translation, successor body frame -> loop
//                                joint frame (RBDA's XS)
//   length <d>                   distance joints only: the imposed distance
//                                between the ptree and stree anchor points
//                                (must be positive)
//
// A loop joint imposes nc constraints: revolute/prismatic 5, cylindrical 4,
// spherical 3, point-on-line 2, distance 1. point-on-line is RBDA sec 8.10's
// sphere-in-cylinder: it constrains only the two translations perpendicular
// to the axis, and exists so planar loops can be closed WITHOUT redundant
// constraints -- a planar loop closed by a revolute has rank(K) = 2 < nc = 5,
// which this pipeline rejects (the validator checks rank(K) = nc and fails
// the build; see README). Close planar loops with point-on-line.
//
// distance is a rigid massless rod of length d between the ptree point on
// the predecessor and the stree point on the successor (ball joints at both
// ends): the single scalar constraint |p_s - p_p| = d. Its position error is
// carried in the squared form delta = (|rho|^2 - d^2)/2, rho = p_s - p_p,
// which is polynomial in the state (no square roots anywhere) and whose
// gradient is exactly the K row, so closure projection works unchanged.
//
// Points of interest and force elements
// --------------------------------------
// After the bodies (interleaving with loop blocks is fine; a 'poi', 'tsda',
// or 'rsda' keyword ends any open body/loop block, exactly as 'body' and
// 'loop' end each other):
//
//   poi <name> <body> <x> <y> <z>   one line; named point at (x y z) in
//                                   <body>'s frame; body -1 = fixed base
//
//   tsda <name>                     translational spring-damper-actuator
//     from <poi-name>               anchor 1 (order defines rho = p_to - p_from)
//     to <poi-name>                 anchor 2
//     k <k>                         spring rate, N/m       (default 0)
//     c <c>                         damping, N/(m/s)       (default 0)
//     l0 <l0>                       free length, m         (default 0)
//     input                         optional: reserve an act[] slot
//
//   rsda <name>                     rotational spring-damper-actuator
//     joint <body>                  acts on the tree joint of body <body>
//     k <k>                         Nm/rad (or N/m on a prismatic joint)
//     c <c>
//     q0 <q0>                       rest coordinate        (default 0)
//     input
//
// Force law (part of the contract):
//
//   TSDA: scalar tension f = k*(l - l0) + c*ldot + act[slot], positive f
//         pulls the two anchors together. l = |rho|, rho = p_to - p_from,
//         ldot = (rho . rhodot) / l.
//   RSDA: generalized force on the joint coordinate
//         tau_e = -k*(q - q0) - c*qd + act[slot], positive act drives +q.
//
// A TSDA anchor pair coincident on the same body at the identical local
// point is a model error: l = 0 for all time and f/l divides by zero.
//
// Validation: POI names must be unique; TSDA names unique among TSDAs and
// RSDA names unique among RSDAs (names are length-limited to
// MODEL_NAME_LEN - 1 and counts to MODEL_MAX_POI / MODEL_MAX_ELEM). A poi's
// body must be -1 or a declared body index. A tsda's anchors must resolve
// to distinct, already-declared POIs. An rsda's joint must be a single-DOF
// tree joint (revolute or prismatic) -- floating/spherical/cylindrical are
// rejected. act[] slots are assigned in declaration order across both
// element kinds.
//
// All joint types parse; which ones a given consumer implements is up to it
// (see modelUnsupported below for the current single-DOF pipeline).

#pragma once
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Bumped from 32 so a generated track band (tens of shoes per side, two sides)
// fits. A ModelBody is small and Model is heap/stack-friendly at this size.
#define MODEL_MAX_BODIES 256
#define MODEL_MAX_Q  ( MODEL_MAX_BODIES * 7 ) // worst case: every body floating
#define MODEL_MAX_V  ( MODEL_MAX_BODIES * 6 )

#define MODEL_NAME_LEN 32
#define MODEL_MAX_POI 24
#define MODEL_MAX_ELEM 24

#define MODEL_MAX_TRACKS  4    // track bands (typically 2: left + right)
#define MODEL_MAX_ROLLERS 24   // rollers/wheels a single band may ride on

typedef enum {
    JOINT_REVOLUTE,
    JOINT_PRISMATIC,
    JOINT_CYLINDRICAL,
    JOINT_SPHERICAL,
    JOINT_FLOATING,    // genuine 6-DOF joint, quaternion orientation
    JOINT_POINTONLINE, // sphere-in-cylinder; loop joints only
    JOINT_DISTANCE,    // rigid rod between two anchor points; loop joints only
} ModelJointType;

typedef struct {
    int parent;            // -1 = base
    ModelJointType jtype;
    char axis;             // 'x' | 'y' | 'z'; 0 for spherical
    double tree[3];        // translation, parent frame -> joint frame
    double mass;
    double cg[3];
    double inertia[6];     // xx yy zz xy yz zx, about the CG
    int hasFext;
} ModelBody;

#define MODEL_MAX_LOOPS 40 // a 4-corner car needs 16 (4/corner); track bands add 1 closure each
#define MODEL_MAX_NC ( MODEL_MAX_LOOPS * 5 ) // worst case: every loop revolute/prismatic

typedef struct {
    int pred;              // predecessor body index; -1 = base
    int succ;              // successor body index; -1 = base
    ModelJointType jtype;  // constraint pattern only; a loop joint has no q/qd
    char axis;             // 'x' | 'y' | 'z'; 0 for spherical and distance
    double ptree[3];       // translation, predecessor body frame -> loop joint frame (XP)
    double stree[3];       // translation, successor body frame -> loop joint frame (XS)
    double length;         // distance joints: the imposed anchor-to-anchor distance
} ModelLoop;

// A named point fixed to a body's frame ("body" -1 = fixed base). A TSDA
// anchor *is* a poi, so anchors and observation points share one
// declaration; see the header comment above for the grammar.
typedef struct {
    char name[ MODEL_NAME_LEN ];
    int body;               // -1 = fixed base
    double p[3];             // local coordinates in body's frame
    int obsIdx;              // offset into out_obs (6 doubles: px py pz vx vy vz)
} ModelPoi;

typedef enum { ELEM_TSDA, ELEM_RSDA } ModelElemType;

// A force element: either a TSDA between two POIs or an RSDA on a tree
// joint. See the header comment above for the force-law signs.
typedef struct {
    char name[ MODEL_NAME_LEN ];
    ModelElemType type;
    int poi1, poi2;          // TSDA: 'from'/'to' indices into Model::poi[]
    int joint;                // RSDA: body index owning the tree joint
    double k, c, ref;         // ref = l0 (TSDA) or q0 (RSDA)
    int hasInput;
    int actIdx;                // slot in act[], -1 if !hasInput
    int obsIdx;                 // offset into out_obs (2 doubles)
} ModelElem;

// A cylinder the band rides on / contacts. 'isEnd' marks the two path-defining
// end pulleys (sprocket + idler) whose centers set the belt geometry; 'isDrive'
// marks the sprocket's positive tooth engagement. Body centers are read at the
// design pose from the referenced body's tree offsets (see trkBodyPosInRoot).
typedef struct {
    int body;                // roller/wheel body index (declared before the track)
    double radius;
    int isEnd;               // 1 = path-defining end pulley
    int isDrive;             // 1 = sprocket, positive tooth engagement
    int teeth;               // tooth count when isDrive
} ModelTrackRoller;

// A 'track' block: a closed planar band of identical shoes riding on the listed
// rollers. modelLoad expands it into shoe BODIES (a floating anchor shoe + a
// chain of revolute-y pins) plus one point-on-line closure loop; the contact
// data below is parsed and kept for the (separate) contact force engine, which
// this loader does not implement. The band lies in the root body's x-z plane at
// the end pulleys' lateral offset; every pin spins about y.
typedef struct {
    char name[ MODEL_NAME_LEN ];
    int root;                // hull body: geometry reference frame for the band
    int shoes;               // N
    double pitch;            // declared pin-to-pin (advisory; see linkLen)
    double mass;             // per shoe
    double cg[3];
    double inertia[6];
    int hasFext;
    double kn, kt, mu;       // contact material (for the force engine)
    int nroller;
    ModelTrackRoller roller[ MODEL_MAX_ROLLERS ];
    int hasGround; char groundAxis; double groundOffset;
    // filled during expansion:
    int shoe0;               // body index of the floating anchor shoe
    double linkLen;          // actual shoe length used = perimeter / N
} ModelTrack;

typedef struct {
    int nbod;
    int nq;                // total position coordinates (7 per floating body, 1 otherwise)
    int nv;                // total velocity coordinates (6 per floating body, 1 otherwise)
    int qidx[ MODEL_MAX_BODIES ]; // body i's first slot in q
    int vidx[ MODEL_MAX_BODIES ]; // body i's first slot in qd / tau / out_qdd
    double gravity[3];     // physical gravity
    ModelBody body[ MODEL_MAX_BODIES ];
    int nloop;
    int nc;                // total loop constraints, sum of modelLoopNc
    int cidx[ MODEL_MAX_LOOPS ];  // loop l's first constraint row in K / lambda
    double stab;           // Baumgarte time constant Tstab (s); alpha = beta = 1/stab
    ModelLoop loop[ MODEL_MAX_LOOPS ];
    int npoi;
    ModelPoi poi[ MODEL_MAX_POI ];
    int nelem;
    ModelElem elem[ MODEL_MAX_ELEM ];
    int nact;               // total act[] slots, sum of elem[].hasInput
    int obsLen;              // fwd_dyn_observe_len = 6*npoi + 2*nelem
    int ntrack;
    ModelTrack track[ MODEL_MAX_TRACKS ];
    double q0[ MODEL_MAX_Q ]; // design-pose coordinates: 0 / quaternion-identity
                              // for hand-written bodies, the curved band pose for
                              // generated shoes. Consumers seed initial q from this.
} Model;

static int modelBodyNq( const ModelBody *b ) { return b->jtype == JOINT_FLOATING ? 7 : 1; }
static int modelBodyNv( const ModelBody *b ) { return b->jtype == JOINT_FLOATING ? 6 : 1; }

// Constraints imposed by a loop joint (6 minus its motion freedoms).
static int modelLoopNc( const ModelLoop *l ) {
    switch ( l->jtype ) {
        case JOINT_REVOLUTE:
        case JOINT_PRISMATIC:   return 5;
        case JOINT_CYLINDRICAL: return 4;
        case JOINT_SPHERICAL:   return 3;
        case JOINT_POINTONLINE: return 2;
        case JOINT_DISTANCE:    return 1;
        default:                return 0; // unreachable; floating rejected at parse
    }
}

// Index of a poi by name, or -1 if not found.
static int modelFindPoi( const Model *m, const char *name ) {
    for ( int i = 0 ; i < m->npoi ; i++ )
        if ( !strcmp( m->poi[ i ].name, name ) ) return i;
    return -1;
}

// ---------------------------------------------------------------------------
// Track band geometry: an open (external-tangent) belt around two end pulleys.
// The band is the two straight runs tangent to both circles plus the two outer
// wrap arcs; shoes are laid evenly by arc length. Everything is 2D in the band
// plane, coordinates (x,z) in the root frame. See models notes for the design.
// ---------------------------------------------------------------------------
typedef struct { double x, z; } TrkV2;
typedef struct {
    int arc;                     // 1 = circular arc, 0 = straight run
    TrkV2 p0, dir;               // straight: p(t) = p0 + t*dir (dir unit)
    TrkV2 c; double r, a0, dw;   // arc: center, radius, start angle, signed sweep
    double len;                  // segment arc length
} TrkSeg;

// Build the belt around pulleys A and B (band on the OUTSIDE of both); fill
// seg[0..3] in circulation order (bottom run, B arc, top run, A arc); return
// the perimeter. Caller guarantees |rA-rB| < |B-A|.
static double trkBuildBelt( TrkV2 cA, double rA, TrkV2 cB, double rB, TrkSeg seg[4] ) {
    TrkV2 d; d.x = cB.x - cA.x; d.z = cB.z - cA.z;
    double D = sqrt( d.x*d.x + d.z*d.z );
    TrkV2 u; u.x = d.x / D; u.z = d.z / D;
    TrkV2 perp; perp.x = -u.z; perp.z = u.x;         // +90 deg in (x,z)
    double alpha = asin( ( rA - rB ) / D ), ca = cos( alpha ), sa = sin( alpha );
    TrkV2 mup, mlo;                                   // tangent-line outward normals
    mup.x = sa*u.x + ca*perp.x; mup.z = sa*u.z + ca*perp.z;
    mlo.x = sa*u.x - ca*perp.x; mlo.z = sa*u.z - ca*perp.z;
    TrkV2 Aup, Bup, Alo, Blo;
    Aup.x = cA.x + rA*mup.x; Aup.z = cA.z + rA*mup.z;
    Bup.x = cB.x + rB*mup.x; Bup.z = cB.z + rB*mup.z;
    Alo.x = cA.x + rA*mlo.x; Alo.z = cA.z + rA*mlo.z;
    Blo.x = cB.x + rB*mlo.x; Blo.z = cB.z + rB*mlo.z;
    double Lt = sqrt( D*D - ( rA - rB )*( rA - rB ) );

    seg[0].arc = 0; seg[0].p0 = Alo;                  // bottom run Alo -> Blo
    seg[0].dir.x = ( Blo.x - Alo.x )/Lt; seg[0].dir.z = ( Blo.z - Alo.z )/Lt; seg[0].len = Lt;
    seg[1].arc = 1; seg[1].c = cB; seg[1].r = rB;     // B arc, CCW, wrap pi-2a
    seg[1].a0 = atan2( mlo.z, mlo.x ); seg[1].dw = M_PI - 2*alpha; seg[1].len = rB*seg[1].dw;
    seg[2].arc = 0; seg[2].p0 = Bup;                  // top run Bup -> Aup
    seg[2].dir.x = ( Aup.x - Bup.x )/Lt; seg[2].dir.z = ( Aup.z - Bup.z )/Lt; seg[2].len = Lt;
    seg[3].arc = 1; seg[3].c = cA; seg[3].r = rA;     // A arc, CCW, wrap pi+2a
    seg[3].a0 = atan2( mup.z, mup.x ); seg[3].dw = M_PI + 2*alpha; seg[3].len = rA*seg[3].dw;

    return seg[0].len + seg[1].len + seg[2].len + seg[3].len;
}

// Point at arc length s along the closed 4-segment path; also its tangent angle.
static TrkV2 trkPoint( const TrkSeg seg[4], double s, double *phi ) {
    for ( int k = 0 ; k < 4 ; k++ ) {
        if ( s < seg[k].len || k == 3 ) {
            TrkV2 p;
            if ( !seg[k].arc ) {
                p.x = seg[k].p0.x + s*seg[k].dir.x;
                p.z = seg[k].p0.z + s*seg[k].dir.z;
                *phi = atan2( seg[k].dir.z, seg[k].dir.x );
            } else {
                double sgn = seg[k].dw >= 0 ? 1.0 : -1.0;
                double a = seg[k].a0 + sgn*( s / seg[k].r );
                p.x = seg[k].c.x + seg[k].r*cos( a );
                p.z = seg[k].c.z + seg[k].r*sin( a );
                *phi = atan2( sgn*cos( a ), -sgn*sin( a ) );
            }
            return p;
        }
        s -= seg[k].len;
    }
    TrkV2 z; z.x = 0; z.z = 0; *phi = 0; return z;    // unreachable
}

// Design-pose position of body b in body 'root's frame. Valid because every
// joint between them is at q=0 (identity rotation) in the design pose, so the
// position is just the sum of tree translations. Returns 0 on success, nonzero
// if root is not an ancestor or a floating body sits in between.
static int trkBodyPosInRoot( const Model *m, int b, int root, double out[3] ) {
    out[0] = out[1] = out[2] = 0;
    int guard = 0;
    while ( b != root ) {
        if ( b < 0 ) return -1;                        // hit base without root
        const ModelBody *bb = &m->body[ b ];
        if ( bb->jtype == JOINT_FLOATING ) return -2;  // breaks the q=0 sum
        out[0] += bb->tree[0]; out[1] += bb->tree[1]; out[2] += bb->tree[2];
        b = bb->parent;
        if ( ++guard > MODEL_MAX_BODIES ) return -1;
    }
    return 0;
}

// Loads a mechanism file. Returns 0 on success; on failure returns nonzero
// and writes a message (with file/line context) into err.
static int modelLoad( const char *path, Model *m, char *err, size_t errSize ) {
    FILE *f = fopen( path, "r" );
    if ( !f ) {
        snprintf( err, errSize, "%s: cannot open file", path );
        return 1;
    }

    memset( m, 0, sizeof *m );
    m->stab = 0.1;
    ModelBody *cur = NULL;
    ModelLoop *curLoop = NULL;
    ModelElem *curElem = NULL;
    ModelTrack *curTrack = NULL;
    char line[ 512 ];
    int lineNo = 0;

    while ( fgets( line, sizeof line, f ) ) {
        lineNo++;
        char *hash = strchr( line, '#' );
        if ( hash ) *hash = 0;

        char key[ 64 ];
        int pos = 0;
        if ( sscanf( line, " %63s%n", key, &pos ) != 1 ) continue; // blank line
        const char *rest = line + pos;

        if ( !strcmp( key, "gravity" ) ) {
            if ( sscanf( rest, "%lf %lf %lf", &m->gravity[0], &m->gravity[1], &m->gravity[2] ) != 3 ) goto badArgs;
            continue;
        }
        if ( !strcmp( key, "stab" ) ) {
            if ( sscanf( rest, "%lf", &m->stab ) != 1 ) goto badArgs;
            continue;
        }
        if ( !strcmp( key, "body" ) ) {
            if ( m->nbod >= MODEL_MAX_BODIES ) {
                snprintf( err, errSize, "%s:%d: too many bodies (max %d)", path, lineNo, MODEL_MAX_BODIES );
                goto fail;
            }
            cur = &m->body[ m->nbod++ ];
            curLoop = NULL;
            curElem = NULL;
            curTrack = NULL;
            memset( cur, 0, sizeof *cur );
            cur->parent = -1;
            cur->jtype = JOINT_REVOLUTE;
            cur->axis = 'z';
            continue;
        }
        if ( !strcmp( key, "loop" ) ) {
            if ( m->nloop >= MODEL_MAX_LOOPS ) {
                snprintf( err, errSize, "%s:%d: too many loops (max %d)", path, lineNo, MODEL_MAX_LOOPS );
                goto fail;
            }
            curLoop = &m->loop[ m->nloop++ ];
            cur = NULL;
            curElem = NULL;
            curTrack = NULL;
            memset( curLoop, 0, sizeof *curLoop );
            curLoop->pred = -2; // "unset" so validation can insist both ends are given
            curLoop->succ = -2;
            curLoop->jtype = JOINT_REVOLUTE;
            curLoop->axis = 'z';
            continue;
        }
        if ( !strcmp( key, "poi" ) ) {
            cur = NULL;
            curLoop = NULL;
            curElem = NULL;
            curTrack = NULL;
            if ( m->npoi >= MODEL_MAX_POI ) {
                snprintf( err, errSize, "%s:%d: too many pois (max %d)", path, lineNo, MODEL_MAX_POI );
                goto fail;
            }
            char name[ 64 ] = "";
            int body;
            double x, y, z;
            if ( sscanf( rest, "%63s %d %lf %lf %lf", name, &body, &x, &y, &z ) != 5 ) goto badArgs;
            if ( strlen( name ) >= MODEL_NAME_LEN ) {
                snprintf( err, errSize, "%s:%d: poi name '%s' too long (max %d chars)", path, lineNo, name, MODEL_NAME_LEN - 1 );
                goto fail;
            }
            ModelPoi *poi = &m->poi[ m->npoi++ ];
            memset( poi, 0, sizeof *poi );
            strcpy( poi->name, name );
            poi->body = body;
            poi->p[0] = x; poi->p[1] = y; poi->p[2] = z;
            continue;
        }
        if ( !strcmp( key, "tsda" ) || !strcmp( key, "rsda" ) ) {
            cur = NULL;
            curLoop = NULL;
            curTrack = NULL;
            if ( m->nelem >= MODEL_MAX_ELEM ) {
                snprintf( err, errSize, "%s:%d: too many elements (max %d)", path, lineNo, MODEL_MAX_ELEM );
                goto fail;
            }
            char name[ 64 ] = "";
            if ( sscanf( rest, "%63s", name ) != 1 ) goto badArgs;
            if ( strlen( name ) >= MODEL_NAME_LEN ) {
                snprintf( err, errSize, "%s:%d: element name '%s' too long (max %d chars)", path, lineNo, name, MODEL_NAME_LEN - 1 );
                goto fail;
            }
            curElem = &m->elem[ m->nelem++ ];
            memset( curElem, 0, sizeof *curElem );
            strcpy( curElem->name, name );
            curElem->type = !strcmp( key, "tsda" ) ? ELEM_TSDA : ELEM_RSDA;
            curElem->poi1 = curElem->poi2 = -1;
            curElem->joint = -1;
            curElem->actIdx = -1;
            continue;
        }
        if ( !strcmp( key, "track" ) ) {
            cur = NULL;
            curLoop = NULL;
            curElem = NULL;
            if ( m->ntrack >= MODEL_MAX_TRACKS ) {
                snprintf( err, errSize, "%s:%d: too many tracks (max %d)", path, lineNo, MODEL_MAX_TRACKS );
                goto fail;
            }
            curTrack = &m->track[ m->ntrack++ ];
            memset( curTrack, 0, sizeof *curTrack );
            curTrack->root = -1;
            curTrack->mu = 0.8;             // sensible contact defaults
            curTrack->shoe0 = -1;
            char name[ 64 ] = "";
            if ( sscanf( rest, "%63s", name ) == 1 ) {
                if ( strlen( name ) >= MODEL_NAME_LEN ) {
                    snprintf( err, errSize, "%s:%d: track name '%s' too long (max %d chars)", path, lineNo, name, MODEL_NAME_LEN - 1 );
                    goto fail;
                }
                strcpy( curTrack->name, name );
            }
            continue;
        }

        // Track-block keys.
        if ( curTrack ) {
            if ( !strcmp( key, "root" ) ) {
                if ( sscanf( rest, "%d", &curTrack->root ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "plane" ) ) {
                char pl[ 8 ] = "";
                if ( sscanf( rest, "%7s", pl ) != 1 ) goto badArgs;
                if ( strcmp( pl, "xz" ) ) {
                    snprintf( err, errSize, "%s:%d: track: only 'plane xz' is supported (pin axis y)", path, lineNo );
                    goto fail;
                }
            } else if ( !strcmp( key, "pin" ) ) {
                char type[ 16 ] = "", ax[ 8 ] = "";
                if ( sscanf( rest, "%15s %7s", type, ax ) != 2 || strcmp( type, "revolute" ) || strcmp( ax, "y" ) ) {
                    snprintf( err, errSize, "%s:%d: track: pin must be 'revolute y'", path, lineNo );
                    goto fail;
                }
            } else if ( !strcmp( key, "shoes" ) ) {
                if ( sscanf( rest, "%d", &curTrack->shoes ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "pitch" ) ) {
                if ( sscanf( rest, "%lf", &curTrack->pitch ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "mass" ) ) {
                if ( sscanf( rest, "%lf", &curTrack->mass ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "cg" ) ) {
                if ( sscanf( rest, "%lf %lf %lf", &curTrack->cg[0], &curTrack->cg[1], &curTrack->cg[2] ) != 3 ) goto badArgs;
            } else if ( !strcmp( key, "inertia" ) ) {
                if ( sscanf( rest, "%lf %lf %lf %lf %lf %lf",
                             &curTrack->inertia[0], &curTrack->inertia[1], &curTrack->inertia[2],
                             &curTrack->inertia[3], &curTrack->inertia[4], &curTrack->inertia[5] ) != 6 ) goto badArgs;
            } else if ( !strcmp( key, "fext" ) ) {
                if ( sscanf( rest, "%d", &curTrack->hasFext ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "kn" ) ) {
                if ( sscanf( rest, "%lf", &curTrack->kn ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "kt" ) ) {
                if ( sscanf( rest, "%lf", &curTrack->kt ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "mu" ) ) {
                if ( sscanf( rest, "%lf", &curTrack->mu ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "roller" ) ) {
                if ( curTrack->nroller >= MODEL_MAX_ROLLERS ) {
                    snprintf( err, errSize, "%s:%d: track '%s': too many rollers (max %d)", path, lineNo, curTrack->name, MODEL_MAX_ROLLERS );
                    goto fail;
                }
                ModelTrackRoller *rr = &curTrack->roller[ curTrack->nroller++ ];
                memset( rr, 0, sizeof *rr );
                if ( sscanf( rest, "%d %lf", &rr->body, &rr->radius ) != 2 ) goto badArgs;
                if ( strstr( rest, "end" ) )   rr->isEnd = 1;
                if ( strstr( rest, "drive" ) ) rr->isDrive = 1;
                const char *tp = strstr( rest, "teeth" );
                if ( tp ) sscanf( tp, "teeth %d", &rr->teeth );
            } else if ( !strcmp( key, "ground" ) ) {
                char ax[ 8 ] = "";
                if ( sscanf( rest, "%7s %lf", ax, &curTrack->groundOffset ) != 2 ) goto badArgs;
                curTrack->groundAxis = ax[0];
                curTrack->hasGround = 1;
            } else {
                snprintf( err, errSize, "%s:%d: unknown keyword '%s' in track block", path, lineNo, key );
                goto fail;
            }
            continue;
        }

        // Element-block keys (tsda/rsda).
        if ( curElem ) {
            if ( !strcmp( key, "from" ) || !strcmp( key, "to" ) ) {
                if ( curElem->type != ELEM_TSDA ) {
                    snprintf( err, errSize, "%s:%d: '%s' is only valid in a tsda block", path, lineNo, key );
                    goto fail;
                }
                char name[ 64 ];
                if ( sscanf( rest, "%63s", name ) != 1 ) goto badArgs;
                int idx = modelFindPoi( m, name );
                if ( idx < 0 ) {
                    snprintf( err, errSize, "%s:%d: tsda '%s': unknown poi '%s'", path, lineNo, curElem->name, name );
                    goto fail;
                }
                if ( !strcmp( key, "from" ) ) curElem->poi1 = idx;
                else                          curElem->poi2 = idx;
            } else if ( !strcmp( key, "joint" ) ) {
                if ( curElem->type != ELEM_RSDA ) {
                    snprintf( err, errSize, "%s:%d: 'joint' is only valid in an rsda block", path, lineNo );
                    goto fail;
                }
                if ( sscanf( rest, "%d", &curElem->joint ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "k" ) ) {
                if ( sscanf( rest, "%lf", &curElem->k ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "c" ) ) {
                if ( sscanf( rest, "%lf", &curElem->c ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "l0" ) ) {
                if ( curElem->type != ELEM_TSDA ) {
                    snprintf( err, errSize, "%s:%d: 'l0' is only valid in a tsda block", path, lineNo );
                    goto fail;
                }
                if ( sscanf( rest, "%lf", &curElem->ref ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "q0" ) ) {
                if ( curElem->type != ELEM_RSDA ) {
                    snprintf( err, errSize, "%s:%d: 'q0' is only valid in an rsda block", path, lineNo );
                    goto fail;
                }
                if ( sscanf( rest, "%lf", &curElem->ref ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "input" ) ) {
                curElem->hasInput = 1;
            } else {
                snprintf( err, errSize, "%s:%d: unknown keyword '%s' in %s block", path, lineNo, key,
                          curElem->type == ELEM_TSDA ? "tsda" : "rsda" );
                goto fail;
            }
            continue;
        }

        // Loop-block keys.
        if ( curLoop ) {
            if ( !strcmp( key, "predecessor" ) ) {
                if ( sscanf( rest, "%d", &curLoop->pred ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "successor" ) ) {
                if ( sscanf( rest, "%d", &curLoop->succ ) != 1 ) goto badArgs;
            } else if ( !strcmp( key, "joint" ) ) {
                char type[ 32 ] = "", axis[ 8 ] = "";
                int n = sscanf( rest, "%31s %7s", type, axis );
                if ( n < 1 ) goto badArgs;
                if      ( !strcmp( type, "revolute" ) )      curLoop->jtype = JOINT_REVOLUTE;
                else if ( !strcmp( type, "prismatic" ) )     curLoop->jtype = JOINT_PRISMATIC;
                else if ( !strcmp( type, "cylindrical" ) )   curLoop->jtype = JOINT_CYLINDRICAL;
                else if ( !strcmp( type, "spherical" ) )     curLoop->jtype = JOINT_SPHERICAL;
                else if ( !strcmp( type, "point-on-line" ) ) curLoop->jtype = JOINT_POINTONLINE;
                else if ( !strcmp( type, "distance" ) )      curLoop->jtype = JOINT_DISTANCE;
                else {
                    snprintf( err, errSize, "%s:%d: unknown loop joint type '%s'", path, lineNo, type );
                    goto fail;
                }
                if ( curLoop->jtype == JOINT_SPHERICAL || curLoop->jtype == JOINT_DISTANCE ) {
                    curLoop->axis = 0;
                } else {
                    if ( n < 2 || axis[1] != 0 || ( axis[0] != 'x' && axis[0] != 'y' && axis[0] != 'z' ) ) {
                        snprintf( err, errSize, "%s:%d: loop joint '%s' needs an axis of x, y, or z", path, lineNo, type );
                        goto fail;
                    }
                    curLoop->axis = axis[0];
                }
            } else if ( !strcmp( key, "ptree" ) ) {
                if ( sscanf( rest, "%lf %lf %lf", &curLoop->ptree[0], &curLoop->ptree[1], &curLoop->ptree[2] ) != 3 ) goto badArgs;
            } else if ( !strcmp( key, "stree" ) ) {
                if ( sscanf( rest, "%lf %lf %lf", &curLoop->stree[0], &curLoop->stree[1], &curLoop->stree[2] ) != 3 ) goto badArgs;
            } else if ( !strcmp( key, "length" ) ) {
                if ( sscanf( rest, "%lf", &curLoop->length ) != 1 ) goto badArgs;
            } else {
                snprintf( err, errSize, "%s:%d: unknown keyword '%s' in loop block", path, lineNo, key );
                goto fail;
            }
            continue;
        }

        // Everything below applies to the current body.
        if ( !cur ) {
            snprintf( err, errSize, "%s:%d: '%s' before any 'body'", path, lineNo, key );
            goto fail;
        }
        if ( !strcmp( key, "parent" ) ) {
            if ( sscanf( rest, "%d", &cur->parent ) != 1 ) goto badArgs;
        } else if ( !strcmp( key, "joint" ) ) {
            char type[ 32 ] = "", axis[ 8 ] = "";
            int n = sscanf( rest, "%31s %7s", type, axis );
            if ( n < 1 ) goto badArgs;
            if      ( !strcmp( type, "revolute" ) )    cur->jtype = JOINT_REVOLUTE;
            else if ( !strcmp( type, "prismatic" ) )   cur->jtype = JOINT_PRISMATIC;
            else if ( !strcmp( type, "cylindrical" ) ) cur->jtype = JOINT_CYLINDRICAL;
            else if ( !strcmp( type, "spherical" ) )   cur->jtype = JOINT_SPHERICAL;
            else if ( !strcmp( type, "floating" ) )    cur->jtype = JOINT_FLOATING;
            else {
                snprintf( err, errSize, "%s:%d: unknown joint type '%s'", path, lineNo, type );
                goto fail;
            }
            if ( cur->jtype == JOINT_SPHERICAL || cur->jtype == JOINT_FLOATING ) {
                cur->axis = 0;
            } else {
                if ( n < 2 || axis[1] != 0 || ( axis[0] != 'x' && axis[0] != 'y' && axis[0] != 'z' ) ) {
                    snprintf( err, errSize, "%s:%d: joint '%s' needs an axis of x, y, or z", path, lineNo, type );
                    goto fail;
                }
                cur->axis = axis[0];
            }
        } else if ( !strcmp( key, "tree" ) ) {
            if ( sscanf( rest, "%lf %lf %lf", &cur->tree[0], &cur->tree[1], &cur->tree[2] ) != 3 ) goto badArgs;
        } else if ( !strcmp( key, "mass" ) ) {
            if ( sscanf( rest, "%lf", &cur->mass ) != 1 ) goto badArgs;
        } else if ( !strcmp( key, "cg" ) ) {
            if ( sscanf( rest, "%lf %lf %lf", &cur->cg[0], &cur->cg[1], &cur->cg[2] ) != 3 ) goto badArgs;
        } else if ( !strcmp( key, "inertia" ) ) {
            if ( sscanf( rest, "%lf %lf %lf %lf %lf %lf",
                         &cur->inertia[0], &cur->inertia[1], &cur->inertia[2],
                         &cur->inertia[3], &cur->inertia[4], &cur->inertia[5] ) != 6 ) goto badArgs;
        } else if ( !strcmp( key, "fext" ) ) {
            if ( sscanf( rest, "%d", &cur->hasFext ) != 1 ) goto badArgs;
        } else {
            snprintf( err, errSize, "%s:%d: unknown keyword '%s'", path, lineNo, key );
            goto fail;
        }
        continue;

    badArgs:
        snprintf( err, errSize, "%s:%d: bad arguments for '%s'", path, lineNo, key );
        goto fail;
    }
    fclose( f );

    if ( m->nbod == 0 ) {
        snprintf( err, errSize, "%s: no bodies defined", path );
        return 1;
    }
    for ( int i = 0 ; i < m->nbod ; i++ ) {
        ModelBody *b = &m->body[ i ];
        if ( b->parent < -1 || b->parent >= i ) {
            snprintf( err, errSize, "%s: body %d: parent %d out of range (must be -1..%d, i.e. defined earlier)",
                      path, i, b->parent, i - 1 );
            return 1;
        }
        if ( b->mass <= 0 ) {
            snprintf( err, errSize, "%s: body %d: mass must be positive", path, i );
            return 1;
        }
        if ( b->jtype == JOINT_FLOATING && b->parent != -1 ) {
            snprintf( err, errSize, "%s: body %d: a floating joint must attach to the base (parent -1)",
                      path, i );
            return 1;
        }
        if ( b->jtype == JOINT_POINTONLINE || b->jtype == JOINT_DISTANCE ) {
            snprintf( err, errSize, "%s: body %d: %s is a loop joint type, not a tree joint",
                      path, i, b->jtype == JOINT_POINTONLINE ? "point-on-line" : "distance" );
            return 1;
        }
    }
    for ( int l = 0 ; l < m->nloop ; l++ ) {
        ModelLoop *lp = &m->loop[ l ];
        if ( lp->pred < -1 || lp->pred >= m->nbod || lp->succ < -1 || lp->succ >= m->nbod ) {
            snprintf( err, errSize, "%s: loop %d: predecessor/successor must both be set, to -1 (base) or a body index 0..%d",
                      path, l, m->nbod - 1 );
            return 1;
        }
        if ( lp->pred == lp->succ ) {
            snprintf( err, errSize, "%s: loop %d: predecessor and successor must differ", path, l );
            return 1;
        }
        if ( lp->jtype == JOINT_DISTANCE && lp->length <= 0 ) {
            snprintf( err, errSize, "%s: loop %d: a distance joint needs 'length <d>' with d > 0", path, l );
            return 1;
        }
    }
    if ( m->stab <= 0 ) {
        snprintf( err, errSize, "%s: stab must be positive", path );
        return 1;
    }
    for ( int p = 0 ; p < m->npoi ; p++ ) {
        ModelPoi *poi = &m->poi[ p ];
        if ( poi->body < -1 || poi->body >= m->nbod ) {
            snprintf( err, errSize, "%s: poi '%s': body %d out of range (must be -1..%d)",
                      path, poi->name, poi->body, m->nbod - 1 );
            return 1;
        }
        for ( int q = 0 ; q < p ; q++ )
            if ( !strcmp( m->poi[ q ].name, poi->name ) ) {
                snprintf( err, errSize, "%s: poi name '%s' declared more than once", path, poi->name );
                return 1;
            }
    }
    for ( int e = 0 ; e < m->nelem ; e++ ) {
        ModelElem *el = &m->elem[ e ];
        for ( int q = 0 ; q < e ; q++ )
            if ( m->elem[ q ].type == el->type && !strcmp( m->elem[ q ].name, el->name ) ) {
                snprintf( err, errSize, "%s: %s name '%s' declared more than once", path,
                          el->type == ELEM_TSDA ? "tsda" : "rsda", el->name );
                return 1;
            }
        if ( el->type == ELEM_TSDA ) {
            if ( el->poi1 < 0 || el->poi2 < 0 ) {
                snprintf( err, errSize, "%s: tsda '%s': needs both 'from' and 'to'", path, el->name );
                return 1;
            }
            if ( el->poi1 == el->poi2 ) {
                snprintf( err, errSize, "%s: tsda '%s': 'from' and 'to' must differ", path, el->name );
                return 1;
            }
            ModelPoi *p1 = &m->poi[ el->poi1 ], *p2 = &m->poi[ el->poi2 ];
            if ( p1->body == p2->body && p1->p[0] == p2->p[0] && p1->p[1] == p2->p[1] && p1->p[2] == p2->p[2] ) {
                snprintf( err, errSize, "%s: tsda '%s': anchors are coincident (same body, same point) -- l is always 0",
                          path, el->name );
                return 1;
            }
        } else { // ELEM_RSDA
            if ( el->joint < 0 || el->joint >= m->nbod ) {
                snprintf( err, errSize, "%s: rsda '%s': joint %d out of range (must be 0..%d)",
                          path, el->name, el->joint, m->nbod - 1 );
                return 1;
            }
            ModelJointType jt = m->body[ el->joint ].jtype;
            if ( jt != JOINT_REVOLUTE && jt != JOINT_PRISMATIC ) {
                snprintf( err, errSize, "%s: rsda '%s': joint on body %d must be revolute or prismatic (single-DOF)",
                          path, el->name, el->joint );
                return 1;
            }
        }
    }

    // Expand track bands into shoe bodies + a closure loop. Done after the
    // hand-written bodies/loops are validated and before coordinate layout, so
    // the generated shoes get counted into nq/nv. The curved design pose is
    // recorded here (the band cannot sit at q=0) and written to q0 after layout.
    // Arrays are declared without initializers so the forward goto to fail is
    // legal in C++.
    double trkAngle[ MODEL_MAX_BODIES ];               // pin design turn angle
    double trkShoe0Pose[ MODEL_MAX_TRACKS ][ 7 ];      // floating anchor pose
    TrkV2 trkPin[ MODEL_MAX_BODIES ];                  // scratch: pin positions
    for ( int i = 0 ; i < MODEL_MAX_BODIES ; i++ ) trkAngle[ i ] = 0;

    for ( int t = 0 ; t < m->ntrack ; t++ ) {
        ModelTrack *tr = &m->track[ t ];
        if ( tr->root < 0 || tr->root >= m->nbod ) {
            snprintf( err, errSize, "%s: track '%s': root %d out of range (0..%d)", path, tr->name, tr->root, m->nbod - 1 );
            goto fail;
        }
        if ( tr->shoes < 3 ) {
            snprintf( err, errSize, "%s: track '%s': need at least 3 shoes (got %d)", path, tr->name, tr->shoes );
            goto fail;
        }
        if ( tr->mass <= 0 ) {
            snprintf( err, errSize, "%s: track '%s': shoe mass must be positive", path, tr->name );
            goto fail;
        }
        int endIdx[ 2 ], nend = 0;
        for ( int r = 0 ; r < tr->nroller ; r++ ) {
            ModelTrackRoller *rr = &tr->roller[ r ];
            if ( rr->body < 0 || rr->body >= m->nbod ) {
                snprintf( err, errSize, "%s: track '%s': roller body %d out of range (0..%d)", path, tr->name, rr->body, m->nbod - 1 );
                goto fail;
            }
            if ( rr->radius <= 0 ) {
                snprintf( err, errSize, "%s: track '%s': roller on body %d needs a positive radius", path, tr->name, rr->body );
                goto fail;
            }
            if ( rr->isEnd ) { if ( nend < 2 ) endIdx[ nend ] = r; nend++; }
        }
        if ( nend != 2 ) {
            snprintf( err, errSize, "%s: track '%s': need exactly 2 'end' rollers (sprocket + idler), got %d", path, tr->name, nend );
            goto fail;
        }
        double pA[ 3 ], pB[ 3 ];
        if ( trkBodyPosInRoot( m, tr->roller[ endIdx[0] ].body, tr->root, pA ) ||
             trkBodyPosInRoot( m, tr->roller[ endIdx[1] ].body, tr->root, pB ) ) {
            snprintf( err, errSize, "%s: track '%s': an end pulley is not rigidly located in the root frame", path, tr->name );
            goto fail;
        }
        TrkV2 cA, cB; cA.x = pA[0]; cA.z = pA[2]; cB.x = pB[0]; cB.z = pB[2];
        double rA = tr->roller[ endIdx[0] ].radius, rB = tr->roller[ endIdx[1] ].radius;
        double Dc = sqrt( ( cB.x - cA.x )*( cB.x - cA.x ) + ( cB.z - cA.z )*( cB.z - cA.z ) );
        if ( Dc <= fabs( rA - rB ) ) {
            snprintf( err, errSize, "%s: track '%s': end pulleys overlap (centers too close for the belt)", path, tr->name );
            goto fail;
        }
        int N = tr->shoes;
        if ( m->nbod + N > MODEL_MAX_BODIES ) {
            snprintf( err, errSize, "%s: track '%s': %d shoes overflow MODEL_MAX_BODIES (%d)", path, tr->name, N, MODEL_MAX_BODIES );
            goto fail;
        }
        if ( m->nloop + 1 > MODEL_MAX_LOOPS ) {
            snprintf( err, errSize, "%s: track '%s': closure loop overflows MODEL_MAX_LOOPS (%d)", path, tr->name, MODEL_MAX_LOOPS );
            goto fail;
        }

        TrkSeg seg[ 4 ];
        double P = trkBuildBelt( cA, rA, cB, rB, seg );
        double ds = P / N;                             // even spacing == shoe length
        double yband = pA[1];
        tr->linkLen = ds;
        for ( int i = 0 ; i < N ; i++ ) { double ph; trkPin[ i ] = trkPoint( seg, i*ds, &ph ); }

        int b0 = m->nbod;
        tr->shoe0 = b0;
        // Anchor shoe: a floating body (held in the plane by contact at runtime,
        // not a rigid mount), parented to the base like the quartercar spindle.
        ModelBody *s0 = &m->body[ m->nbod++ ];
        memset( s0, 0, sizeof *s0 );
        s0->parent = -1; s0->jtype = JOINT_FLOATING; s0->axis = 0;
        s0->mass = tr->mass;
        s0->cg[0] = tr->cg[0]; s0->cg[1] = tr->cg[1]; s0->cg[2] = tr->cg[2];
        for ( int j = 0 ; j < 6 ; j++ ) s0->inertia[ j ] = tr->inertia[ j ];
        s0->hasFext = tr->hasFext;

        // phi = link direction; joint angle about y is q = -phi (right-handed
        // R_y maps +x toward -z; see the joint-convention note at the top).
        TrkV2 l0; l0.x = trkPin[1].x - trkPin[0].x; l0.z = trkPin[1].z - trkPin[0].z;
        double phiPrev = atan2( l0.z, l0.x );
        double q0abs = -phiPrev;
        trkShoe0Pose[ t ][ 0 ] = cos( q0abs / 2 );     // quaternion about y
        trkShoe0Pose[ t ][ 1 ] = 0;
        trkShoe0Pose[ t ][ 2 ] = sin( q0abs / 2 );
        trkShoe0Pose[ t ][ 3 ] = 0;
        trkShoe0Pose[ t ][ 4 ] = trkPin[0].x;          // base-frame position (root
        trkShoe0Pose[ t ][ 5 ] = yband;                // is at identity in design
        trkShoe0Pose[ t ][ 6 ] = trkPin[0].z;          // pose, so base == root)

        for ( int i = 1 ; i < N ; i++ ) {
            ModelBody *sh = &m->body[ m->nbod++ ];
            memset( sh, 0, sizeof *sh );
            sh->parent = b0 + i - 1; sh->jtype = JOINT_REVOLUTE; sh->axis = 'y';
            sh->tree[0] = ds; sh->tree[1] = 0; sh->tree[2] = 0;
            sh->mass = tr->mass;
            sh->cg[0] = tr->cg[0]; sh->cg[1] = tr->cg[1]; sh->cg[2] = tr->cg[2];
            for ( int j = 0 ; j < 6 ; j++ ) sh->inertia[ j ] = tr->inertia[ j ];
            sh->hasFext = tr->hasFext;
            int inext = ( i + 1 ) % N;
            TrkV2 lk; lk.x = trkPin[inext].x - trkPin[i].x; lk.z = trkPin[inext].z - trkPin[i].z;
            double phi = atan2( lk.z, lk.x );
            double turn = -( phi - phiPrev );           // relative turn (q=-phi)
            while ( turn >   M_PI ) turn -= 2*M_PI;      // principal value: chord
            while ( turn <= -M_PI ) turn += 2*M_PI;      // directions wrap at +/-pi
            trkAngle[ b0 + i ] = turn;
            phiPrev = phi;
        }

        // Close the loop: last shoe's trailing pin back to shoe0's leading pin,
        // point-on-line about y (constrains the 2 in-plane translations, leaves
        // the pin rotation free) -- the planar-loop rule, as in fourbar.mech.
        ModelLoop *lp = &m->loop[ m->nloop++ ];
        memset( lp, 0, sizeof *lp );
        lp->pred = b0 + N - 1; lp->succ = b0;
        lp->jtype = JOINT_POINTONLINE; lp->axis = 'y';
        lp->ptree[0] = ds; lp->ptree[1] = 0; lp->ptree[2] = 0;
        lp->stree[0] = 0;  lp->stree[1] = 0; lp->stree[2] = 0;
    }

    // Position/velocity coordinate layout: floating bodies take 7 q slots
    // and 6 qd slots, everything else takes 1 of each. (Declared without
    // initializers so the error-path goto above is legal in C++ too.)
    int qi, vi;
    qi = 0;
    vi = 0;
    for ( int i = 0 ; i < m->nbod ; i++ ) {
        m->qidx[ i ] = qi;
        m->vidx[ i ] = vi;
        qi += modelBodyNq( &m->body[ i ] );
        vi += modelBodyNv( &m->body[ i ] );
    }
    m->nq = qi;
    m->nv = vi;

    // Design-pose coordinates. Hand-written bodies sit at 0 (floating ones at
    // quaternion identity); generated shoes carry the curved band pose computed
    // above. Consumers seed their initial q from m->q0.
    for ( int k = 0 ; k < m->nq ; k++ ) m->q0[ k ] = 0;
    for ( int i = 0 ; i < m->nbod ; i++ )
        if ( m->body[ i ].jtype == JOINT_FLOATING ) m->q0[ m->qidx[ i ] ] = 1.0;
    for ( int t = 0 ; t < m->ntrack ; t++ ) {
        int b0 = m->track[ t ].shoe0, N = m->track[ t ].shoes;
        if ( b0 < 0 ) continue;
        for ( int j = 0 ; j < 7 ; j++ ) m->q0[ m->qidx[ b0 ] + j ] = trkShoe0Pose[ t ][ j ];
        for ( int i = 1 ; i < N ; i++ ) m->q0[ m->qidx[ b0 + i ] ] = trkAngle[ b0 + i ];
    }

    // Constraint-row layout: loop l owns rows cidx[l] .. cidx[l]+nc-1 of K.
    int ci;
    ci = 0;
    for ( int l = 0 ; l < m->nloop ; l++ ) {
        m->cidx[ l ] = ci;
        ci += modelLoopNc( &m->loop[ l ] );
    }
    m->nc = ci;

    // Observation layout: all POIs in declaration order (6 doubles each),
    // then all TSDAs (2 each: l ldot), then all RSDAs (2 each: angle rate) --
    // grouped by kind, not interleaved declaration order, so runtime lookup
    // tables stay contiguous per kind.
    for ( int p = 0 ; p < m->npoi ; p++ ) m->poi[ p ].obsIdx = 6 * p;
    int obsIdx;
    obsIdx = 6 * m->npoi;
    for ( int e = 0 ; e < m->nelem ; e++ )
        if ( m->elem[ e ].type == ELEM_TSDA ) { m->elem[ e ].obsIdx = obsIdx; obsIdx += 2; }
    for ( int e = 0 ; e < m->nelem ; e++ )
        if ( m->elem[ e ].type == ELEM_RSDA ) { m->elem[ e ].obsIdx = obsIdx; obsIdx += 2; }
    m->obsLen = obsIdx;

    // act[] slots: assigned in declaration order across both element kinds.
    int actIdx;
    actIdx = 0;
    for ( int e = 0 ; e < m->nelem ; e++ )
        if ( m->elem[ e ].hasInput ) m->elem[ e ].actIdx = actIdx++;
    m->nact = actIdx;

    return 0;

fail:
    fclose( f );
    return 1;
}

// Returns NULL if every joint is a type the current single-DOF pipeline
// (one q per body) can handle, otherwise a message. Spherical (3 DOF) and
// cylindrical (2 DOF) parse fine but need multi-DOF q indexing that neither
// the generator nor the dense reference implements yet.
static const char *modelUnsupported( const Model *m ) {
    for ( int i = 0 ; i < m->nbod ; i++ ) {
        if ( m->body[i].jtype == JOINT_SPHERICAL || m->body[i].jtype == JOINT_CYLINDRICAL )
            return "spherical/cylindrical joints are parsed but not implemented yet (multi-DOF)";
    }
    return NULL;
}
