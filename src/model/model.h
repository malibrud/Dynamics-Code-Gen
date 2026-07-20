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

#define MODEL_MAX_BODIES 32
#define MODEL_MAX_Q  ( MODEL_MAX_BODIES * 7 ) // worst case: every body floating
#define MODEL_MAX_V  ( MODEL_MAX_BODIES * 6 )

#define MODEL_NAME_LEN 32
#define MODEL_MAX_POI 16
#define MODEL_MAX_ELEM 16

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

#define MODEL_MAX_LOOPS 8
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
