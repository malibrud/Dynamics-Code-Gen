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
//   joint <type> [axis]          type: revolute | prismatic | cylindrical | spherical
//                                axis: x | y | z (required except for spherical)
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
// All joint types parse; which ones a given consumer implements is up to it
// (see modelUnsupported below for the current single-DOF pipeline).

#pragma once
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MODEL_MAX_BODIES 32

typedef enum {
    JOINT_REVOLUTE,
    JOINT_PRISMATIC,
    JOINT_CYLINDRICAL,
    JOINT_SPHERICAL,
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

typedef struct {
    int nbod;
    double gravity[3];     // physical gravity
    ModelBody body[ MODEL_MAX_BODIES ];
} Model;

// Loads a mechanism file. Returns 0 on success; on failure returns nonzero
// and writes a message (with file/line context) into err.
static int modelLoad( const char *path, Model *m, char *err, size_t errSize ) {
    FILE *f = fopen( path, "r" );
    if ( !f ) {
        snprintf( err, errSize, "%s: cannot open file", path );
        return 1;
    }

    memset( m, 0, sizeof *m );
    ModelBody *cur = NULL;
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
        if ( !strcmp( key, "body" ) ) {
            if ( m->nbod >= MODEL_MAX_BODIES ) {
                snprintf( err, errSize, "%s:%d: too many bodies (max %d)", path, lineNo, MODEL_MAX_BODIES );
                goto fail;
            }
            cur = &m->body[ m->nbod++ ];
            memset( cur, 0, sizeof *cur );
            cur->parent = -1;
            cur->jtype = JOINT_REVOLUTE;
            cur->axis = 'z';
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
            else {
                snprintf( err, errSize, "%s:%d: unknown joint type '%s'", path, lineNo, type );
                goto fail;
            }
            if ( cur->jtype == JOINT_SPHERICAL ) {
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
        if ( m->body[i].parent < -1 || m->body[i].parent >= i ) {
            snprintf( err, errSize, "%s: body %d: parent %d out of range (must be -1..%d, i.e. defined earlier)",
                      path, i, m->body[i].parent, i - 1 );
            return 1;
        }
        if ( m->body[i].mass <= 0 ) {
            snprintf( err, errSize, "%s: body %d: mass must be positive", path, i );
            return 1;
        }
    }
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
