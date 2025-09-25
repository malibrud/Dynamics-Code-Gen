#pragma once
#include <inttypes.h>
#include "svec.h"
#include "vec3.h"
#include "frame.h"
#include "body.h"
#include "joint.h"
#include "xform.h"
#include "arena.h"

typedef struct {
    int nBod;
    int16_t *parent;
    int16_t *children;
    uint8_t *jtype;
    uint8_t *jdof;
    Body *bodies;
    Joint *joints;
    int16_t *qidx;
    int16_t *qdidx;

    double *q;
    double *qd;
    double *qdd;
    Xform **xTree;
    Xform **xJoint;
    Xform **xFromParent;
    Xform **xFromBase;
    Svec *vel;
    Svec *acc;
    Svec *force;

    Arena *arena;
} System;

void sys_updateJointXform( System *s, int jidx ) {
    int qidx = s->qidx[ jidx ];
    x_update( s->xJoint, &s->q[ qidx ] );
}

#define MB ( 1024 * 1024 )
void sys_init( System *sys, int nBod ) {
    assert( sys );
    assert( nBod > 0 );

    int N = nBod;
    Arena *ar = ar_create( 1*MB );

    sys->nBod = nBod;
    sys->parent   = ar_alloc_array( ar, N, sizeof( int16_t ) );
    sys->children = ar_alloc_array( ar, N, sizeof( int16_t ) );

    sys->jtype = ar_alloc_array( ar, N, sizeof( int8_t ) );
    sys->jdof  = ar_alloc_array( ar, N, sizeof( int8_t ) );

    sys->bodies = ar_alloc_array( ar, N, sizeof( Body ) );
    sys->joints = ar_alloc_array( ar, N, sizeof( Joint ) );

    sys->qidx  = ar_alloc_array( ar, N, sizeof( int16_t ) );
    sys->qdidx = ar_alloc_array( ar, N, sizeof( int16_t ) );

    sys->q   = ar_alloc_array( ar, N, sizeof( double ) );
    sys->qd  = ar_alloc_array( ar, N, sizeof( double ) );
    sys->qdd = ar_alloc_array( ar, N, sizeof( double ) );

    sys->xTree       = ar_alloc_array( ar, N, sizeof( Xform* ) );
    sys->xJoint      = ar_alloc_array( ar, N, sizeof( Xform* ) );
    sys->xFromParent = ar_alloc_array( ar, N, sizeof( Xform* ) );
    sys->xFromBase   = ar_alloc_array( ar, N, sizeof( Xform* ) );

    sys->vel   = ar_alloc_array( ar, N, sizeof( Svec ) );
    sys->acc   = ar_alloc_array( ar, N, sizeof( Svec ) );
    sys->force = ar_alloc_array( ar, N, sizeof( Svec ) );

    sys->arena = ar;
}

Xform *sys_newXform( System *sys, int type ) {
    Xform *X;
    switch ( type )
    {
    case XFORM_IDENT:
        X = ar_alloc_aligned( sys->arena, sizeof( Xform), sizeof( void* ) );
        X->type = type;
        return X;
    
    default:
        break;
    }

}