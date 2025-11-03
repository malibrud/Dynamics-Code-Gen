#pragma once
#include <inttypes.h>
#include <assert.h>
#include "svec.h"
#include "vec3.h"
#include "frame.h"
#include "inertia.h"
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
    // Joint *joints;
    int16_t *qidx;
    int16_t *qdidx;

    double *q;
    double *qd;
    double *qdd;
    double *tau;
    double *C;
    double *H;
    Svec *S;
    Xform **xTree;
    Xform **xJoint;
    Xform **xFromParent;
    Xform **xFromBase;
    SpatialInertia *I0; // Spatial inertia of body in body coordinates, see RBDA, pp 33, eq. (2.63)
    SpatialInertia *Ic; // Spatial inertia of the composite rigid body in body coordinates, see RBDA, pp 106, eq. (6.17)
    Svec *vel;
    Svec *acc;
    Svec *force;
    Svec *fext;
    Svec grav;

    Arena *arena;
} System;

void sys_updateJointXform( System *s, int jidx ) {
    int qidx = s->qidx[ jidx ];
    x_update( s->xJoint[ jidx ], &s->q[ qidx ] );
}

#define MB ( 1024 * 1024 )
void sys_init( System *sys, int nBod, Svec gravity ) {
    assert( sys );
    assert( nBod > 0 );

    int N = nBod;
    Arena *ar = ar_create( 1*MB );

    sys->grav = gravity;
    sys->nBod = nBod;
    sys->parent   = ar_alloc_array( ar, N, sizeof( int16_t ) );
    sys->children = ar_alloc_array( ar, N, sizeof( int16_t ) );

    sys->jtype = ar_alloc_array( ar, N, sizeof( int8_t ) );
    sys->jdof  = ar_alloc_array( ar, N, sizeof( int8_t ) );

    sys->bodies = ar_alloc_array( ar, N, sizeof( Body ) );
    //sys->joints = ar_alloc_array( ar, N, sizeof( Joint ) );

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
