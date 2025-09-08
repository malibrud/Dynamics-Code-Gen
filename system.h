#pragma once
#include <inttypes.h>
#include "svec.h"
#include "body.h"
#include "joint.h"
#include "xform.h"

typedef struct {
    int16_t nBod;
    int16_t *parent;
    int16_t *children;
    Body *bodies;
    Xform *xFromParent;
    Xform *xFromBase;
    Svec *vel;
    Svec *acc;
    Svec *force;
} System;