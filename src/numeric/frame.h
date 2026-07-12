#pragma once
#include "vec3.h"

typedef struct {
    Vec3 O;
    Vec3 x;
    Vec3 y;
    Vec3 z;
} Frame;

void fr_setToIdentity( Frame *f ) {
    f->x = (Vec3){ 1, 0, 0 };
    f->y = (Vec3){ 0, 1, 0 };
    f->z = (Vec3){ 0, 0, 1 };
}