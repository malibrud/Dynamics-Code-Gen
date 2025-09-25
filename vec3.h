#pragma once
#include <stdbool.h>
#include <math.h>

typedef struct {
    double x, y, z;
} Vec3;

Vec3 v3_sub( Vec3 *a, Vec3 *b) {
    return (Vec3){ 
        a->x - b->x, 
        a->y - b->y, 
        a->z - b->z
    };
}

bool v3_isZero( Vec3 *v ) {
    return 
        v->x == 0.0 &&
        v->y == 0.0 &&
        v->z == 0.0
    ;
}

Vec3 v3_cross( Vec3 *a, Vec3 *b ) {
    Vec3 c;
    c.x = a->y * b->z - a->z * b->y;
    c.y = a->z * b->x - a->x * b->z;
    c.z = a->x * b->y - a->y * b->x;
    return c;
}

double v3_dot( Vec3 *a, Vec3 *b ) {
    double c;
    c = a->x * b->x + 
        a->y * b->y + 
        a->z * b->z;
    return c;
}

double v3_dist( Vec3 *a, Vec3 *b ) {
    double c = 0;
    double d;
    d = a->x - b->x; c += d*d;
    d = a->y - b->y; c += d*d;
    d = a->z - b->z; c += d*d;
    return sqrt( c );
}