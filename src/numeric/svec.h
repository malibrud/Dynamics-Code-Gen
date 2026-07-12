#pragma once

typedef struct {
    double v[6];
} Svec;

// Computes vc = va + vb;
void sv_add( Svec *va, Svec *vb, Svec *vc ) {
    double *a = va->v;
    double *b = vb->v;
    double *c = vc->v;
    for ( int i = 0 ; i < 6 ; i++ ) c[i] = a[i] + b[i];
}

// Computes vc = va - vb;
void sv_sub( Svec *va, Svec *vb, Svec *vc ) {
    double *a = va->v;
    double *b = vb->v;
    double *c = vc->v;
    for ( int i = 0 ; i < 6 ; i++ ) c[i] = a[i] - b[i];
}

// Computes vc += va;
void sv_addTo( Svec *vc, Svec *va ) {
    sv_add( vc, va, vc );
}

// Computes vc -= va;
void sv_subFrom( Svec *vc, Svec *va ) {
    sv_sub( vc, va, vc );
}

// Computes the motion cross product, vc = va x vb;
// See RBDA, pp 25, eq (2.33)
void sv_mcross( Svec *va, Svec *vb, Svec *vc ) {
    double *w  = va->v;
    double *v0 = va->v + 3;
    double *m  = vb->v;
    double *m0 = vb->v + 3;
    double *c  = vc->v;

    // w x m
    c[0] = w[1]*m[2] - w[2]*m[1];
    c[1] = w[2]*m[0] - w[0]*m[2];
    c[2] = w[0]*m[1] - w[1]*m[0];

    // w x m0 + v0 x m
    c[3] = w[1]*m0[2] - w[2]*m0[1]  +  v0[1]*m[2] - v0[2]*m[1];
    c[4] = w[2]*m0[0] - w[0]*m0[2]  +  v0[2]*m[0] - v0[0]*m[2];
    c[5] = w[0]*m0[1] - w[1]*m0[0]  +  v0[0]*m[1] - v0[1]*m[0];
}

// Computes the force cross product, vc = va x* vb;
// See RBDA, pp 25, eq (2.34)
void sv_fcross( Svec *va, Svec *vb, Svec *vc ) {
    double *w  = va->v;
    double *v0 = va->v + 3;
    double *f0  = vb->v;
    double *f = vb->v + 3;
    double *c  = vc->v;

    // w x f0 + v0 x f
    c[0] = w[1]*f0[2] - w[2]*f0[1]  +  v0[1]*f[2] - v0[2]*f[1];
    c[1] = w[2]*f0[0] - w[0]*f0[2]  +  v0[2]*f[0] - v0[0]*f[2];
    c[2] = w[0]*f0[1] - w[1]*f0[0]  +  v0[0]*f[1] - v0[1]*f[0];

    // w x f
    c[3] = w[1]*f[2] - w[2]*f[1];
    c[4] = w[2]*f[0] - w[0]*f[2];
    c[5] = w[0]*f[1] - w[1]*f[0];
}