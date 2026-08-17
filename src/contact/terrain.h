// Terrain: the ground surface, as a HEIGHT FIELD h(x,y).
//
// This file is the ONLY thing that differs between flat ground and uneven
// ground. The force law in ground_contact.h reads nothing but h and its two
// first derivatives, so a slope, a sinusoid or a sampled grid are new cases
// HERE and change nothing downstream -- not the force law, not the wrench
// bookkeeping, not the energy budget, not the tests. Flat ground is not a
// special case at all; it is a plane with a vertical normal.
//
// WHY A HEIGHT FIELD AND NOT A GENERAL SURFACE. Two reasons, and the second is
// the one that matters.
//
//   1. Lookup is O(1) with no broadphase, no contact manifold and no persistent
//      pair set -- which is what keeps the ground engine stateless, and
//      statelessness is what lets src/contact/main.c exercise it at random
//      off-trajectory states. A triangle mesh would need all three.
//
//   2. It makes the contact law EXACTLY CONSERVATIVE on arbitrary terrain, for
//      free. Take the vertical gap
//
//          g(p) = h(p.x, p.y) - p.z          ( > 0 = penetrating )
//
//      as the strain measure and U = 1/2 k g^2 as the potential. Then
//
//          grad g = ( h_x, h_y, -1 )
//          -grad U = -k g grad g = k g ( -h_x, -h_y, 1 )
//
//      which is the elastic force written down directly -- no planarity
//      assumption, and only FIRST derivatives of h are needed. Note the
//      gradient is the UNNORMALIZED ( -h_x, -h_y, 1 ); normalizing it is
//      exactly what would break the identity.
//
//      Contrast the obvious alternative, "penetration = perpendicular distance
//      to the surface, force along the unit normal". That is exact on a plane
//      and NOT the gradient of anything on curved terrain, so it would leak
//      energy the same way the sprocket tooth law did before it grew its
//      g'(delta)*U0 radial completion (track_contact.h) -- a bug that took two
//      rounds of trajectory-energy debugging to find. Getting it structurally
//      right here costs nothing.
//
// The price is that terrain cannot be vertical or overhanging. For ground under
// a vehicle that is not a restriction; modelLoad rejects nz <= 0 at parse time.

#pragma once
#include <math.h>
#include "../model/model.h"

// Evaluates the terrain under (x,y): height and its two first derivatives.
// Returns 0 if there is no terrain, in which case the outputs are zeroed and
// the caller must not treat h as a surface.
//
// dh/dx and dh/dy are ANALYTIC, never differenced. src/contact/main.c checks
// them against central differences of h itself, which catches a gradient bug in
// every terrain kind for free -- worth having, because an inconsistent gradient
// would break the conservativeness argument above silently.
static int terrainAt( const ModelTerrain *tn, double x, double y,
                      double *h, double *hx, double *hy ) {
    *h = *hx = *hy = 0;
    switch ( tn->kind ) {
        case TERRAIN_PLANE: {
            // n.p = d  =>  h = ( d - nx*x - ny*y ) / nz. Guaranteed nz > 0 by
            // the parser, so no division guard is needed here.
            double inv = 1.0 / tn->n[2];
            *h  = ( tn->d - tn->n[0]*x - tn->n[1]*y ) * inv;
            *hx = -tn->n[0] * inv;
            *hy = -tn->n[1] * inv;
            return 1;
        }
        case TERRAIN_NONE:
        default:
            return 0;
    }
}
