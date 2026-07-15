#pragma once
#include <vector>
#include <inttypes.h>
#include <assert.h>
#include "joint.h"
#include "xform.h"
#include "svec.h"
#include "inertia.h"

using namespace std;

struct System{
    int nBod;
    vector<int16_t> parent;
    vector<int16_t> children;
    vector<uint8_t> jtype;     // J_REVOLUTE, J_PRISMATIC, or J_FLOATING (joint.h)
    vector<char> jaxis;        // 'x', 'y', or 'z' (unused for floating)

    // Coordinate layout: floating joints take 7 q slots (quaternion + position)
    // and 6 velocity slots; single-DOF joints take 1 of each, so nq != nv when
    // any body floats.
    vector<int> qidx;          // body i's first slot in the q vector
    vector<int> vidx;          // body i's first slot in qd / tau / out_qdd
    int nq = 0;
    int nv = 0;

    vector<Svec> S;
    vector<Svec> vJ;
    vector<Svec> v;
    vector<Svec> a;
    vector<Svec> force;            // spatial bias force f_i, accumulated by the RNEA inward pass
    vector<Xform> xTree;
    vector<Xform> xJoint;
    vector<Xform> xFromParent;
    vector<Xform> xFromBase;

    vector<Body> bodies;           // numeric physical description, set by the caller
    vector<SpatialInertia> I;      // body i's own spatial inertia, about its own frame
    vector<SpatialInertia> Ic;     // composite (subtree) spatial inertia, accumulated by CRBA

    // Per-body flag: generate code to apply an external spatial force to this
    // body (read from the fext argument of the generated function). Bodies
    // left at 0 contribute no code and their fext slots are never read --
    // the same structural-sparsity idea as spatial_v2's empty f_ext{i} cells.
    vector<uint8_t> hasFext;

    // Base pseudo-acceleration used to seed the RNEA acceleration recursion,
    // i.e. -a_grav in spatial_v2's HandC.m (avp{i} = Xup{i} * -a_grav for the
    // root body) -- so this field should already hold the *negated* physical
    // gravity vector, set by the caller.
    Svec grav;

    System(): System(0) {}
    System( int N )
        :nBod( N )
        ,parent( N )
        ,children( N )
        ,jtype( N )
        ,jaxis( N, 'z' )
        ,qidx( N )
        ,vidx( N )
        ,S( N )
        ,vJ( N )
        ,v( N )
        ,a( N )
        ,force( N )
        ,xTree( N )
        ,xJoint( N )
        ,xFromParent( N )
        ,xFromBase( N )
        ,bodies( N )
        ,I( N )
        ,Ic( N )
        ,hasFext( N )
        ,grav( "grav"s )
    {
    }
};
