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

    // Kinematic loops (RBDA ch. 8). Loop joints add constraints, never
    // coordinates. gen_loops() fills Kentry/kfulName with the emitted temp
    // names (or "0" for structural zeros) that gen_solve_qdd() assembles
    // into the method-2 solve.
    struct SysLoop {
        int pred, succ;        // body indices; -1 = fixed base
        uint8_t jtype;         // J_REVOLUTE .. J_DISTANCE
        char axis;             // 'x' | 'y' | 'z'; unused for spherical and distance
        double ptree[3];       // predecessor body frame -> loop joint frame
        double stree[3];       // successor body frame -> loop joint frame
        double length;         // distance joints: imposed anchor-to-anchor distance
    };
    vector<SysLoop> loops;
    double stab = 0.1;             // Baumgarte Tstab; alpha = beta = 1/stab
    int nc = 0;                    // total constraint rows, set by gen_loops()
    vector<vector<string>> Kentry; // nc x nv expressions ("0" = structural zero)
    vector<string> kfulName;       // nc RHS temps (k + kstab of RBDA eq 8.35)

    // Points of interest and force elements (TSDA/RSDA). A poi is a named
    // point fixed to a body ("body" -1 = fixed base); a TSDA anchor *is* a
    // poi. See src/model/model.h for the grammar and force-law signs.
    struct SysPoi {
        string name;
        int body;           // -1 = fixed base
        double p[3];
        int obsIdx;          // offset into out_obs (6 doubles)
    };
    struct SysElem {
        enum Kind { TSDA, RSDA };
        string name;
        Kind type;
        int poi1, poi2;     // TSDA: indices into pois[] ("from"/"to")
        int joint;           // RSDA: body index owning the tree joint
        double k, c, ref;    // ref = l0 (TSDA) or q0 (RSDA)
        int hasInput;
        int actIdx;           // slot in act[], -1 if !hasInput
        int obsIdx;            // offset into out_obs (2 doubles)
    };
    vector<SysPoi> pois;
    vector<SysElem> elems;
    int nact = 0;                  // total act[] slots, set by the caller from Model

    // Element-applied generalized force per velocity coordinate, filled in
    // by gen_elems(): "0" (or an empty vector, when there are no elements at
    // all) means no element touches that coordinate -- gen_solve_qdd() skips
    // it entirely so element-free models generate byte-identical rhs code.
    vector<string> Qentry;

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
