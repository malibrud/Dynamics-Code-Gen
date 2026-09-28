r"""Plot the time histories written by sim_<model>.exe --trace.

    build\sim_m113.exe ..\models\m113.mech --trace m113_trace.csv
    python tools\plot_trace.py build\m113_trace.csv

Produces three PNGs beside the CSV:

  <stem>_settling.png   hull height, attitude, ground load, penetration,
                        loaded-contact counts, horizontal creep
  <stem>_band.png       the band in the x-z plane at the first and last
                        sample, with its rollers and the terrain
  <stem>_rollers.png    every roller's height against time, and the band's
                        own height envelope over all shoes

The trace carries its own constants (weight, roller radii, terrain plane) in
the leading '#' lines, so nothing here has to be told which model it came
from -- see traceHeader in src/sim/main.c.

Colors are the first slots of a validated categorical palette, assigned in a
fixed order and never cycled; --dark restages the same hues for a dark
surface. No panel carries two y-scales.
"""

import argparse
import os
import re
import sys

import matplotlib
import numpy as np

# --- palette -------------------------------------------------------------
# Fixed slot order (blue, orange, aqua, yellow, magenta, green, violet, red);
# a series takes the next unused slot and keeps it. The light and dark rows
# are the same eight hues stepped for their own surface.
SERIES_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SERIES_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500",
               "#d55181", "#008300", "#9085e9", "#e66767"]
INK_LIGHT = ("#fcfcfb", "#0b0b0b", "#52514e", "#d8d7d2")
INK_DARK = ("#1a1a19", "#ffffff", "#c3c2b7", "#3a3a38")


def style(dark):
    surface, primary, secondary, gridc = INK_DARK if dark else INK_LIGHT
    matplotlib.rcParams.update({
        "figure.facecolor": surface, "axes.facecolor": surface,
        "savefig.facecolor": surface,
        "text.color": primary, "axes.labelcolor": secondary,
        "xtick.color": secondary, "ytick.color": secondary,
        "axes.edgecolor": gridc, "grid.color": gridc,
        "axes.titlecolor": primary,
        "axes.grid": True, "grid.linewidth": 0.6, "grid.alpha": 0.7,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.6, "font.size": 9,
        "axes.titlesize": 10, "legend.frameon": False, "legend.fontsize": 8,
    })
    return (SERIES_DARK if dark else SERIES_LIGHT), primary, secondary, gridc


# --- the file ------------------------------------------------------------

def load(path):
    """Returns (columns dict, metadata dict). The '#' preamble is metadata,
    the first bare line is the column names, the rest is the run."""
    meta = {"tracks": [], "terrain": None}
    names, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                body = line[1:].strip()
                if body.startswith("track "):
                    t = dict(re.findall(r"(\w+)=(\S+)", body))
                    t["name"] = body.split()[1]
                    t["shoes"] = int(t["shoes"])
                    t["radii"] = [(float(r.rstrip("de")),
                                   r[-1] if r[-1] in "de" else "")
                                  for r in t["radii"].split(";")]
                    for k in ("linkLen", "padZ", "thickness"):
                        t[k] = float(t[k])
                    meta["tracks"].append(t)
                elif body.startswith("terrain plane"):
                    kv = dict(re.findall(r"(\w+)=(\S+)", body))
                    meta["terrain"] = ([float(v) for v in kv["n"].split(";")],
                                       float(kv["d"]))
                else:
                    # "tracks=2" on the model line is a count; the track list
                    # is built from the per-track lines below it.
                    for k, v in re.findall(r"(\w+)=(\S+)", body):
                        meta["ntracks" if k == "tracks" else k] = v
            elif names is None:
                names = line.strip().split(",")
            else:
                rows.append(line)
    if names is None or not rows:
        sys.exit("%s: no data rows -- was the traced run a track model?" % path)
    data = np.loadtxt(rows, delimiter=",")
    return {n: data[:, i] for i, n in enumerate(names)}, meta


def rollerLabel(i, flag):
    return {"d": "sprocket", "e": "idler"}.get(flag, "road wheel %d" % (i - 1))


def terrainZ(meta, x, y):
    if not meta["terrain"]:
        return None
    n, d = meta["terrain"]
    return (d - n[0] * x - n[1] * y) / n[2]


# --- figures -------------------------------------------------------------

def settling(c, meta, out, dark):
    import matplotlib.pyplot as plt
    col, primary, secondary, gridc = style(dark)
    t = c["t"]
    fig, ax = plt.subplots(2, 4, figsize=(17, 7))
    fig.suptitle("Settling run -- %s" % meta.get("model", ""), y=0.99)

    ax[0][0].plot(t, c["hull_z"], color=col[0])
    ax[0][0].set_title("Hull height")
    ax[0][0].set_ylabel("z (m)")

    for i, (k, lab) in enumerate([("hull_roll_deg", "roll"),
                                  ("hull_pitch_deg", "pitch"),
                                  ("hull_yaw_deg", "yaw")]):
        ax[0][1].plot(t, c[k], color=col[i], label=lab)
    ax[0][1].set_title("Hull attitude")
    ax[0][1].set_ylabel("deg")
    ax[0][1].legend()

    # Vertical and horizontal on separate axes, not one panel with two scales:
    # the vertical load is the vehicle's weight and the lateral one is a couple
    # of percent of it, so together they would be one line and a flat trace.
    # The weight is the number the settled load is asserted against, so it is
    # drawn rather than left to be inferred off the axis.
    w = float(meta.get("weight_N", 0))
    ax[0][2].plot(t, c["fgnd_z"], color=col[0])
    if w:
        ax[0][2].axhline(w, color=secondary, lw=1.0, ls="--")
        ax[0][2].annotate("weight %.4g N" % w, (t[-1], w), ha="right",
                          va="bottom", fontsize=8, color=secondary)
    ax[0][2].set_title("Ground reaction, vertical")
    ax[0][2].set_ylabel("N")

    ax[0][3].plot(t, c["fgnd_x"], color=col[0], label="fore/aft")
    ax[0][3].plot(t, c["fgnd_y"], color=col[1], label="lateral")
    ax[0][3].set_title("Ground reaction, horizontal")
    ax[0][3].set_ylabel("N")
    ax[0][3].legend()

    ax[1][0].plot(t, c["pen_band"], color=col[0], label="pin into roller")
    ax[1][0].plot(t, c["pen_gnd"], color=col[1], label="pad into terrain")
    ax[1][0].set_title("Worst penetration")
    ax[1][0].set_ylabel("m")
    ax[1][0].legend()

    ax[1][1].plot(t, c["n_band"], color=col[0], label="loaded shoe-roller pairs")
    ax[1][1].plot(t, c["n_gnd"], color=col[1], label="loaded pad corners")
    ax[1][1].set_title("Contacts carrying load")
    ax[1][1].set_ylabel("count")
    ax[1][1].legend()

    ax[1][2].plot(t, c["hull_x"] - c["hull_x"][0], color=col[0], label="x")
    ax[1][2].plot(t, c["hull_y"] - c["hull_y"][0], color=col[1], label="y")
    ax[1][2].set_title("Hull horizontal drift (creep)")
    ax[1][2].set_ylabel("m")
    ax[1][2].legend()

    # The gain ramp is context for every panel rather than a panel of its own:
    # before it closes the contact stiffness is still easing in and an
    # artificial viscosity is still draining velocity, so nothing before the
    # mark is the force law the rest of the run is showing.
    fig.delaxes(ax[1][3])
    ramp = t[np.argmax(c["gain"] >= 0.999)] if (c["gain"] >= 0.999).any() else 0
    if ramp > 0:
        for a in ax.ravel()[:-1]:
            a.axvline(ramp, color=gridc, lw=1.0)
        ax[0][0].annotate("gain ramp ends\n%.3g s" % ramp, (ramp, 1),
                          xycoords=("data", "axes fraction"), xytext=(5, -16),
                          textcoords="offset points", fontsize=8, color=secondary)

    for a in ax.ravel():
        a.set_xlabel("t (s)")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    return fig


def band(c, meta, out, dark):
    """The band itself, in the plane it lives in. Two samples per track: where
    the design pose put it, and where settling left it."""
    import matplotlib.pyplot as plt
    col, primary, secondary, gridc = style(dark)
    tracks = meta["tracks"]
    # Height per row is set by the band's own aspect ratio (about 5:1 over a
    # 6.5 in column) plus room for the legend under each panel: equal aspect on
    # a taller figure just pads it with whitespace.
    fig, ax = plt.subplots(len(tracks), 2, figsize=(13, 2.9 * len(tracks)),
                           squeeze=False)
    for r, tr in enumerate(tracks):
        for k, row in enumerate([0, len(c["t"]) - 1]):
            a = ax[r][k]
            sx = np.array([c["%s_shoe%d_x" % (tr["name"], i)][row]
                           for i in range(tr["shoes"])])
            sz = np.array([c["%s_shoe%d_z" % (tr["name"], i)][row]
                           for i in range(tr["shoes"])])
            a.plot(np.append(sx, sx[0]), np.append(sz, sz[0]), color=col[0],
                   lw=1.2, marker="o", ms=3, label="shoe pins")
            for i, (rad, flag) in enumerate(tr["radii"]):
                cx = c["%s_roller%d_x" % (tr["name"], i)][row]
                cz = c["%s_roller%d_z" % (tr["name"], i)][row]
                a.add_patch(plt.Circle((cx, cz), rad, fill=False, lw=1.6,
                                       color=col[1] if flag else col[2]))
            gy = c["%s_roller0_y" % tr["name"]][row]
            xs = np.linspace(sx.min() - 0.6, sx.max() + 0.6, 2)
            gz = terrainZ(meta, xs, gy)
            if gz is not None:
                a.plot(xs, gz, color=secondary, lw=1.2, ls="--")
                a.annotate("terrain", (xs[-1], gz[-1]), ha="right", va="top",
                           fontsize=8, color=secondary)
            # Two roller identities, both named in the legend rather than left
            # to the color.
            a.plot([], [], color=col[1], label="end pulleys")
            a.plot([], [], color=col[2], label="road wheels")
            a.set_aspect("equal")
            a.set_title("%s at t = %.3f s" % (tr["name"], c["t"][row]))
            a.set_xlabel("x (m)")
            a.set_ylabel("z (m)")
    # One legend for the figure: every panel draws the same three identities,
    # and at equal aspect the axes box is exactly the band, so an inset legend
    # would sit on top of it.
    h, l = ax[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out, dpi=140)
    return fig


def rollers(c, meta, out, dark):
    import matplotlib.pyplot as plt
    col, primary, secondary, gridc = style(dark)
    t = c["t"]
    tracks = meta["tracks"]
    fig, ax = plt.subplots(3, len(tracks), figsize=(6.8 * len(tracks), 10),
                           squeeze=False)
    for k, tr in enumerate(tracks):
        a = ax[0][k]
        for i, (rad, flag) in enumerate(tr["radii"]):
            a.plot(t, c["%s_roller%d_z" % (tr["name"], i)], color=col[i % len(col)],
                   label=rollerLabel(i, flag))
        a.set_title("%s -- roller heights" % tr["name"])
        a.set_ylabel("z (m)")
        a.set_xlabel("t (s)")
        a.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.18))

        # 63 shoes is never 63 colored series. What each shoe has DONE is the
        # readable statement -- absolute heights would only redraw the two runs
        # of the band, which the snapshot figure already shows -- so both rows
        # below are displacements from the design pose.
        for r, axis in enumerate("zx"):
            p = np.array([c["%s_shoe%d_%s" % (tr["name"], i, axis)]
                          for i in range(tr["shoes"])])
            d = p - p[:, :1]
            a = ax[1 + r][k]
            a.fill_between(t, d.min(axis=0), d.max(axis=0), color=col[0],
                           alpha=0.18, lw=0, label="all %d shoes" % tr["shoes"])
            a.plot(t, d.mean(axis=0), color=col[0], label="mean")
            a.set_title("%s -- shoe %s travel from the design pose"
                        % (tr["name"], "height" if axis == "z" else "fore/aft"))
            a.set_ylabel("d%s (m)" % axis)
            a.set_xlabel("t (s)")
            a.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    return fig


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="trace written by sim_<model>.exe --trace")
    ap.add_argument("--out", help="output directory (default: beside the CSV)")
    ap.add_argument("--dark", action="store_true", help="dark surface")
    ap.add_argument("--show", action="store_true", help="open the figures too")
    args = ap.parse_args()

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    c, meta = load(args.csv)
    stem = os.path.splitext(os.path.basename(args.csv))[0]
    outdir = args.out or os.path.dirname(os.path.abspath(args.csv))
    made = []
    for name, fn in [("settling", settling), ("band", band), ("rollers", rollers)]:
        if name != "settling" and not meta["tracks"]:
            continue
        p = os.path.join(outdir, "%s_%s.png" % (stem, name))
        fn(c, meta, p, args.dark)
        made.append(p)
    print("%d rows, %.3f s of run" % (len(c["t"]), c["t"][-1] - c["t"][0]))
    print("\n".join(made))
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
