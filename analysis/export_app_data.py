"""
@file export_app_data.py
@brief Writes results/app_data.json, the data behind the blob explorer page: every 80 mrad lone-atom kernel (both probes),
cut down to what the page draws.

Per kernel: the centre-column depth profile, the same for the clean kernel blurred as the physics says (deconv_test's
convolution model, not aligned, so a depth shift shows), the x-z cut through the centre (+-1.5 A, every second pixel), the
layer through the clean atom's peak (+-1.2 A, every second pixel), the tilt curve C(theta), and the blob numbers (peak,
contrast, depth FWHM, sideways width S, share on the next atom's site C0, theta_10 / theta_5, the share the convolution
model leaves unexplained). Depth is relative to the clean atom's peak on the same probe. Run after blobs.py and
deconv_test.py (it reads their TSVs).

    ~/hyperspy-bundle/bin/python analysis/export_app_data.py ~/Desktop/zblur/g0d ~/Desktop/zblur/g1 ~/Desktop/zblur/g1b
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import blobs  # noqa: E402
import deconv_test  # noqa: E402

HALF_XZ, HALF_XY, EVERY = 1.5, 1.2, 2


def tsv(name):
    with open(os.path.join(blobs.HERE, "results", name)) as fh:
        return {r["name"].replace("  [NOT CLEAN]", ""): r for r in csv.DictReader(fh, delimiter="\t")}


def num(v, nd=4):
    v = float(v)
    return None if not np.isfinite(v) else round(v, nd)


def main(*packs):
    tag = "_".join(os.path.basename(os.path.normpath(p)) for p in packs)
    bt, dt = tsv(f"blobs_{tag}.tsv"), tsv(f"deconv_{tag}.tsv")
    rs = [r for r in blobs.load(packs) if float(r["alpha"]) == 80 and str(r["nl"]) == "20"]
    bases = {r["aber"]: r for r in rs if not r["fs"] and not r["src"] and r["T"] is None}
    out = dict(spacing=blobs.SPACING, thetas=[round(float(t), 5) for t in blobs.THETAS], probes={}, kernels=[])
    for aber, b in bases.items():
        cb = blobs.core(b)
        nz, ny, nx = b["K"].shape
        out["probes"][aber] = dict(z=[round(b["z0"] + (k + 0.5) * b["dz"] - cb["z_peak"], 3) for k in range(nz)],
                                   peak=round(cb["peak"], 4), k_peak=cb["k_peak"], dx=b["dx"], dz=b["dz"],
                                   x_xz=round(HALF_XZ, 3), x_xy=round(HALF_XY, 3))
    for r in rs:
        b = bases[r["aber"]]
        name = blobs.describe(r)
        c = blobs.core(r)
        K = deconv_test.bg(r["K"])
        nz, ny, nx = K.shape
        wx, wy = int(round(HALF_XZ / r["dx"])), int(round(HALF_XY / r["dx"]))
        xz = K[:, ny // 2, nx // 2 - wx:nx // 2 + wx + 1:EVERY]
        kp = out["probes"][r["aber"]]["k_peak"]
        xy = K[kp, ny // 2 - wy:ny // 2 + wy + 1:EVERY, nx // 2 - wy:nx // 2 + wy + 1:EVERY]
        row = bt.get(name, {})
        item = dict(aber=r["aber"], fs=r["fs"], src=r["src"], T=int(r["T"]) if r["T"] is not None else 0,
                    profile=[round(float(v), 4) for v in blobs.zp.centre_profile(r["K"])],
                    xz=np.round(xz, 3).tolist(), xy=np.round(xy, 3).tolist(),
                    tilt=[round(float(v), 4) for v in blobs.tilt_curve(r, c)],
                    peak=round(c["peak"], 4), fwhm_z=num(c["fwhm_z"], 3), dip=num(c["dip"], 3), side=num(c["side"], 3),
                    contrast=round(c["peak"] / blobs.core(b)["peak"], 4), S=num(row.get("S", "nan")),
                    C0=num(row.get("C0", "nan")), theta10=num(row.get("theta10", "nan")),
                    theta5=num(row.get("theta5", "nan")), clean=bool(np.isfinite(c["fwhm_z"]) and c["dip"] > -0.9))
        if r is not b:
            item["pred"] = [round(float(v), 4) for v in blobs.zp.centre_profile(deconv_test.predict(b, r))]
            d = dt.get(name)
            item["unexplained"] = num(d["resid_conv"], 3) if d else None
        out["kernels"].append(item)
    path = os.path.join(blobs.HERE, "results", "app_data.json")
    with open(path, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    print(f"wrote {path}: {len(out['kernels'])} kernels, {os.path.getsize(path) / 1e3:.0f} kB")


if __name__ == "__main__":
    main(*(sys.argv[1:] or [os.path.expanduser(f"~/Desktop/zblur/{p}") for p in ("g0d", "g1", "g1b")]))
