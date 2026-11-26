"""
@file make_figure_pptx.py
@brief The paper's figure list as a plain PowerPoint file (2026-11-26): the default Office template, white slides, the default title,
       one figure per slide and a plain text box of caption notes. Figures come from analysis/fig_npj.py.

    <python with python-pptx> analysis/make_figure_pptx.py [figs_dir] [out.pptx]
"""
from __future__ import annotations

import os
import sys

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.dml import MSO_LINE_DASH_STYLE
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(os.path.dirname(HERE), "ptychoshelves-clean", "report", "aberrations_paper")

# (title, figure file or None, to-make text or None, notes)
SLIDES = [
    ("Fig. 1  The atom against aperture", "npj_atoms_aperture", None, [
        "a: perfect probe, aperture 30, 40, 60 and 80 mrad",
        "b: 80 mrad aperture over a corrector flat only to 30, 40, 60 and 80 mrad (dashed circle); the phase across the aperture above, "
        "the atom below; each column has the same x as panel a above it. Flat to 80 is the perfect 80 mrad probe: its phase is flat and "
        "its atom is the last one in a",
        "x–z through each atom's peak, scaled to its own peak; probe refined in two passes (C1 and A1)",
        "W80: the depth holding 80 % of the atom's phase. Peak: the atom's highest phase",
        "The perfect 30 and 40 mrad peaks are low (0.51, 0.63 rad) because the atom's total phase, 0.97 and 1.07 rad·Å, is spread over "
        "9.1 and 7.3 Å; flat to 30 and 40 opened to 80 hold 1.09 and 1.16 rad·Å over 5.5 and 2.6 Å (peaks 1.02, 1.43 rad)"]),
    ("Fig. 2  Depth extent against x", "npj_w80_vs_x", None, [
        "W80 against x for both lines; solid: probe refined in two passes; dotted: probe known",
        "x = 80 is the same probe on both lines",
        "Flat to 40 and 60 opened to 80: 2.57 and 2.57 Å refined, within 6 % of the probe known; perfect 80: 2.43 Å",
        "Flat to 30 opened to 80: 5.47 Å refined against 3.56 Å with the probe known (54 % longer)"]),
    ("Fig. 3  Depth extent against dose", "npj_dose_line", None, [
        "Poisson likelihood on the full patterns (712 pixels; 914 for flat to 30)",
        "Solid: probe known; dashed: A1 refined, C1 held at the search's value; ×: no atom or broken",
        "Perfect 80 keeps its atom to 10⁶ e Å⁻² with the probe known (4.13 Å); the refinement fails there",
        "Flat to 30 at 80 is lost at 10⁶; perfect 30 diverges at all three doses on 1 Å slices (Supp. Fig. 10)",
        "To do: the refined line with the two-pass C1 + A1 recipe"]),
    ("Fig. 4  The atoms at each dose", "npj_dose_atoms", None, [
        "The atoms behind Fig. 3: x–z through the peak, each scaled to its own peak",
        "Rows: perfect 80 and flat to 30 at 80, each with the probe known and with A1 refined",
        "No atom: noise or streaks along the column. Broken: an X-shaped artefact in place of the atom"]),
    ("Fig. 5  Focal spread", "npj_focal_spread", None, [
        "W80 against x at an rms focal spread of 0, 2.5, 5 and 7.29 Å (0.10, 0.21 and 0.30 eV FWHM for Cc = 1.4 mm)",
        "One coherent probe mode, started from the search's C1 and C3 errors with the true A1",
        "From 5 Å, flat to 30 or 40 opened to 80 is longer than before opening (13.51 and 12.75 Å against 8.47 and 6.49 Å at 5 Å)",
        "At 7.29 Å the perfect 80 atom is longer than the perfect 60 one (7.59 against 6.25 Å)",
        "To do: 8 probe modes (Supp. Fig. 12); source size on these probes"]),
    ("Fig. 6  How far a flat-to-30 corrector can be opened", None,
     "W80 against aperture, 40 to 80 mrad: the perfect probe at x beside the flat-to-30 corrector opened to x. No noise, and Poisson at "
     "10⁷ and 10⁸ e Å⁻². Inset: the share of the probe outside the solver's window.", [
        "Runs: gopen30 (P23) submitted; the simulations are waiting for slots",
        "Prediction: more than 1 Å from the perfect probe at 60–70 mrad under noise",
        "90 and 100 mrad not run: they need windows of 80 Å and over 105 Å"]),
    ("Fig. 7  The crystal", None,
     "PTO/STO superlattice: found atoms drawn over the truth (x–z through a Pb column and a Ti–O column) for flat to 30 used at 30, opened "
     "to 50 and 80, and the perfect 80 mrad probe. Below: the share of Pb and Ti found in the right place.", [
        "Score: the right element within 0.6 Å in plane and 0.5 Å in depth (Supp. Note 7)",
        "With the coherence wall and a monochromated beam",
        "Runs: to set up"]),
    ("Fig. 8  Vibration and cooling", None,
     "The crystal at room temperature, liquid nitrogen and liquid helium, probe opened to 80 mrad: found atoms over the truth, and Pb and "
     "Ti found in the right place against temperature.", [
        "Frozen phonons, independent (Einstein) displacements",
        "Needs new runs after the √3 phonon fix (new simulation names)"]),
    ("Fig. 9  Today's microscopes", None,
     "Instrument classes, by their numbers only (energy spread, chromatic aberration, source size), placed on the focal-spread and source "
     "results, with the aperture each class can use.", [
        "No brand names: numbers only",
        "Needs Fig. 5 and Supp. Figs 5 and 6 finished: probe modes, source size on the flat probes"]),
    ("Supp. Fig. 1  Ronchigrams of the flat probes", None,
     "Ronchigrams of the probes flat to 30, 40 and 60 mrad opened to 80 mrad (amorphous carbon, abTEM), beside the phase across each "
     "aperture.", [
        "The Results text points to a Ronchigram figure",
        "Perfect probes left out: their Ronchigram is a plain disc",
        "Needs a short simulation"]),
    ("Supp. Fig. 2  The refinement's gain", None,
     "The perfect 80 mrad probe's phase error across the aperture, at the start and after 60 + 60 iterations at gains of 10, 30 and 100.", [
        "A1 error left: 2.22, 0.12 and 0.00 Å of 3.6 Å",
        "A draft exists (fig_probe.py), not yet in npj style",
        "Could be a sentence in the SI instead"]),
    ("Supp. Fig. 3  The search's A1 error: the atoms", "npj_a1_atoms", None, [
        "a: start probe with the search's C1 and C3 errors and the true A1; b: A1 also 3.6 Å off; probe not refined",
        "Waves: that A1 error at each probe's own aperture edge",
        "Every probe where it reaches 1/3 wave or more breaks: perfect 60 (0.33 waves), perfect 80 and the three flat probes (0.59 waves)",
        "Perfect 30 and 40 (0.08 and 0.15 waves) keep one atom; W80 9.00 to 7.19 Å and 7.08 to 5.38 Å, peak 0.48 to 0.34 and 0.61 to 0.69 rad"]),
    ("Supp. Fig. 4  The A1 ladder", "npj_a1_ladder", None, [
        "A1 error of 0, 0.25, 0.5, 1.0 and 3.6 Å in the start probe; the search's C1 and C3 errors kept; probe not refined",
        "a: W80; b: peak phase; perfect 80 mrad and flat to 30 opened to 80",
        "Dotted: the π/4 criterion for twofold astigmatism at the 80 mrad edge, λ/4α² = 0.77 Å",
        "W80 hides the damage to perfect 80 at 3.6 Å (3.05 Å); its peak falls from 1.80 to 0.66 rad"]),
    ("Supp. Fig. 5  Source size", "npj_source", None, [
        "Perfect 80 mrad, true probe held, no noise; source 0, 0.05, 0.1, 0.2 and 0.4 Å FWHM",
        "a: W80 2.51, 2.57, 3.15, 3.72 and 4.21 Å; b: peak 1.86, 1.76, 1.27, 0.38 and 0.23 rad",
        "The atom stays single but spreads sideways and fades (Supp. Fig. 6)",
        "16 random position shifts leave the beam more coherent than the stated source, so these are optimistic",
        "To do: source size on the flat probes"]),
    ("Supp. Fig. 6  Source size: the atoms", "npj_source_atoms", None, [
        "Perfect 80 mrad, true probe held, no noise: source 0, 0.1, 0.2 and 0.4 Å FWHM",
        "x–z through the peak, each scaled to its own peak",
        "The atom spreads sideways and fades as the source grows"]),
    ("Supp. Fig. 7  Flat to 30: the faint tail after refinement", "npj_ghost", None, [
        "a: depth profiles (core phase) with the probe known, refined in two passes (C1 and A1), and A1 refined alone",
        "b: the three atoms",
        "Two passes: the core matches the probe known, plus a faint tail 3–8 Å away; A1 alone: no far tail, peak 0.58 against 1.15 rad",
        "To do: true probe held with 2.0 and 2.5 Å of defocus added (P26), for the tail's cause"]),
    ("Supp. Fig. 8  Camera pixels", "npj_camera", None, [
        "No noise, true probe held; patterns summed 4 × 4 into a camera",
        "Perfect 80: 2.51 Å at 712 px, 2.27 Å at 178 px. Flat to 30: 3.56 Å at 914 px, 4.78 Å at 228 px",
        "The aberrated probe's fine pattern features need fine pixels (Nguyen et al. 2024)"]),
    ("Supp. Fig. 9  W80 from the histogram, not a KDE", "npj_kde", None, [
        "Change in W80 when the depth profile is smoothed with a Gaussian KDE of bandwidth 0.25, 0.5 and 1.0 Å; true-probe atoms",
        "At 0.5 Å: perfect 80 +18 %, perfect 30 +3 %",
        "A KDE adds its own width and squeezes the differences measured; the histogram matches what a slice is in the model"]),
    ("Supp. Fig. 10  Perfect 30 mrad under noise", "npj_p30noise", None, [
        "1 Å slices, probe known, Poisson: the phase runs to ±π at 10⁸, 10⁷ and 10⁶ e Å⁻²",
        "2 Å slices, no noise: W80 11.25 Å and peak 0.90 rad, against 8.42 Å and 0.46 rad on 1 Å slices",
        "To do: 2 Å slices under noise (P25)"]),
    ("Supp. Fig. 11  Higher-order aberrations the fit leaves", None,
     "Flat to 40 at 80 mrad, no noise, C1 and A1 fitted: the control, A2 off by 1/8 and 1/2 wave, A3 and S3 off by 1/2 wave, and A2 off by "
     "1/2 wave with A2 fitted. The atom and W80 for each.", [
        "Runs: P24, 3 of 6 pulled, not yet scored",
        "Tests whether the 1/8-wave rule holds beyond A1",
        "Full test later: tableau-sized errors in every term"]),
    ("Supp. Fig. 12  Probe modes under partial coherence", None,
     "W80 against x at 5 and 7.29 Å of focal spread with 8 probe modes, beside the single-mode result of Fig. 5, and the occupancies the "
     "modes reached.", [
        "Runs: P22; the 20-slice points are done, the 40-slice points are rerunning at 140 + 140 iterations",
        "Methods must say the iteration counts differ"]),
    ("Supp. Fig. 13  Finding atoms in the crystal", None,
     "The finder on one column: preprocessing, the column seed, matching pursuit placing single-atom responses, Gauss–Newton refinement, "
     "element assignment; then the found atoms with their 95 % depth intervals over the truth.", [
        "Goes with Supp. Note 7",
        "Tools exist: analysis/atomfind"]),
]

W, H = 13.333, 7.5
L, R, TOP, BOT = 0.5, 12.83, 1.35, 7.15


def set_title(slide, text):
    t = slide.shapes.title
    t.text = text
    t.left, t.top, t.width, t.height = Inches(L), Inches(0.3), Inches(R - L), Inches(0.9)
    p = t.text_frame.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    for r in p.runs:
        r.font.size = Pt(30)


def notes_box(slide, x, y, w, h, items, size=16, head=None, head_size=None):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    first = True
    if head:
        p = tf.paragraphs[0]; first = False
        r = p.add_run(); r.text = head; r.font.size = Pt(head_size or size); r.font.bold = True
        p.space_after = Pt(10)
    for it in items:
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        r = p.add_run(); r.text = it; r.font.size = Pt(size)
        p.space_after = Pt(6)
        pPr = p._p.get_or_add_pPr()
        pPr.set("marL", str(int(Inches(0.22)))); pPr.set("indent", str(-int(Inches(0.22))))
        pPr.append(pPr.makeelement(qn("a:buChar"), {"char": "•"}))
    return tb


def plain_box(slide, x, y, w, h, head, text):
    box = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    box.fill.background()
    box.line.color.rgb = RGBColor(0xA6, 0xA6, 0xA6); box.line.width = Pt(1); box.line.dash_style = MSO_LINE_DASH_STYLE.DASH
    box.text_frame.text = ""
    tb = slide.shapes.add_textbox(Inches(x + 0.4), Inches(y + 0.4), Inches(w - 0.8), Inches(h - 0.8))
    tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; r = p.add_run(); r.text = head; r.font.size = Pt(20); r.font.bold = True; p.space_after = Pt(12)
    p = tf.add_paragraph(); r = p.add_run(); r.text = text; r.font.size = Pt(18)


def main(figs=None, out=None):
    figs = figs or os.path.join(PAPER, "figs", "npj")
    out = out or os.path.join(PAPER, "FigPPT.pptx")
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    s = prs.slides.add_slide(prs.slide_layouts[0])
    s.shapes.title.text = "Figures for the depth paper"
    s.shapes.title.left, s.shapes.title.top, s.shapes.title.width, s.shapes.title.height = Inches(L), Inches(2.2), Inches(R - L), Inches(1.4)
    sub = s.placeholders[1]
    sub.left, sub.top, sub.width, sub.height = Inches(L), Inches(3.8), Inches(R - L), Inches(1.6)
    sub.text = "Working figure list, 26 November 2026. One figure per slide: a first attempt where the data exist, a description where they don't."
    for slide_def in SLIDES:
        title, fig, todo, notes = slide_def
        s = prs.slides.add_slide(prs.slide_layouts[5])
        set_title(s, title)
        if fig:
            png = os.path.join(figs, fig + ".png")
            with Image.open(png) as im:
                a = im.width / im.height
            if a >= 2.0:                                   # wide: figure across the top, notes below
                w = R - L; h = w / a
                if h > 3.6:
                    h = 3.6; w = h * a
                pic = s.shapes.add_picture(png, Inches(L + (R - L - w) / 2), Inches(TOP), Inches(w), Inches(h))
                notes_box(s, L, TOP + h + 0.25, R - L, BOT - (TOP + h + 0.25), notes, size=15)
            else:                                          # figure left, notes right
                maxw, maxh = 8.0, BOT - TOP
                w = min(maxw, maxh * a); h = w / a
                pic = s.shapes.add_picture(png, Inches(L), Inches(TOP), Inches(w), Inches(h))
                x = L + w + 0.35
                notes_box(s, x, TOP, R - x, BOT - TOP, notes, size=15)
            pic._element.nvPicPr.cNvPr.set("descr", title)
        else:
            plain_box(s, L, TOP, 7.6, BOT - TOP, "To make", todo)
            notes_box(s, L + 7.95, TOP, R - (L + 7.95), BOT - TOP, notes, size=15)
    prs.save(out)
    print("wrote", out)


if __name__ == "__main__":
    main(*sys.argv[1:])
