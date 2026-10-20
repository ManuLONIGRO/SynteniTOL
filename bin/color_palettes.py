#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Categorical color palette.

Uses the 32-color Glasbey palette (Glasbey et al. 2007) for up to 32
categories.  For larger sets the palette is extended algorithmically
following the same Glasbey criterion: candidate colors are taken from a
fixed RGB grid and the one maximizing the minimum CIE94 color difference
(dE) to the colors already chosen is selected.  Candidates that would be
unreadable on iTOL's white background are excluded (too light, or pale: low
contrast against white combined with low chroma).  Selection is
deterministic, so the same n always yields the same colors.

Reference:
    Glasbey CA, van der Heijden GWAM, Toh VFK, Gray AJ (2007).
    Colour Displays for Categorical Images. Color Research and Application,
    32, 304-309.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np

# 32-color Glasbey palette (sourced from the R 'pals' package).
# Optimised for maximum perceptual distinctiveness in CIE L*u*v* space.
GLASBEY_PALETTE: List[str] = [
    "#0000FF", "#FF0000", "#00FF00", "#000033", "#FF00B6", "#005300",
    "#FFD300", "#009FFF", "#9A4D42", "#00FFBE", "#783FC1", "#1F9698",
    "#FFACFD", "#B1CC71", "#F1085C", "#FE8F42", "#DD00FF", "#201A01",
    "#720055", "#766C95", "#02AD24", "#C8FF00", "#886C00", "#FFB79F",
    "#858567", "#A10300", "#14F9FF", "#00479E", "#DC5E93", "#93D4FF",
    "#004CFF", "#F2F318",
]

# Candidate grid: every 8th level plus 255 -> 33 levels -> 33**3 candidates.
_GRID_STEP = 8
# Visibility on iTOL's white background.  A candidate must
#   - fall in the luminosity band [MIN_LIGHTNESS, MAX_LIGHTNESS], and
#   - not be "pale": low contrast against white *and* low chroma at the same
#     time (light but strongly saturated colors are kept: they separate by hue,
#     light grays/pastels are rejected because they vanish on the white plot).
_MIN_LIGHTNESS = 15.0
_MAX_LIGHTNESS = 85.0
_MIN_WHITE_CONTRAST = 2.5   # WCAG contrast ratio against white
_MIN_CHROMA = 25.0          # CIELAB chroma

# Cache of previously generated palettes: length -> colors.
_palette_cache: Dict[int, List[str]] = {}
# Candidate grid, shared by every extension call.
_candidates: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None


def _hex_to_rgb(hex_colors: List[str]) -> np.ndarray:
    """'#RRGGBB' strings -> (n, 3) uint8 array."""
    return np.array(
        [[int(h[1:3], 16), int(h[3:5], 16), int(h[5:7], 16)] for h in hex_colors],
        dtype=np.uint8,
    )


def _srgb_to_lab(rgb: np.ndarray) -> np.ndarray:
    """(n, 3) sRGB values (0-255) -> (n, 3) CIELAB, D65 reference white."""
    x = rgb.astype(np.float64) / 255.0
    lin = np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)
    r, g, b = lin[:, 0], lin[:, 1], lin[:, 2]

    # sRGB (D65) -> XYZ
    X = 0.4124564 * r + 0.3575761 * g + 0.1804375 * b
    Y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b
    Z = 0.0193339 * r + 0.1191920 * g + 0.9503041 * b

    xn, yn, zn = 0.95047, 1.0, 1.08883
    epsilon = 216.0 / 24389.0
    kappa = 24389.0 / 27.0

    def _f(t: np.ndarray, t0: float) -> np.ndarray:
        ratio = t / t0
        return np.where(ratio > epsilon, np.cbrt(ratio), (kappa * ratio + 16.0) / 116.0)

    fx = _f(X, xn)
    fy = _f(Y, yn)
    fz = _f(Z, zn)

    L = 116.0 * fy - 16.0
    a = 500.0 * (fx - fy)
    b_lab = 200.0 * (fy - fz)
    return np.stack([L, a, b_lab], axis=1)


def _delta_e94(lab_candidates: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """CIE94 color difference (graphic arts constants) to a single reference."""
    L1, a1, b1 = lab_candidates[:, 0], lab_candidates[:, 1], lab_candidates[:, 2]
    L2, a2, b2 = reference[0], reference[1], reference[2]

    dL = L1 - L2
    C1 = np.hypot(a1, b1)
    C2 = np.hypot(a2, b2)
    dC = C1 - C2
    da = a1 - a2
    db = b1 - b2
    # dH^2 = da^2 + db^2 - dC^2 (avoids computing hue angles)
    dH2 = np.maximum(da * da + db * db - dC * dC, 0.0)

    sC = 1.0 + 0.045 * C2  # Cb = chroma of the reference color
    sH = 1.0 + 0.015 * C2
    return np.sqrt(dL * dL + (dC / sC) ** 2 + dH2 / (sH * sH))


def _relative_luminance(rgb: np.ndarray) -> np.ndarray:
    """(n, 3) sRGB values (0-255) -> WCAG relative luminance."""
    x = rgb.astype(np.float64) / 255.0
    lin = np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)
    return 0.2126 * lin[:, 0] + 0.7152 * lin[:, 1] + 0.0722 * lin[:, 2]


def _visible_mask(rgb: np.ndarray, lab: np.ndarray) -> np.ndarray:
    """Which colors remain readable on iTOL's white background (see constants)."""
    contrast = 1.05 / (_relative_luminance(rgb) + 0.05)
    chroma = np.hypot(lab[:, 1], lab[:, 2])
    pale = (contrast < _MIN_WHITE_CONTRAST) & (chroma < _MIN_CHROMA)
    lightness = lab[:, 0]
    return (lightness >= _MIN_LIGHTNESS) & (lightness <= _MAX_LIGHTNESS) & ~pale


def _candidate_grid() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (rgb, lab, visible) for the fixed candidate grid (cached)."""
    global _candidates
    if _candidates is None:
        levels = np.unique(
            np.append(np.arange(0, 256, _GRID_STEP, dtype=np.int16), 255)
        )
        mesh = np.stack(
            np.meshgrid(levels, levels, levels, indexing="ij"), axis=-1
        ).reshape(-1, 3)
        rgb = mesh.astype(np.uint8)
        lab = _srgb_to_lab(rgb)
        visible = _visible_mask(rgb, lab)
        _candidates = (rgb, lab, visible)
    return _candidates


def _glasbey_extend(base_colors: List[str], n: int) -> List[str]:
    """Extend *base_colors* to *n* colors with the Glasbey maximin criterion."""
    rgb, lab, visible = _candidate_grid()

    used = np.zeros(len(rgb), dtype=bool)
    min_dist = np.full(len(rgb), np.inf)
    for hex_color in base_colors:
        ref = _srgb_to_lab(_hex_to_rgb([hex_color]))[0]
        min_dist = np.minimum(min_dist, _delta_e94(lab, ref))

    colors = list(base_colors)
    while len(colors) < n:
        scores = np.where(visible & ~used, min_dist, -1.0)
        idx = int(np.argmax(scores))
        if scores[idx] <= 0.0:  # no candidate left (practically unreachable)
            break
        colors.append("#%02X%02X%02X" % tuple(int(v) for v in rgb[idx]))
        used[idx] = True
        min_dist = np.minimum(min_dist, _delta_e94(lab, lab[idx]))
    return colors


def categorical_colors(n: int) -> List[str]:
    """Return *n* maximally-distinct categorical colors.

    For n <= 32 the hardcoded Glasbey palette is used.  Larger n extends it
    algorithmically (Glasbey maximin selection in CIELAB with CIE94 dE),
    deterministically and with caching between calls.
    """
    if n <= 0:
        return []
    if n <= len(GLASBEY_PALETTE):
        return list(GLASBEY_PALETTE[:n])

    for cached_n in sorted(_palette_cache):
        if cached_n >= n:
            return list(_palette_cache[cached_n][:n])

    # Seed with the longest palette already generated (or the hardcoded one).
    seed = GLASBEY_PALETTE
    for cached_n in sorted(_palette_cache):
        if cached_n < n:
            seed = _palette_cache[cached_n]
    colors = _glasbey_extend(seed, n)
    _palette_cache[n] = colors
    return list(colors)
