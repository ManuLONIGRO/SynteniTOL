#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Color palette for categorical plots.

Uses the 32-color Glasbey palette (Glasbey et al. 2007) for up to 32
categories. Falls back to an evenly-spaced HLS rainbow for larger sets.

Reference:
    Glasbey CA, van der Heijden GWAM, Toh VFK, Gray AJ (2007).
    Colour Displays for Categorical Images. Color Research and Application,
    32, 304-309.
"""

import colorsys
from typing import List

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


def rainbow_colors(n: int) -> List[str]:
    """Evenly-spaced HLS rainbow palette (fallback for n > 32)."""
    colors: List[str] = []
    for i in range(n):
        h = (i / max(1, n)) % 1.0
        r, g, b = colorsys.hls_to_rgb(h, 0.55, 0.95)
        colors.append(f"#{int(r * 255):02X}{int(g * 255):02X}{int(b * 255):02X}")
    return colors


def categorical_colors(n: int) -> List[str]:
    """Return *n* maximally-distinct categorical colors.

    For n <= 32 the Glasbey palette is used.  For larger n the function
    falls back to an evenly-spaced HLS rainbow (colors will repeat).
    """
    if n <= len(GLASBEY_PALETTE):
        return list(GLASBEY_PALETTE[:n])
    return rainbow_colors(n)
