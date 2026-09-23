"""Tiny shim providing only what gel_annotator.py needs (binary_fill_holes),
implemented with OpenCV since real scipy isn't installable in this sandbox."""
import numpy as np
import cv2


def binary_fill_holes(input):
    arr = np.asarray(input, dtype=np.uint8)
    h, w = arr.shape
    inv = ((arr == 0).astype(np.uint8)) * 255
    floodfilled = inv.copy()
    mask = np.zeros((h + 2, w + 2), np.uint8)
    cv2.floodFill(floodfilled, mask, (0, 0), 0)
    holes = floodfilled > 0
    return (arr.astype(bool) | holes)
