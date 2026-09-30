"""Image augmentations for optional patch preprocessing.

Images use channels-last arrays with intensities in [0, 1]. HSV operations
expect BGR channel order. Each callable applies its specified random transform.
"""

import random

import cv2
import numpy as np
from PIL import Image


def do_random_revolve(image, s=0.5):
    """Apply a random resized crop and horizontal flip to a square output.

    ``s`` is accepted for caller compatibility and does not affect the transform.
    """
    from torchvision.transforms import RandomHorizontalFlip, RandomResizedCrop

    size = image.shape[1]
    image_uint8 = np.uint8(np.clip(image, 0, 1) * 255)
    cropped = RandomResizedCrop(size=size)(Image.fromarray(image_uint8))
    flipped = RandomHorizontalFlip()(cropped)
    return np.asarray(flipped, dtype=np.float32) / 255.0


def do_random_flip(image):
    """Independently flip rows, columns, and spatial axes with probability 0.5."""
    if np.random.rand() > 0.5:
        image = cv2.flip(image, 0)
    if np.random.rand() > 0.5:
        image = cv2.flip(image, 1)
    if np.random.rand() > 0.5:
        image = image.transpose(1, 0, 2)
    return np.ascontiguousarray(image)


def do_random_rot90(image):
    """Uniformly sample a rotation by zero, one, two, or three quarter turns."""
    # NumPy quarter turns avoid overlap between OpenCV's zero code and identity.
    return np.ascontiguousarray(np.rot90(image, k=int(np.random.randint(4))))


def do_random_contrast(image, mag=0.3):
    """Scale intensity by a random factor centered on one and clip its range."""
    factor = 1 + random.uniform(-1, 1) * mag
    return np.clip(image * factor, 0, 1)


# The compatibility alias accepts the established public spelling.
do_random_contast = do_random_contrast


def do_random_hsv(image, mag=(0.15, 0.25, 0.25)):
    """Perturb hue, saturation, and value of a BGR image independently."""
    image_uint8 = (np.clip(image, 0, 1) * 255).astype(np.uint8)
    hsv = cv2.cvtColor(image_uint8, cv2.COLOR_BGR2HSV)
    hue, saturation, value = np.moveaxis(hsv.astype(np.float32), -1, 0)
    hue = (hue * (1 + random.uniform(-1, 1) * mag[0])) % 180
    saturation *= 1 + random.uniform(-1, 1) * mag[1]
    value *= 1 + random.uniform(-1, 1) * mag[2]
    hsv[:, :, 0] = np.clip(hue, 0, 179).astype(np.uint8)
    hsv[:, :, 1] = np.clip(saturation, 0, 255).astype(np.uint8)
    hsv[:, :, 2] = np.clip(value, 0, 255).astype(np.uint8)
    return cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR).astype(np.float32) / 255.0


def do_random_noise(image, mag=0.1):
    """Add spatially independent uniform noise shared across color channels."""
    height, width = image.shape[:2]
    noise = np.random.uniform(-1, 1, (height, width, 1)) * mag
    return np.clip(image + noise, 0, 1)


def do_random_rotate_scale(image, angle=35, scale=(0.6, 1.4)):
    """Rotate and scale around the image center with zero-padded boundaries."""
    sampled_angle = np.random.uniform(-angle, angle)
    sampled_scale = np.random.uniform(*scale) if scale is not None else 1
    height, width = image.shape[:2]
    center = (width / 2.0, height / 2.0)
    transform = cv2.getRotationMatrix2D(center, sampled_angle, sampled_scale)
    return cv2.warpAffine(image, transform, (width, height), flags=cv2.INTER_LINEAR,
                          borderMode=cv2.BORDER_CONSTANT, borderValue=(0, 0, 0))
