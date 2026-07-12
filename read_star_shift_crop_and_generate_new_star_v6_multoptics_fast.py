#!/usr/bin/env python3
"""
Read a RELION-style STAR file, shift/crop particle images, write new MRC stacks,
and generate a matching STAR file.

This version adds:
  * Multi-optics-group support: pixel size and original box size are selected
    independently for every particle via rlnOpticsGroup.
  * Fast source I/O: source MRC files are memory-mapped through an LRU cache,
    and particles are processed by optics group and source file for locality.
  * Robust MRC layout handling, including:
      - standard 3-D stacks:              (N, box, box)
      - a single 2-D MRC image:           (box, box)
      - flattened-row stacks:             (N, box*box)
      - vertically concatenated stacks:   (N*box, box)
      - horizontally concatenated stacks: (box, N*box)
      - fully flattened data:             (N*box*box,)
  * Per-output-stack voxel size and image-stack header metadata.
  * Exact center cropping for both even and odd crop sizes.

Dependencies:
    pip install numpy mrcfile

Example:
    python read_star_shift_crop_and_generate_new_star_v6_multoptics_fast.py \
        --star_name particles.star \
        --output_root_name cropped \
        --newboxsize 128 \
        --batchsize 5000
"""

import argparse
import os
import shlex
import sys
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import mrcfile
import numpy as np


PROGRAM_VERSION = "6.0"


# -----------------------------------------------------------------------------
# Lightweight STAR reader/writer
# -----------------------------------------------------------------------------


def _split_star_line(line: str) -> List[str]:
    """Split one ordinary STAR line while respecting single/double quotes."""
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return []
    if stripped.startswith(";"):
        raise ValueError(
            "Semicolon-delimited multiline STAR values are not supported by this "
            "lightweight parser. RELION optics/particles tables normally do not use them."
        )
    lexer = shlex.shlex(stripped, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    return list(lexer)


def _normalise_block_name(name: str) -> str:
    name = name.strip()
    if name.lower().startswith("data_"):
        name = name[5:]
    return name.lower()


def _normalise_label(label: str) -> str:
    return label.strip().lstrip("_").lower()


def read_star_text(filename: str) -> dict:
    """Read a RELION-style STAR file without pandas/starfile."""
    with open(filename, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    star = {"preamble": [], "blocks": OrderedDict()}
    current_block = None
    i = 0

    def ensure_block(block_name: str) -> None:
        if block_name not in star["blocks"]:
            star["blocks"][block_name] = {"items": [], "loops": []}

    while i < len(lines):
        raw = lines[i]
        stripped = raw.strip()
        lower = stripped.lower()

        if lower.startswith("data_"):
            current_block = stripped[5:].strip()
            ensure_block(current_block)
            i += 1
            continue

        if current_block is None:
            star["preamble"].append(raw.rstrip("\n"))
            i += 1
            continue

        if not stripped or stripped.startswith("#"):
            i += 1
            continue

        if lower == "loop_":
            i += 1
            columns: List[str] = []

            while i < len(lines):
                header_line = lines[i].strip()
                if not header_line or header_line.startswith("#"):
                    i += 1
                    continue
                if not header_line.startswith("_"):
                    break
                tokens = _split_star_line(lines[i])
                if tokens:
                    columns.append(tokens[0])
                i += 1

            if not columns:
                raise ValueError(
                    f"Found loop_ without column labels in {filename} near line {i + 1}."
                )

            rows: List[List[str]] = []
            pending: List[str] = []
            ncols = len(columns)

            while i < len(lines):
                value_line = lines[i].strip()
                value_lower = value_line.lower()

                if (
                    value_lower.startswith("data_")
                    or value_lower == "loop_"
                    or value_lower.startswith("save_")
                ):
                    break

                if value_line.startswith("_") and not pending:
                    break

                if not value_line or value_line.startswith("#"):
                    i += 1
                    continue

                pending.extend(_split_star_line(lines[i]))
                while len(pending) >= ncols:
                    rows.append(pending[:ncols])
                    pending = pending[ncols:]
                i += 1

            if pending:
                raise ValueError(
                    f"Incomplete STAR loop row in {filename} near line {i + 1}: "
                    f"expected {ncols} values, found {len(pending)} trailing values."
                )

            star["blocks"][current_block]["loops"].append(
                {"columns": columns, "rows": rows}
            )
            continue

        if stripped.startswith("_"):
            tokens = _split_star_line(raw)
            if len(tokens) < 2:
                raise ValueError(
                    f"STAR scalar item has no value in {filename}, line {i + 1}: {stripped}"
                )
            star["blocks"][current_block]["items"].append((tokens[0], tokens[1]))
            i += 1
            continue

        raise ValueError(
            f"Unsupported or malformed STAR syntax in {filename}, line {i + 1}: {stripped}"
        )

    return star


def _quote_star_value(value: object) -> str:
    """Convert a value to a safe single-line STAR token."""
    text = str(value)
    if text == "":
        return "''"
    if (
        any(char.isspace() for char in text)
        or text.startswith("#")
        or "'" in text
        or '"' in text
    ):
        if "'" not in text:
            return f"'{text}'"
        if '"' not in text:
            return f'"{text}"'
        raise ValueError(
            f"Cannot safely write STAR token containing both quote types: {text!r}"
        )
    return text


def write_star_text(star: dict, filename: str) -> None:
    """Write the in-memory STAR representation in canonical RELION-style text."""
    output_parent = Path(filename).expanduser().parent
    output_parent.mkdir(parents=True, exist_ok=True)

    with open(filename, "w", encoding="utf-8", newline="\n") as handle:
        preamble = star.get("preamble", [])
        for line in preamble:
            handle.write(line.rstrip("\n") + "\n")
        if preamble and preamble[-1].strip():
            handle.write("\n")

        for block_name, block in star["blocks"].items():
            handle.write(f"data_{block_name}\n\n")

            for label, value in block.get("items", []):
                handle.write(f"{label} {_quote_star_value(value)}\n")
            if block.get("items"):
                handle.write("\n")

            for loop in block.get("loops", []):
                columns = loop["columns"]
                rows = loop["rows"]
                handle.write("loop_\n")
                for column_number, label in enumerate(columns, start=1):
                    handle.write(f"{label} #{column_number}\n")
                for row in rows:
                    if len(row) != len(columns):
                        raise ValueError(
                            f"Cannot write data_{block_name}: row has {len(row)} values "
                            f"but table has {len(columns)} columns."
                        )
                    handle.write(
                        " ".join(_quote_star_value(value) for value in row) + "\n"
                    )
                handle.write("\n")


def find_star_loop(star: dict, block_name: str, required_labels: Sequence[str]) -> dict:
    """Find a loop by block name, with a label-based fallback."""
    required = {_normalise_label(label) for label in required_labels}
    requested_block = _normalise_block_name(block_name)

    for name, block in star["blocks"].items():
        if _normalise_block_name(name) != requested_block:
            continue
        for loop in block["loops"]:
            labels = {_normalise_label(label) for label in loop["columns"]}
            if required.issubset(labels):
                return loop

    for block in star["blocks"].values():
        for loop in block["loops"]:
            labels = {_normalise_label(label) for label in loop["columns"]}
            if required.issubset(labels):
                return loop

    missing = ", ".join(sorted(required_labels))
    raise KeyError(
        f"Could not find a STAR loop containing [{missing}] in data_{block_name}."
    )


def star_column_index(loop: dict, label: str) -> int:
    wanted = _normalise_label(label)
    for index, existing_label in enumerate(loop["columns"]):
        if _normalise_label(existing_label) == wanted:
            return index
    raise KeyError(f"STAR column {label} was not found.")


def has_star_column(loop: dict, label: str) -> bool:
    wanted = _normalise_label(label)
    return any(_normalise_label(column) == wanted for column in loop["columns"])


def ensure_star_column(loop: dict, label: str, default_value: str = "?") -> int:
    """Return a column index, adding the column when it is missing."""
    if has_star_column(loop, label):
        return star_column_index(loop, label)
    output_label = label if label.startswith("_") else f"_{label}"
    loop["columns"].append(output_label)
    for row in loop["rows"]:
        row.append(str(default_value))
    return len(loop["columns"]) - 1


def get_star_column(loop: dict, label: str) -> List[str]:
    index = star_column_index(loop, label)
    return [row[index] for row in loop["rows"]]


def set_star_column(loop: dict, label: str, values: Iterable[object]) -> None:
    index = star_column_index(loop, label)
    values = list(values)
    if len(values) != len(loop["rows"]):
        raise ValueError(
            f"Cannot set {label}: got {len(values)} values for {len(loop['rows'])} rows."
        )
    for row, value in zip(loop["rows"], values):
        row[index] = str(value)


def set_star_column_constant(loop: dict, label: str, value: object) -> None:
    index = star_column_index(loop, label)
    text = str(value)
    for row in loop["rows"]:
        row[index] = text


# -----------------------------------------------------------------------------
# Optics-group metadata
# -----------------------------------------------------------------------------


def canonical_optics_group(value: object) -> str:
    """Normalise numeric optics-group values such as '1' and '1.0' to '1'."""
    text = str(value).strip()
    try:
        number = float(text)
    except ValueError:
        return text
    if np.isfinite(number) and number.is_integer():
        return str(int(number))
    return text


@dataclass(frozen=True)
class OpticsInfo:
    group_key: str
    row_index: int
    pixel_size: float
    image_size: Optional[int]


def build_optics_mapping(optics_loop: dict) -> Tuple[OrderedDict, bool]:
    if not optics_loop["rows"]:
        raise ValueError("The data_optics table contains no rows.")

    pixel_index = star_column_index(optics_loop, "rlnImagePixelSize")
    has_group_column = has_star_column(optics_loop, "rlnOpticsGroup")
    group_index = (
        star_column_index(optics_loop, "rlnOpticsGroup")
        if has_group_column
        else None
    )
    size_index = (
        star_column_index(optics_loop, "rlnImageSize")
        if has_star_column(optics_loop, "rlnImageSize")
        else None
    )

    if not has_group_column and len(optics_loop["rows"]) != 1:
        raise ValueError(
            "data_optics has multiple rows but no rlnOpticsGroup column, so particles "
            "cannot be mapped unambiguously to optics groups."
        )

    mapping: OrderedDict[str, OpticsInfo] = OrderedDict()
    for row_index, row in enumerate(optics_loop["rows"]):
        group_key = (
            canonical_optics_group(row[group_index])
            if group_index is not None
            else "1"
        )
        if group_key in mapping:
            raise ValueError(f"Duplicate optics group {group_key!r} in data_optics.")

        try:
            pixel_size = float(row[pixel_index])
        except ValueError as exc:
            raise ValueError(
                f"Invalid rlnImagePixelSize in optics group {group_key}: "
                f"{row[pixel_index]!r}"
            ) from exc
        if not np.isfinite(pixel_size) or pixel_size <= 0:
            raise ValueError(
                f"rlnImagePixelSize must be positive in optics group {group_key}; "
                f"got {pixel_size}."
            )

        image_size: Optional[int] = None
        if size_index is not None:
            size_text = str(row[size_index]).strip()
            if size_text not in {"", "?", "."}:
                try:
                    image_size = int(round(float(size_text)))
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid rlnImageSize in optics group {group_key}: {size_text!r}"
                    ) from exc
                if image_size <= 0:
                    raise ValueError(
                        f"rlnImageSize must be positive in optics group {group_key}; "
                        f"got {image_size}."
                    )

        mapping[group_key] = OpticsInfo(
            group_key=group_key,
            row_index=row_index,
            pixel_size=pixel_size,
            image_size=image_size,
        )

    return mapping, has_group_column


def map_particles_to_optics(
    particles_loop: dict, optics_mapping: OrderedDict
) -> List[str]:
    if has_star_column(particles_loop, "rlnOpticsGroup"):
        particle_groups = [
            canonical_optics_group(value)
            for value in get_star_column(particles_loop, "rlnOpticsGroup")
        ]
    else:
        if len(optics_mapping) != 1:
            raise ValueError(
                "data_particles has no rlnOpticsGroup column while data_optics contains "
                "multiple groups."
            )
        only_group = next(iter(optics_mapping))
        particle_groups = [only_group] * len(particles_loop["rows"])

    unknown = sorted({group for group in particle_groups if group not in optics_mapping})
    if unknown:
        raise ValueError(
            "Particles reference optics groups missing from data_optics: "
            + ", ".join(unknown)
        )
    return particle_groups


def update_optics_image_sizes(
    optics_loop: dict,
    optics_mapping: OrderedDict,
    output_sizes: Dict[str, int],
    forced_crop_size: int,
) -> None:
    size_index = ensure_star_column(optics_loop, "rlnImageSize", "?")
    for group_key, info in optics_mapping.items():
        size: Optional[int]
        if forced_crop_size > 0:
            size = forced_crop_size
        else:
            size = output_sizes.get(group_key, info.image_size)
        if size is not None:
            optics_loop["rows"][info.row_index][size_index] = str(int(size))


# -----------------------------------------------------------------------------
# Source image references and cached MRC access
# -----------------------------------------------------------------------------


@dataclass
class ParticleRecord:
    row_index: int
    image_name: str
    group_key: str
    source_index: Optional[int] = None
    source_path: Optional[str] = None
    parse_error: Optional[str] = None


def resolve_source_path(path_text: str, star_directory: Path) -> str:
    """
    Resolve a source path robustly.

    Existing paths relative to the current working directory retain their normal
    RELION behaviour. If not found there, a path relative to the STAR directory is
    tried as a fallback.
    """
    expanded = Path(os.path.expandvars(os.path.expanduser(path_text)))
    if expanded.is_absolute():
        return str(expanded)

    cwd_candidate = expanded
    if cwd_candidate.exists():
        return str(cwd_candidate.resolve())

    star_candidate = star_directory / expanded
    if star_candidate.exists():
        return str(star_candidate.resolve())

    # Preserve the conventional current-working-directory interpretation in the
    # eventual error message if neither candidate exists.
    return str(cwd_candidate.absolute())


def parse_image_reference(image_name: str, star_directory: Path) -> Tuple[int, str]:
    """Parse RELION's 1-based 'index@path' reference; plain paths mean index 0."""
    text = str(image_name).strip()
    if not text:
        raise ValueError("empty rlnImageName")

    if "@" in text:
        index_text, path_text = text.split("@", 1)
        try:
            relion_index = int(index_text)
        except ValueError as exc:
            raise ValueError(f"invalid image index {index_text!r}") from exc
        if relion_index < 0:
            raise ValueError(f"image index must be >= 0, got {relion_index}")
        # RELION is 1-based. A zero index is accepted as a compatibility fallback.
        source_index = 0 if relion_index == 0 else relion_index - 1
    else:
        source_index = 0
        path_text = text

    if not path_text:
        raise ValueError("empty MRC path after '@'")
    return source_index, resolve_source_path(path_text, star_directory)


def _detach_float32(image: np.ndarray) -> np.ndarray:
    """Detach an image from a memmap and convert it to native float32."""
    return np.array(image, dtype=np.float32, copy=True, order="C")


class MRCImageCache:
    """LRU cache of open mrcfile memory maps."""

    def __init__(self, max_open: int = 16):
        if max_open <= 0:
            raise ValueError("max_open must be positive")
        self.max_open = max_open
        self._files: OrderedDict[str, object] = OrderedDict()
        self.opens = 0
        self.cache_hits = 0
        self.evictions = 0
        self.reads = 0
        self.layout_counts: Counter[str] = Counter()

    def _get(self, filename: str):
        if filename in self._files:
            mrc = self._files.pop(filename)
            self._files[filename] = mrc
            self.cache_hits += 1
            return mrc

        if not os.path.isfile(filename):
            raise FileNotFoundError(f"MRC file does not exist: {filename}")

        mrc = mrcfile.mmap(filename, mode="r", permissive=True)
        if mrc.data is None:
            mrc.close()
            raise ValueError(f"MRC file contains no readable data: {filename}")

        self._files[filename] = mrc
        self.opens += 1

        while len(self._files) > self.max_open:
            _, old_mrc = self._files.popitem(last=False)
            old_mrc.close()
            self.evictions += 1
        return mrc

    @staticmethod
    def _header_dims(mrc) -> Tuple[int, int, int]:
        return (
            int(mrc.header.nz),
            int(mrc.header.ny),
            int(mrc.header.nx),
        )

    @staticmethod
    def _extract_image(
        data: np.ndarray,
        index: int,
        expected_box: Optional[int],
        filename: str,
        header_dims: Tuple[int, int, int],
    ) -> Tuple[np.ndarray, str]:
        if index < 0:
            raise IndexError(f"negative image index {index}")

        arr = np.asarray(data)
        # Accommodate unusual wrappers such as (N, 1, box, box).
        while arr.ndim > 3 and 1 in arr.shape:
            singleton_axis = next(axis for axis, size in enumerate(arr.shape) if size == 1)
            arr = np.squeeze(arr, axis=singleton_axis)

        box = int(expected_box) if expected_box is not None else None
        if box is not None and box <= 0:
            box = None

        if arr.ndim == 3:
            shape = arr.shape

            # Standard MRC stack orientation.
            if index < shape[0] and (
                box is None or shape[1:] == (box, box)
            ):
                return arr[index, :, :], "3d-standard"

            # Nonstandard axis order fallbacks.
            if box is not None and shape[:2] == (box, box) and index < shape[2]:
                return arr[:, :, index], "3d-stack-last-axis"
            if (
                box is not None
                and (shape[0], shape[2]) == (box, box)
                and index < shape[1]
            ):
                return arr[:, index, :], "3d-stack-middle-axis"

            # Header/image-size metadata can occasionally be stale. Prefer the
            # standard MRC axis if it still yields a 2-D image.
            if index < shape[0]:
                return arr[index, :, :], "3d-standard-box-mismatch"

            raise IndexError(
                f"Image index {index} is outside MRC data shape {shape}; "
                f"header(nz,ny,nx)={header_dims}, file={filename}"
            )

        if arr.ndim == 2:
            rows, cols = arr.shape

            # Crucial cryoSPARC/single-image case: mrcfile exposes nz=1 as a 2-D
            # array. Do not index it again, otherwise one row (shape=(box,)) is read.
            if box is not None and (rows, cols) == (box, box):
                if index != 0:
                    raise IndexError(
                        f"{filename} is a single 2-D image of shape {arr.shape}, but "
                        f"rlnImageName requests zero-based index {index}."
                    )
                return arr, "2d-single-image"

            if box is not None:
                pixels_per_image = box * box

                # One flattened image per row: (N, box*box).
                if cols == pixels_per_image:
                    if index >= rows:
                        raise IndexError(
                            f"Image index {index} outside flattened-row stack with {rows} images"
                        )
                    return arr[index, :].reshape(box, box), "2d-flattened-rows"

                # Vertically concatenated images: (N*box, box).
                if cols == box and rows % box == 0:
                    n_images = rows // box
                    if index >= n_images:
                        raise IndexError(
                            f"Image index {index} outside vertical stack with {n_images} images"
                        )
                    start = index * box
                    return arr[start : start + box, :], "2d-vertical-concatenation"

                # Horizontally concatenated images: (box, N*box).
                if rows == box and cols % box == 0:
                    n_images = cols // box
                    if index >= n_images:
                        raise IndexError(
                            f"Image index {index} outside horizontal stack with {n_images} images"
                        )
                    start = index * box
                    return arr[:, start : start + box], "2d-horizontal-concatenation"

                # One flattened image per column: (box*box, N).
                if rows == pixels_per_image:
                    if index >= cols:
                        raise IndexError(
                            f"Image index {index} outside flattened-column stack with {cols} images"
                        )
                    return arr[:, index].reshape(box, box), "2d-flattened-columns"

                # Last-resort contiguous flattening for an unusual 2-D storage layout.
                if arr.size % pixels_per_image == 0:
                    n_images = arr.size // pixels_per_image
                    if index >= n_images:
                        raise IndexError(
                            f"Image index {index} outside generic flattened stack with "
                            f"{n_images} images"
                        )
                    start = index * pixels_per_image
                    return (
                        arr.reshape(-1)[start : start + pixels_per_image].reshape(box, box),
                        "2d-generic-flattened",
                    )

            # Without usable optics metadata, a square 2-D array can still be
            # identified safely as one image.
            if rows == cols and index == 0:
                return arr, "2d-single-image-inferred"

            raise ValueError(
                f"Cannot interpret 2-D MRC data shape {arr.shape} as particle images. "
                f"expected_box={expected_box}, header(nz,ny,nx)={header_dims}, "
                f"file={filename}"
            )

        if arr.ndim == 1:
            if box is None:
                raise ValueError(
                    f"Cannot reshape 1-D MRC data without rlnImageSize: shape={arr.shape}, "
                    f"header(nz,ny,nx)={header_dims}, file={filename}"
                )
            pixels_per_image = box * box
            if arr.size % pixels_per_image != 0:
                raise ValueError(
                    f"1-D MRC data length {arr.size} is not divisible by box^2 "
                    f"({box}^2={pixels_per_image}); file={filename}"
                )
            n_images = arr.size // pixels_per_image
            if index >= n_images:
                raise IndexError(
                    f"Image index {index} outside 1-D flattened stack with {n_images} images"
                )
            start = index * pixels_per_image
            return (
                arr[start : start + pixels_per_image].reshape(box, box),
                "1d-flattened",
            )

        raise ValueError(
            f"Unsupported MRC data dimensionality ndim={arr.ndim}, shape={arr.shape}; "
            f"header(nz,ny,nx)={header_dims}, file={filename}"
        )

    def read_image(
        self, filename: str, index: int, expected_box: Optional[int]
    ) -> np.ndarray:
        mrc = self._get(filename)
        image, layout = self._extract_image(
            mrc.data,
            index=index,
            expected_box=expected_box,
            filename=filename,
            header_dims=self._header_dims(mrc),
        )
        if image.ndim != 2:
            raise ValueError(
                f"Internal error: interpreted image is not 2-D; shape={image.shape}, "
                f"file={filename}"
            )
        self.reads += 1
        self.layout_counts[layout] += 1
        return _detach_float32(image)

    def close(self) -> None:
        while self._files:
            _, mrc = self._files.popitem(last=False)
            mrc.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    def print_summary(self) -> None:
        print(
            f"MRC I/O: {self.reads} image reads, {self.opens} file opens, "
            f"{self.cache_hits} cache hits, {self.evictions} evictions."
        )
        if self.layout_counts:
            layouts = ", ".join(
                f"{name}={count}" for name, count in sorted(self.layout_counts.items())
            )
            print(f"Detected MRC layouts: {layouts}")


# -----------------------------------------------------------------------------
# Image processing
# -----------------------------------------------------------------------------


def translation_twoD_image(image: np.ndarray, trans_X: float, trans_Y: float) -> np.ndarray:
    """Apply an integer translation with zero padding and no wraparound."""
    trans_x = int(np.floor(float(trans_X) + 0.5))
    trans_y = int(np.floor(float(trans_Y) + 0.5))
    ysize, xsize = image.shape

    translated = np.zeros_like(image)

    dst_x0 = max(0, trans_x)
    dst_y0 = max(0, trans_y)
    dst_x1 = min(xsize, xsize + trans_x)
    dst_y1 = min(ysize, ysize + trans_y)

    if dst_x0 >= dst_x1 or dst_y0 >= dst_y1:
        return translated

    src_x0 = max(0, -trans_x)
    src_y0 = max(0, -trans_y)
    src_x1 = src_x0 + (dst_x1 - dst_x0)
    src_y1 = src_y0 + (dst_y1 - dst_y0)

    translated[dst_y0:dst_y1, dst_x0:dst_x1] = image[
        src_y0:src_y1, src_x0:src_x1
    ]
    return translated


def center_crop(image: np.ndarray, crop_size: int) -> np.ndarray:
    if crop_size <= 0:
        return np.array(image, copy=True, order="C")
    ysize, xsize = image.shape
    if crop_size > ysize or crop_size > xsize:
        raise ValueError(
            f"Requested crop size {crop_size} exceeds image shape {image.shape}."
        )
    start_y = (ysize - crop_size) // 2
    start_x = (xsize - crop_size) // 2
    return np.array(
        image[start_y : start_y + crop_size, start_x : start_x + crop_size],
        copy=True,
        order="C",
    )


@lru_cache(maxsize=16)
def _radial_integer_bins(shape: Tuple[int, int]) -> np.ndarray:
    ysize, xsize = shape
    center_y, center_x = ysize // 2, xsize // 2
    yy, xx = np.meshgrid(np.arange(ysize), np.arange(xsize), indexing="ij")
    return np.rint(
        np.sqrt((yy - center_y) ** 2 + (xx - center_x) ** 2)
    ).astype(np.int32)


def getSpectrum_divideBySpectrum(img: np.ndarray) -> np.ndarray:
    img_fft = np.fft.fftshift(np.fft.fft2(img))
    dist = _radial_integer_bins(tuple(img.shape))
    max_radius = int(dist.max()) + 1
    spectrum = np.bincount(
        dist.ravel(), weights=np.abs(img_fft).ravel(), minlength=max_radius
    )
    counts = np.bincount(dist.ravel(), minlength=max_radius)

    radial_mean = np.zeros_like(spectrum, dtype=np.float64)
    valid_counts = counts > 0
    radial_mean[valid_counts] = spectrum[valid_counts] / counts[valid_counts]

    inverse = np.ones_like(radial_mean, dtype=np.float64)
    valid_power = radial_mean > np.finfo(np.float64).eps
    inverse[valid_power] = 1.0 / radial_mean[valid_power]
    inverse[0] = 1.0

    img_fft *= inverse[dist]
    return np.fft.ifft2(np.fft.ifftshift(img_fft)).real


@lru_cache(maxsize=16)
def _coordinate_grids(shape: Tuple[int, int]) -> Tuple[np.ndarray, np.ndarray]:
    ysize, xsize = shape
    yy, xx = np.meshgrid(
        np.arange(ysize, dtype=np.float64),
        np.arange(xsize, dtype=np.float64),
        indexing="ij",
    )
    return yy, xx


def normalize(img: np.ndarray) -> np.ndarray:
    """Subtract the least-squares plane and scale to zero mean/unit standard deviation."""
    image = np.asarray(img, dtype=np.float64)
    yy, xx = _coordinate_grids(tuple(image.shape))

    # Solve the 3-parameter least-squares plane z = a*x + b*y + c without
    # building a large (N, 3) design matrix.
    sx = float(np.sum(xx))
    sy = float(np.sum(yy))
    sxx = float(np.sum(xx * xx))
    syy = float(np.sum(yy * yy))
    sxy = float(np.sum(xx * yy))
    n = float(image.size)

    matrix = np.array(
        [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n]], dtype=np.float64
    )
    rhs = np.array(
        [np.sum(xx * image), np.sum(yy * image), np.sum(image)], dtype=np.float64
    )
    try:
        plane_a, plane_b, plane_c = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError:
        plane_a = plane_b = 0.0
        plane_c = float(np.mean(image))

    image = image - (plane_a * xx + plane_b * yy + plane_c)
    mean = float(np.mean(image))
    std = float(np.std(image))
    if std > np.finfo(np.float64).eps:
        image = (image - mean) / std
    else:
        image = image - mean
    return image


@lru_cache(maxsize=64)
def _highpass_mask(
    shape: Tuple[int, int], high_pass: float, angpix: float
) -> np.ndarray:
    ysize, xsize = shape
    # Preserve the original script's box-size convention while supporting
    # rectangular images if encountered.
    ori_size = ysize
    filter_edge_width = 4.0
    ires_filter = round((ori_size * angpix) / high_pass)
    half_width = filter_edge_width / 2.0

    edge_low = max(0.0, (ires_filter - half_width) / ori_size)
    edge_high = max(edge_low, (ires_filter + half_width) / ori_size)
    edge_width = edge_high - edge_low

    y = np.arange(ysize, dtype=np.float64) - ysize // 2
    x = np.arange(xsize, dtype=np.float64) - xsize // 2
    yy, xx = np.meshgrid(y, x, indexing="ij")
    radius = np.sqrt(xx * xx + yy * yy) / ori_size

    mask = np.ones(shape, dtype=np.float64)
    mask[radius < edge_low] = 0.0
    transition = (radius >= edge_low) & (radius <= edge_high)
    if edge_width > 0:
        mask[transition] = 0.5 - 0.5 * np.cos(
            np.pi * (radius[transition] - edge_low) / edge_width
        )
    return mask


def highpassfilter(img: np.ndarray, high_pass: float, angpix: float) -> np.ndarray:
    if high_pass <= 0:
        raise ValueError(f"highpassFreq must be positive, got {high_pass}")
    if angpix <= 0:
        raise ValueError(f"Pixel size must be positive, got {angpix}")
    fft = np.fft.fftshift(np.fft.fft2(img))
    fft *= _highpass_mask(tuple(img.shape), float(high_pass), float(angpix))
    return np.fft.ifft2(np.fft.ifftshift(fft)).real


# -----------------------------------------------------------------------------
# Output stack writer
# -----------------------------------------------------------------------------


class OutputStackWriter:
    """
    Write batches while keeping each output stack within one optics group.

    Processing is group-by-group, so only one batch is resident in RAM. Global file
    numbering preserves the familiar output_0001.mrcs naming pattern.
    """

    def __init__(self, output_root: str, batch_size: int, dry_run: bool = False):
        if batch_size <= 0:
            raise ValueError("batch_size must be positive")
        self.output_root = output_root
        self.batch_size = batch_size
        self.dry_run = dry_run
        self.next_stack_number = 1
        self.current_group: Optional[str] = None
        self.current_pixel_size: Optional[float] = None
        self.current_filename: Optional[str] = None
        self.images: List[np.ndarray] = []
        self.placeholder_count = 0
        self.files_written = 0
        self.images_written = 0

    def _reserve_filename(self) -> str:
        filename = f"{self.output_root}_{self.next_stack_number:04d}.mrcs"
        self.next_stack_number += 1
        Path(filename).expanduser().parent.mkdir(parents=True, exist_ok=True)
        return filename

    def begin_group(self, group_key: str, pixel_size: float) -> None:
        if self.current_group is not None and self.current_group != group_key:
            self.flush()
        self.current_group = group_key
        self.current_pixel_size = float(pixel_size)

    def _current_count(self) -> int:
        return self.placeholder_count if self.dry_run else len(self.images)

    def append(self, image: Optional[np.ndarray] = None) -> str:
        if self.current_group is None or self.current_pixel_size is None:
            raise RuntimeError("begin_group() must be called before append().")
        if self.current_filename is None:
            self.current_filename = self._reserve_filename()

        if self.dry_run:
            self.placeholder_count += 1
        else:
            if image is None or image.ndim != 2:
                raise ValueError("Output image must be a 2-D NumPy array.")
            image32 = np.asarray(image, dtype=np.float32)
            if self.images and image32.shape != self.images[0].shape:
                raise ValueError(
                    f"Output image shape changed within optics group {self.current_group}: "
                    f"{self.images[0].shape} -> {image32.shape}"
                )
            self.images.append(np.array(image32, copy=True, order="C"))

        position = self._current_count()
        reference = f"{position}@{self.current_filename}"
        if position >= self.batch_size:
            self.flush()
        return reference

    def flush(self) -> None:
        count = self._current_count()
        if count == 0:
            self.current_filename = None
            return

        filename = self.current_filename
        if filename is None:
            raise RuntimeError("Internal error: non-empty output batch has no filename.")

        if not self.dry_run:
            stack = np.stack(self.images, axis=0).astype(np.float32, copy=False)
            with mrcfile.new(filename, overwrite=True) as output_mrc:
                output_mrc.set_data(stack)
                output_mrc.set_image_stack()
                output_mrc.voxel_size = self.current_pixel_size
                output_mrc.update_header_stats()
            self.files_written += 1
            self.images_written += count
            print(
                f"Saved {filename}: {count} images, shape={stack.shape[1:]}, "
                f"optics_group={self.current_group}, "
                f"pixel_size={self.current_pixel_size:g} Å/pixel"
            )

        self.images = []
        self.placeholder_count = 0
        self.current_filename = None

    def finish(self) -> None:
        self.flush()
        self.current_group = None
        self.current_pixel_size = None


# -----------------------------------------------------------------------------
# Command line and main workflow
# -----------------------------------------------------------------------------


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Shift and crop particle images referenced by a RELION STAR file, with "
            "multi-optics-group support and cached MRC I/O."
        )
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {PROGRAM_VERSION}")
    parser.add_argument("--star_name", type=str, required=True, help="Input STAR file")
    parser.add_argument(
        "--output_root_name",
        type=str,
        default="output",
        help=(
            "Output root. Stacks are ROOT_0001.mrcs, ROOT_0002.mrcs, ... and "
            "the STAR file is ROOT.star. Default: output"
        ),
    )
    parser.add_argument(
        "--newboxsize",
        type=int,
        default=128,
        help="Output crop size in pixels. <=0 keeps the full input image. Default: 128",
    )
    parser.add_argument(
        "--batchsize",
        type=int,
        default=5000,
        help="Maximum images per output stack. Default: 5000",
    )
    parser.add_argument(
        "--max_open_mrc",
        type=int,
        default=16,
        help="Maximum source MRC memory maps kept in the LRU cache. Default: 16",
    )
    parser.add_argument(
        "--progress_every",
        type=int,
        default=1000,
        help="Print progress every N attempted particles; 0 disables. Default: 1000",
    )
    parser.add_argument(
        "--dohighpass",
        action="store_true",
        help="Apply a high-pass filter to output particles",
    )
    parser.add_argument(
        "--highpassFreq",
        type=float,
        default=50.0,
        help="High-pass frequency in Angstrom. Default: 50",
    )
    parser.add_argument(
        "--dowhitening",
        action="store_true",
        help="Apply radial-spectrum whitening followed by normalization",
    )
    parser.add_argument(
        "--doSkipShifting",
        action="store_true",
        help="Do not apply rlnOriginXAngst/rlnOriginYAngst shifts",
    )
    parser.add_argument(
        "--doOnlyMakeStar",
        action="store_true",
        help="Generate only the rewritten STAR file; do not read/write MRC data",
    )
    parser.add_argument(
        "--fail_fast",
        action="store_true",
        help="Stop on the first bad particle instead of skipping it",
    )
    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    if args.batchsize <= 0:
        raise ValueError("--batchsize must be positive")
    if args.max_open_mrc <= 0:
        raise ValueError("--max_open_mrc must be positive")
    if args.progress_every < 0:
        raise ValueError("--progress_every must be >= 0")
    if args.dohighpass and args.highpassFreq <= 0:
        raise ValueError("--highpassFreq must be positive when --dohighpass is used")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    validate_arguments(args)

    star_name = str(Path(args.star_name).expanduser())
    output_root = str(Path(args.output_root_name).expanduser())
    crop_size = int(args.newboxsize)

    print(f"Reading STAR file: {star_name}")
    star_data = read_star_text(star_name)
    particles_loop = find_star_loop(
        star_data,
        "particles",
        ["rlnOriginXAngst", "rlnOriginYAngst", "rlnImageName"],
    )
    optics_loop = find_star_loop(
        star_data,
        "optics",
        ["rlnImagePixelSize"],
    )

    optics_mapping, _ = build_optics_mapping(optics_loop)
    particle_groups = map_particles_to_optics(particles_loop, optics_mapping)

    n_particles = len(particles_loop["rows"])
    origin_x_angst = np.asarray(
        get_star_column(particles_loop, "rlnOriginXAngst"), dtype=np.float64
    )
    origin_y_angst = np.asarray(
        get_star_column(particles_loop, "rlnOriginYAngst"), dtype=np.float64
    )
    image_names = get_star_column(particles_loop, "rlnImageName")

    if not (
        len(origin_x_angst)
        == len(origin_y_angst)
        == len(image_names)
        == len(particle_groups)
        == n_particles
    ):
        raise RuntimeError("STAR particle columns have inconsistent lengths.")

    print(
        f"Particles: {n_particles}; optics groups: {len(optics_mapping)} "
        f"({', '.join(optics_mapping.keys())})"
    )

    star_directory = Path(star_name).absolute().parent
    records_by_group: Dict[str, List[ParticleRecord]] = defaultdict(list)

    for row_index, (image_name, group_key) in enumerate(
        zip(image_names, particle_groups)
    ):
        record = ParticleRecord(
            row_index=row_index,
            image_name=str(image_name),
            group_key=group_key,
        )
        if not args.doOnlyMakeStar:
            try:
                source_index, source_path = parse_image_reference(
                    record.image_name, star_directory
                )
                record.source_index = source_index
                record.source_path = source_path
            except Exception as exc:  # stored and handled in the normal skip path
                record.parse_error = str(exc)
        records_by_group[group_key].append(record)

    # Sorting within each optics group improves locality dramatically for STAR files
    # whose rows interleave many source stacks. STAR row order itself is preserved by
    # assigning references back to their original row indices.
    if not args.doOnlyMakeStar:
        for records in records_by_group.values():
            records.sort(
                key=lambda record: (
                    record.source_path is None,
                    record.source_path or "",
                    record.source_index if record.source_index is not None else -1,
                    record.row_index,
                )
            )

    new_image_names: List[Optional[str]] = [None] * n_particles
    output_sizes: Dict[str, int] = {}
    failed = 0
    attempted = 0

    writer = OutputStackWriter(
        output_root=output_root,
        batch_size=args.batchsize,
        dry_run=args.doOnlyMakeStar,
    )
    cache = MRCImageCache(max_open=args.max_open_mrc)

    try:
        for group_key, optics_info in optics_mapping.items():
            group_records = records_by_group.get(group_key, [])
            if not group_records:
                continue

            writer.begin_group(group_key, optics_info.pixel_size)
            print(
                f"Processing optics group {group_key}: {len(group_records)} particles, "
                f"pixel_size={optics_info.pixel_size:g} Å/pixel, "
                f"input_box={optics_info.image_size if optics_info.image_size is not None else 'unknown'}"
            )

            for record in group_records:
                attempted += 1
                try:
                    if args.doOnlyMakeStar:
                        if crop_size <= 0 and optics_info.image_size is None:
                            raise ValueError(
                                "--doOnlyMakeStar with --newboxsize <= 0 requires "
                                "rlnImageSize in data_optics."
                            )
                        output_size = (
                            crop_size if crop_size > 0 else int(optics_info.image_size)
                        )
                        existing_size = output_sizes.get(group_key)
                        if existing_size is not None and existing_size != output_size:
                            raise ValueError(
                                f"Inconsistent output size in optics group {group_key}: "
                                f"{existing_size} vs {output_size}"
                            )
                        output_sizes[group_key] = output_size
                        new_image_names[record.row_index] = writer.append(None)
                        continue

                    if record.parse_error is not None:
                        raise ValueError(record.parse_error)
                    if record.source_path is None or record.source_index is None:
                        raise RuntimeError("Image reference was not parsed.")

                    image = cache.read_image(
                        record.source_path,
                        record.source_index,
                        optics_info.image_size,
                    )
                    if image.size == 0 or image.ndim != 2:
                        raise ValueError(
                            f"Bad source image shape {getattr(image, 'shape', None)}"
                        )

                    processed = image
                    if not args.doSkipShifting:
                        shift_x_pixels = (
                            origin_x_angst[record.row_index] / optics_info.pixel_size
                        )
                        shift_y_pixels = (
                            origin_y_angst[record.row_index] / optics_info.pixel_size
                        )
                        processed = translation_twoD_image(
                            processed, shift_x_pixels, shift_y_pixels
                        )

                    processed = center_crop(processed, crop_size)

                    if args.dowhitening:
                        processed = normalize(
                            getSpectrum_divideBySpectrum(processed)
                        )
                    if args.dohighpass:
                        processed = highpassfilter(
                            processed,
                            args.highpassFreq,
                            optics_info.pixel_size,
                        )

                    processed = np.asarray(processed, dtype=np.float32)
                    if processed.ndim != 2 or processed.size == 0:
                        raise ValueError(
                            f"Processed image is invalid: shape={processed.shape}"
                        )
                    if processed.shape[0] != processed.shape[1]:
                        raise ValueError(
                            f"RELION rlnImageSize expects square particles, but output "
                            f"shape is {processed.shape}."
                        )

                    output_size = int(processed.shape[0])
                    existing_size = output_sizes.get(group_key)
                    if existing_size is not None and existing_size != output_size:
                        raise ValueError(
                            f"Output size changed within optics group {group_key}: "
                            f"{existing_size} -> {output_size}"
                        )
                    output_sizes[group_key] = output_size

                    new_image_names[record.row_index] = writer.append(processed)

                except Exception as exc:
                    failed += 1
                    message = (
                        f"Failed particle row {record.row_index + 1}: "
                        f"{record.image_name}; error: {exc}"
                    )
                    if args.fail_fast:
                        raise RuntimeError(message) from exc
                    print("WARNING:", message, file=sys.stderr)

                if args.progress_every and attempted % args.progress_every == 0:
                    success_so_far = attempted - failed
                    print(
                        f"Progress: attempted {attempted}/{n_particles}, "
                        f"successful {success_so_far}, failed {failed}"
                    )

            writer.flush()

    finally:
        writer.finish()
        cache.close()

    valid_indices = [
        index for index, image_name in enumerate(new_image_names) if image_name is not None
    ]
    valid_names = [str(new_image_names[index]) for index in valid_indices]

    if not args.doOnlyMakeStar:
        particles_loop["rows"] = [
            particles_loop["rows"][index] for index in valid_indices
        ]
    elif len(valid_indices) != n_particles:
        # This can only happen if a STAR-only metadata error was skipped.
        particles_loop["rows"] = [
            particles_loop["rows"][index] for index in valid_indices
        ]

    set_star_column(particles_loop, "rlnImageName", valid_names)

    if not args.doSkipShifting:
        set_star_column_constant(particles_loop, "rlnOriginXAngst", "0.0")
        set_star_column_constant(particles_loop, "rlnOriginYAngst", "0.0")

    update_optics_image_sizes(
        optics_loop,
        optics_mapping,
        output_sizes=output_sizes,
        forced_crop_size=crop_size,
    )

    output_star_name = f"{output_root}.star"
    write_star_text(star_data, output_star_name)

    successful = len(valid_indices)
    print(f"New STAR file saved as {output_star_name}")
    print(
        f"Completed: {successful} successful particles, {failed} failed particles, "
        f"{writer.files_written} output stack(s)."
    )
    if not args.doOnlyMakeStar:
        cache.print_summary()

    if successful == 0 and n_particles > 0:
        print("ERROR: no particles were processed successfully.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
