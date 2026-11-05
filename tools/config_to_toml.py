#!/usr/bin/env python3
"""Convert legacy BeAR .config files to the new TOML format.

Usage:
    python config_to_toml.py <retrieval_dir> [<retrieval_dir> ...]

For each directory this converts (if present):
    retrieval.config     -> retrieval.toml
    forward_model*.config -> forward_model*.toml   (all variants)
    post_process*.config  -> post_process*.toml    (all variants)

priors.config is intentionally left untouched (it already uses a name-keyed
row format). The original .config files are kept; delete them once you have
verified the .toml files.

The legacy files are positional; sections are identified here by their comment
header text, which is robust to the per-forward-model ordering differences.
"""

import sys
import os
import glob


def read_blocks(path):
    """Parse a legacy config file into a list of (header, [data_lines]).

    A data block is a maximal run of consecutive non-comment, non-blank lines.
    `header` is the text of the last comment line before it. A line counts as a
    comment if it starts with '#' or ends with '#' (the latter catches the
    decorative 'General retrieval config#' banner lines)."""
    blocks = []
    last_comment = ""
    current = []

    with open(path) as f:
        for raw in f:
            stripped = raw.strip()

            is_comment = stripped.startswith("#") or stripped.endswith("#")

            if stripped == "" or is_comment:
                if current:
                    blocks.append((last_comment, current))
                    current = []
                if is_comment:
                    last_comment = stripped.strip("#").strip()
                continue

            current.append(stripped)

    if current:
        blocks.append((last_comment, current))

    return blocks


def to_bool(token):
    return token.lower() in ("y", "yes", "1", "true")


def quote(s):
    return '"' + s.replace('\\', '\\\\').replace('"', '\\"') + '"'


def str_array(tokens):
    return "[" + ", ".join(quote(t) for t in tokens) + "]"


# ---------------------------------------------------------------------------
# retrieval.config -> retrieval.toml
# ---------------------------------------------------------------------------

def convert_retrieval(path):
    #ordered data blocks: use_gpu, omp, fwd_model, discretisation+res, opacity, error_infl, [highres]
    values = [lines[0] for _, lines in read_blocks(path)]

    use_gpu = to_bool(values[0].split()[0])
    nb_omp = int(values[1].split()[0])
    fwd_model = values[2].split()[0]
    disc_tokens = values[3].split()
    discretisation = disc_tokens[0]
    resolution = disc_tokens[1]
    opacity_folder = values[4].split()[0]
    use_error_inflation = to_bool(values[5].split()[0])

    highres = None
    for v in values[6:]:
        t = v.split()
        if t and t[0] == "spectral_resolution_highres":
            highres = t[1]

    out = []
    out.append("[general]")
    out.append("use_gpu = " + ("true" if use_gpu else "false"))
    out.append("nb_omp_threads = " + str(nb_omp) + "     # 0 = auto")
    out.append("")
    out.append("[retrieval]")
    out.append("forward_model_type = " + quote(fwd_model))
    out.append("spectral_discretisation = " + quote(discretisation))
    out.append("spectral_resolution = " + resolution)
    out.append("opacity_data_folder = " + quote(opacity_folder))
    out.append("use_error_inflation = " + ("true" if use_error_inflation else "false"))
    if highres is not None:
        out.append("spectral_resolution_highres = " + highres)
    out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# forward_model.config -> forward_model.toml
# ---------------------------------------------------------------------------

def hdr_has(header, *keys):
    h = header.lower()
    return any(k in h for k in keys)


def convert_forward_model(path):
    blocks = read_blocks(path)

    atmosphere = {}
    temperature = None            # (model, params)
    stellar = None                # (model, params)
    stellar_smooth = None         # value
    radiative_transfer = None     # (model, params)
    clouds = []                   # list of (model, params)
    modules = None                # list of (model, params); None = no module section
    chemistry = []                # list of (model, params)
    opacity_lists = []            # list of [(species, folder), ...]
    fit_mode = None               # transmission: mmw|sh|no
    variable_gravity = None       # transmission: bool

    def model_params(lines):
        toks = lines[0].split()
        return toks[0], toks[1:]

    def model_list(lines):
        out = []
        for ln in lines:
            toks = ln.split()
            if not toks:
                continue
            if toks[0].lower() in ("none", "no", "n"):
                continue
            out.append((toks[0], toks[1:]))
        return out

    for header, lines in blocks:
        if hdr_has(header, "number of levels", "grid points"):
            atmosphere["nb_grid_points"] = lines[0].split()[0]
        elif hdr_has(header, "bottom of atmosphere", "bottom pressure"):
            atmosphere["bottom_pressure"] = lines[0].split()[0]
        elif hdr_has(header, "top of atmosphere", "top pressure"):
            atmosphere["top_pressure"] = lines[0].split()[0]
        elif hdr_has(header, "mean molecular weight", "scale height"):
            fit_mode = lines[0].split()[0]
        elif hdr_has(header, "variable gravity"):
            variable_gravity = to_bool(lines[0].split()[0])
        elif hdr_has(header, "temperature"):
            temperature = model_params(lines)
        elif hdr_has(header, "stellar"):
            #first data line = model + params; an optional following line is the
            #high-res stellar smoothing sigma
            toks = lines[0].split()
            stellar = (toks[0], toks[1:])
            for ln in lines[1:]:
                t = ln.split()
                if t and t[0] == "highres_stellar_smooth_sigma":
                    stellar_smooth = t[1]
        elif hdr_has(header, "cloud"):
            clouds = model_list(lines)
        elif hdr_has(header, "radiative transfer"):
            radiative_transfer = model_params(lines)
        elif hdr_has(header, "module"):
            modules = model_list(lines)
        elif hdr_has(header, "chemical species", "chemistry"):
            chemistry = model_list(lines)
        elif hdr_has(header, "opacity"):
            entries = []
            for ln in lines:
                toks = ln.split()
                if len(toks) >= 2:
                    entries.append((toks[0], toks[1]))
            opacity_lists.append(entries)

    #Everything is emitted as top-level inline tables / arrays of inline tables,
    #in the original section order. Inline forms are the SAME TOML data model as
    #the [section]/[[array]] forms, so the C++ parser is unchanged; they are just
    #more compact and list-like.
    out = []

    def inline_model(model, params):
        return "{ model = " + quote(model) + ", params = " + str_array(params) + " }"

    def model_array(key, items):
        #an array of {model, params} inline tables, or [] for "none"
        if not items:
            out.append(key + " = []")
            return
        out.append(key + " = [")
        for m, p in items:
            out.append("  " + inline_model(m, p) + ",")
        out.append("]")

    #transmission-only scalars
    if fit_mode is not None:
        out.append("fit_mode = " + quote(fit_mode))
    if variable_gravity is not None:
        out.append("use_variable_gravity = " + ("true" if variable_gravity else "false"))
    if fit_mode is not None or variable_gravity is not None:
        out.append("")

    if atmosphere:
        body = ", ".join([
            "nb_grid_points = " + atmosphere.get("nb_grid_points", "0"),
            "bottom_pressure = " + atmosphere.get("bottom_pressure", "0"),
            "top_pressure = " + atmosphere.get("top_pressure", "0")])
        out.append("atmosphere = { " + body + " }")
        out.append("")

    if temperature is not None:
        out.append("temperature = " + inline_model(temperature[0], temperature[1]))
        out.append("")

    if stellar is not None:
        body = "model = " + quote(stellar[0]) + ", params = " + str_array(stellar[1])
        if stellar_smooth is not None:
            body += ", highres_stellar_smooth_sigma = " + stellar_smooth
        out.append("stellar_spectrum = { " + body + " }")
        out.append("")

    model_array("clouds", clouds)
    out.append("")

    if radiative_transfer is not None:
        out.append("radiative_transfer = " + inline_model(radiative_transfer[0], radiative_transfer[1]))
        out.append("")

    if modules is not None:
        model_array("modules", modules)
        out.append("")

    model_array("chemistry", chemistry)
    out.append("")

    opacity_keys = ["opacity", "opacity_highres"]
    for i, entries in enumerate(opacity_lists):
        key = opacity_keys[i] if i < len(opacity_keys) else "opacity_%d" % i
        if not entries:
            out.append(key + " = []")
        else:
            out.append(key + " = [")
            for species, folder in entries:
                out.append("  { species = " + quote(species) + ", folder = " + quote(folder) + " },")
            out.append("]")
        out.append("")

    return "\n".join(out)


# ---------------------------------------------------------------------------
# post_process.config -> post_process.toml
# ---------------------------------------------------------------------------

def convert_post_process(path):
    out = []
    for header, lines in read_blocks(path):
        val = lines[0].split()

        if hdr_has(header, "delete"):
            out.append("delete_sampler_files = " + ("true" if to_bool(val[0]) else "false"))
        elif hdr_has(header, "save spectra", "posterior spectra"):
            out.append("save_spectra = " + ("true" if to_bool(val[0]) else "false"))
        elif hdr_has(header, "effective temperature"):
            out.append("save_effective_temperatures = " + ("true" if to_bool(val[0]) else "false"))
        elif hdr_has(header, "temperature"):
            out.append("save_temperatures = " + ("true" if to_bool(val[0]) else "false"))
        elif hdr_has(header, "contribution"):
            out.append("save_contribution_functions = " + ("true" if to_bool(val[0]) else "false"))
        elif hdr_has(header, "chemical species", "chemistry"):
            species = [] if (len(val) == 1 and val[0].lower() == "none") else val
            out.append("species_to_save = " + str_array(species))

    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------

def convert_dir(folder):
    def do(src, dst, fn):
        text = fn(src)
        with open(dst, "w") as f:
            f.write(text if text.endswith("\n") else text + "\n")
        print("  wrote", os.path.relpath(dst))

    r = os.path.join(folder, "retrieval.config")
    if os.path.exists(r):
        do(r, os.path.join(folder, "retrieval.toml"), convert_retrieval)

    for src in sorted(glob.glob(os.path.join(folder, "forward_model*.config"))):
        do(src, src[:-len(".config")] + ".toml", convert_forward_model)

    for src in sorted(glob.glob(os.path.join(folder, "post_process*.config"))):
        do(src, src[:-len(".config")] + ".toml", convert_post_process)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    for folder in sys.argv[1:]:
        print("Converting", folder)
        convert_dir(folder)


if __name__ == "__main__":
    main()
