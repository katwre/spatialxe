#!/usr/bin/env python3

import json
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile


# Xenium full-resolution image: 1 pixel = 0.2125 µm (10x Genomics spec).
# https://kb.10xgenomics.com/s/article/11636252598925-What-are-the-Xenium-image-scale-factors
# Transcript x_location / y_location are in microns.
# To overlay transcripts on the full-res image: pixel = micron / pixel_size.
XENIUM_DEFAULT_PIXEL_SIZE_UM = 0.2125


# Mock run:
#nextflow run . -profile docker --input DATA/samplesheet.csv --method scs --mode image -c ./DATA/local.config --outdir results_scs_pipeline_test -ansi-log false
#nextflow run . -profile docker --input DATA/samplesheet.csv --method scs --mode image -c ./DATA/test_scs_local.config --outdir results_scs_pipeline_test -ansi-log false -resume
# nextflow run . -profile docker --input DATA/samplesheet.csv --method scs --mode image -c ./DATA/test_scs_local.config --outdir results_scs_pipeline_test -ansi-log false -resume
# nextflow run . -profile docker --input DATA/samplesheet_full.csv --method scs --mode image -c ./DATA/local.config --outdir results_scs_pipeline_test -ansi-log false

def _pick_column(df: pd.DataFrame,
                 candidates: list[str],
                 required: bool = True):
    for name in candidates:
        if name in df.columns:
            return name
    if required:
        raise ValueError(f"Could not find any of the required columns: {candidates}")
    return None


def _read_pixel_size(experiment_xenium_path: str) -> float:
    """Read pixel_size (µm/px) from experiment.xenium; fall back to 0.2125."""
    try:
        with open(experiment_xenium_path) as fh:
            meta = json.load(fh)
        return float(meta.get("pixel_size", XENIUM_DEFAULT_PIXEL_SIZE_UM))
    except Exception:
        return XENIUM_DEFAULT_PIXEL_SIZE_UM


def _build_density_grid(df: pd.DataFrame) -> np.ndarray:
    """Builds a 2D density grid from the scs_input DataFrame."""
    grid = np.zeros((df["row"].max() + 1, df["column"].max() + 1), dtype=np.float32)

    for row, column, counts in df[["row", "column", "counts"]].itertuples(index=False):
        grid[int(row), int(column)] += float(counts)

    return grid
def _plot_density_map(df, output_path) -> None:
    """Plots a density map of all transcripts."""
    import matplotlib.pyplot as plt
    import numpy as np
    
    density = _build_density_grid(df)

    fig, ax = plt.subplots(figsize=(8, 8))
    im = ax.imshow(np.log1p(density), cmap="magma", origin="lower")
    ax.set_title("All transcripts density map (log1p counts)")
    ax.set_xlabel("column")
    ax.set_ylabel("row")
    fig.colorbar(im, ax=ax, label="log1p(counts)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_hexbin_density(df, output_path) -> None:
    import pandas as pd
    import matplotlib.pyplot as plt

    x_col, y_col = "x_location", "y_location"
    pts = df[[x_col, y_col]].dropna()

    fig, ax = plt.subplots(figsize=(8,8))
    hb = ax.hexbin(pts[x_col], pts[y_col], gridsize=200, bins="log", 
                   mincnt=1, cmap="magma")
    ax.set_title("Transcript molecule density (hexbin, log)")
    ax.set_xlabel("x_location (µm)")
    ax.set_ylabel("y_location (µm)")
    ax.set_aspect("equal")
    fig.colorbar(hb, ax=ax, label="log10(count)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=250)
    plt.close(fig)


def convert_xenium_to_scs(parquet_path: str,
                          output_tsv: str,
                          output_bgi_tsv: str,
                          morphology_image_path: str,
                          output_morph2d_tif: str,
                          metrics_tsv: str,
                          experiment_xenium_path: str = "",
                          bin_size: float = 5.0):
    """
    Convert Xenium transcripts to SCS/BGI format with correct pixel-space coordinates.

    Xenium x_location / y_location are in microns.
    The morphology image full resolution is 0.2125 µm/px (Xenium spec).
    Coordinates are converted to pixels: pixel = micron / pixel_size.
    The morphology image is cropped to the pixel ROI covered by the transcripts.
    """

    pixel_size = _read_pixel_size(experiment_xenium_path) if experiment_xenium_path else XENIUM_DEFAULT_PIXEL_SIZE_UM
    print(f"[xenium2scs] pixel_size_um: {pixel_size}")

    transcripts = pd.read_parquet(parquet_path, engine="pyarrow")
    #print(transcripts.head())
    #     transcript_id     cell_id  ...  codeword_category is_gene
    #0  281474976711277  mpafcjmo-1  ...   predesigned_gene    True
    #1  281474976711278  mpafcjmo-1  ...   predesigned_gene    True
    #2  281474976711281  nljdpiah-1  ...   predesigned_gene    True
    # ...
    # transcripts rows=1985 cols=13
    #print(f"[xenium2scs] transcripts rows={len(transcripts)} cols={len(transcripts.columns)}")

    gene_col = _pick_column(transcripts, ["feature_name", "gene", "gene_id", "geneID"])
    x_col    = _pick_column(transcripts, ["x_location", "x", "x_global_px", "x_centroid"])
    y_col    = _pick_column(transcripts, ["y_location", "y", "y_global_px", "y_centroid"])
    count_col = _pick_column(transcripts, ["counts", "count", "n_counts"], required=False)

    table = transcripts[[gene_col, x_col, y_col]].copy()
    table = table.dropna(subset=[gene_col, x_col, y_col])
    #print(table.head())
    #  feature_name  x_location  y_location
    #0        Defa5  409.250000  301.453125
    #1        Defa5  409.843750  304.281250
    #2        Defa5  418.171875  316.671875
    #3        Defa5  395.000000  290.468750
    #4        Defa5  408.031250  303.218750

    #plot_hexbin_density(table, 
    #                    output_path="/home/katwre/projects/spatialxe_fork/spatialxe/DATA/table_pixel.png")

    # Convert micron coordinates → full-resolution pixel coordinates.
    # Xenium: x_location is along image width (columns), y_location along height (rows).
    table["row_px"]    = (table[y_col].astype(float) / pixel_size).round().astype(int)
    table["column_px"] = (table[x_col].astype(float) / pixel_size).round().astype(int)

    # Optionally merge into user-specified bins (bin_size in pixels, default 1 = no binning).
    if bin_size > 1:
        table["row"]    = (table["row_px"]    / bin_size).astype(int)
        table["column"] = (table["column_px"] / bin_size).astype(int)
    else:
        # Zero-base pixel coordinates so the BGI file starts at (0, 0).
        r0 = table["row_px"].min()
        c0 = table["column_px"].min()
        table["row"]    = table["row_px"]    - r0
        table["column"] = table["column_px"] - c0

    if count_col is None:
        table["counts"] = 1
    else:
        table["counts"] = transcripts.loc[table.index, count_col].fillna(1).astype(int)

    table = table.rename(columns={gene_col: "geneID"})[["geneID", "row", "column", "counts"]]
    table = table.groupby(["geneID", "row", "column"], as_index=False)["counts"].sum()

    #_plot_density_map(table,
    #                 output_path="/home/katwre/projects/spatialxe_fork/spatialxe/DATA/table_spot.png")

    out_tsv = Path(output_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_tsv, sep="\t", index=False)

    # ── Morphology image ────────────────────────────────────────────────────────
    # Load and collapse to 2D (max projection across z/channels).
    image = tifffile.imread(morphology_image_path)
    image = np.squeeze(np.asarray(image))
    if image.ndim == 2:
        image2d = image
    elif image.ndim >= 3:
        h, w = image.shape[-2], image.shape[-1]
        image2d = image.reshape((-1, h, w)).max(axis=0)
    else:
        raise ValueError(f"Unsupported morphology image shape: {image.shape}")

    # contrast stretch for visibility
    p1, p99 = np.percentile(image2d, [1, 99])
    disp = np.clip(image2d, p1, p99)

    import matplotlib.pyplot as plt
    plt.figure(figsize=(8, 8))
    plt.imshow(disp, cmap="gray")
    plt.title(f"Morphology 2D (shape={image2d.shape})")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig("/home/katwre/projects/spatialxe_fork/spatialxe/DATA/morphology.png", 
                dpi=250)
    plt.close()

    # Crop to the pixel ROI covered by transcripts.
    # Derive absolute pixel bounds directly from physical coords in the parquet.
    r_min_abs = int(round(float(transcripts[y_col].min()) / pixel_size))
    r_max_abs = int(round(float(transcripts[y_col].max()) / pixel_size))
    c_min_abs = int(round(float(transcripts[x_col].min()) / pixel_size))
    c_max_abs = int(round(float(transcripts[x_col].max()) / pixel_size))

    # Clamp to image bounds.
    H, W = image2d.shape
    r_min_abs = max(0, r_min_abs)
    r_max_abs = min(H - 1, r_max_abs)
    c_min_abs = max(0, c_min_abs)
    c_max_abs = min(W - 1, c_max_abs)

    cropped = image2d[r_min_abs:r_max_abs + 1, c_min_abs:c_max_abs + 1]

    out_morph2d = Path(output_morph2d_tif)
    out_morph2d.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(out_morph2d, cropped)

    # ── BGI file (SCS/spateo format) ────────────────────────────────────────────
    # spateo read_bgi_agg: x → AnnData dim-0 (height/rows), y → dim-1 (width/cols).
    # Our table["row"] = height direction, table["column"] = width direction.
    bgi = pd.DataFrame({
        "geneID":    table["geneID"],
        "x":         table["row"].astype(int),
        "y":         table["column"].astype(int),
        "MIDCounts": table["counts"].astype(int),
    })

    out_bgi_tsv = Path(output_bgi_tsv)
    out_bgi_tsv.parent.mkdir(parents=True, exist_ok=True)
    bgi.to_csv(out_bgi_tsv, sep="\t", index=False)

    metrics = {
        "n_rows":         int(len(table)),
        "n_unique_genes": int(table["geneID"].nunique()),
        "row_min":        int(table["row"].min())    if len(table) else 0,
        "row_max":        int(table["row"].max())    if len(table) else 0,
        "column_min":     int(table["column"].min()) if len(table) else 0,
        "column_max":     int(table["column"].max()) if len(table) else 0,
        "pixel_size_um":  float(pixel_size),
        "bin_size":       float(bin_size),
        "bin_size_um":    float(bin_size) * float(pixel_size),
        "morph2d_H":      int(cropped.shape[0]),
        "morph2d_W":      int(cropped.shape[1]),
    }

    pd.DataFrame(
        {"metric": list(metrics.keys()), "value": list(metrics.values())}
    ).to_csv(metrics_tsv, sep="\t", index=False)


if __name__ == "__main__":
    transcripts_parquet: str    = "${transcripts_parquet}"
    morphology_image: str       = "${morphology_image}"
    experiment_xenium: str      = "${experiment_xenium}"
    prefix: str                 = "${prefix}"
    bin_size: float             = float("${task.ext.bin_size ?: 5.0}")

    output_tsv        = f"{prefix}/scs_input.tsv"
    output_bgi_tsv    = f"{prefix}/scs_input_bgi.tsv"
    output_morph2d_tif = f"{prefix}/morph2d.tif"
    metrics_tsv       = f"{prefix}/xenium2scs_metrics.tsv"

    convert_xenium_to_scs(
        parquet_path=transcripts_parquet,
        output_tsv=output_tsv,
        output_bgi_tsv=output_bgi_tsv,
        morphology_image_path=morphology_image,
        output_morph2d_tif=output_morph2d_tif,
        metrics_tsv=metrics_tsv,
        experiment_xenium_path=experiment_xenium,
        bin_size=bin_size,
    )

    with open("versions.yml", "w", encoding="utf-8") as fobj:
        fobj.write('"${task.process}":\\n')
        fobj.write('xenium2scs: "1.0.0"\\n')
