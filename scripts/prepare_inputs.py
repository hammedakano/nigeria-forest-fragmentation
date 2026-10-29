"""
Prepare GEE upload files from the original QGIS reference libraries.

Inputs (not redistributed in this repository, see data/README.md):
  - 16 GeoPackages of photo-interpreted reference points, one per
    reserve and epoch (<RESERVE>-18-*.gpkg for 2019, <RESERVE>-23-*.gpkg for 2024)
  - 8 reserve boundary shapefiles exported from the WDPA (December 2025 release)

Outputs:
  - data/reference/original_reference_points.csv   (upload to GEE as a table asset)
  - gee_upload/reserves.zip                        (zipped shapefile, upload to GEE)

Usage:
  python scripts/prepare_inputs.py --points <dir with gpkg> --boundaries <dir with shp> --out <repo root>

Requires geopandas >= 0.14.
"""
import argparse
import glob
import os
import zipfile

import geopandas as gpd
import pandas as pd

# Prefix of the original file names -> reserve identifier used everywhere in the workflow
POINT_PREFIX = {
    "ADORU": "Adoru",
    "AKANTO": "Akanto",
    "GAMBARI": "Gambari",
    "MARADUN": "Maradun",
    "MAYO": "Mayo_Ndaga",
    "STUBBS-CREEK": "Stubbs_Creek",
    "UPPER-ORASHI": "Upper_Orashi",
    "YO": "Yo",
}
BOUNDARY_FILE = {
    "Adoru": "ADORU.shp",
    "Akanto": "AKANTO.shp",
    "Gambari": "GAMBARI.shp",
    "Maradun": "MARADUN.shp",
    "Mayo_Ndaga": "MAYO_NDAGA.shp",
    "Stubbs_Creek": "STUBBS_CREEK.shp",
    "Upper_Orashi": "UPPER_ORASHI_FORESTS.shp",
    "Yo": "YO.shp",
}
BIOME = {
    "Adoru": "Savannah", "Akanto": "Rainforest", "Gambari": "Rainforest",
    "Maradun": "Savannah", "Mayo_Ndaga": "Savannah", "Stubbs_Creek": "Mangrove/swamp",
    "Upper_Orashi": "Mangrove/swamp", "Yo": "Savannah",
}
CLASS_CODE = {"No vegetation": 0, "Low vegetation": 1, "High vegetation": 2}
# Correct the reversed High/Low labels of the source files (see load_points)
SWAP_HIGH_LOW = True


def load_points(point_dir):
    frames = []
    for f in sorted(glob.glob(os.path.join(point_dir, "*.gpkg"))):
        base = os.path.basename(f).upper()
        prefix = [k for k in POINT_PREFIX if base.startswith(k)]
        if not prefix:
            continue
        # The longest matching prefix wins (YO must not capture other names)
        prefix = max(prefix, key=len)
        epoch = 2019 if "-18" in base or "_18" in base else 2024
        g = gpd.read_file(f).to_crs(4326)
        g["reserve"] = POINT_PREFIX[prefix]
        g["year"] = epoch
        g["source_file"] = os.path.basename(f)
        frames.append(g)
    pts = pd.concat(frames, ignore_index=True)

    # Harmonise labels. The original files mix "High vegetation" and "High Vegetation"
    # and contain three unlabelled points, which are removed here.
    pts["class_name"] = (
        pts["name"].fillna("").str.strip().str.lower()
        .map({k.lower(): k for k in CLASS_CODE})
    )
    n_unlabelled = int(pts["class_name"].isna().sum())
    pts = pts[pts["class_name"].notna()].copy()
    # The High and Low labels in the source libraries are reversed. Points named
    # "Low vegetation" fall on closed canopy and points named "High vegetation"
    # on open farmland and fallow mosaic. This was confirmed by visual checks over
    # high-resolution imagery and by Dynamic World, ESA WorldCover, ETH canopy
    # height and Sentinel-2 NDVI at every point (agreement with Dynamic World 20 %
    # as named, 70 % after swapping). The original name is kept in name_original.
    pts["name_original"] = pts["class_name"]
    if SWAP_HIGH_LOW:
        pts["class_name"] = pts["class_name"].map(
            {"High vegetation": "Low vegetation", "Low vegetation": "High vegetation",
             "No vegetation": "No vegetation"})
    pts["class"] = pts["class_name"].map(CLASS_CODE).astype(int)
    pts["longitude"] = pts.geometry.x.round(7)
    pts["latitude"] = pts.geometry.y.round(7)
    n_dup = int(pts.duplicated(["reserve", "year", "longitude", "latitude"]).sum())
    pts = pts.drop_duplicates(["reserve", "year", "longitude", "latitude"])
    pts = pts.sort_values(["reserve", "year", "latitude", "longitude"]).reset_index(drop=True)
    # Identifier restarts within each reserve-year, e.g. REF_Adoru_2019_001
    pts["ref_id"] = (pts.groupby(["reserve", "year"]).cumcount() + 1).map("{:03d}".format)
    pts["ref_id"] = "REF_" + pts.reserve + "_" + pts.year.astype(str) + "_" + pts["ref_id"]
    print(f"Removed {n_unlabelled} unlabelled and {n_dup} duplicate points; kept {len(pts)}")
    return pts


def load_boundaries(boundary_dir):
    rows = []
    for reserve, fname in BOUNDARY_FILE.items():
        g = gpd.read_file(os.path.join(boundary_dir, fname))
        if g.crs is None:
            g = g.set_crs(4326)
        g = g.to_crs(4326).dissolve()
        rows.append({
            "reserve": reserve,
            "biome": BIOME[reserve],
            "wdpaid": int(g["WDPAID"].iloc[0]),
            "wdpa_name": g["NAME"].iloc[0],
            "geometry": g.geometry.iloc[0],
        })
    return gpd.GeoDataFrame(rows, crs=4326)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", required=True)
    ap.add_argument("--boundaries", required=True)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()

    pts = load_points(a.points)
    bnd = load_boundaries(a.boundaries)

    # Every point must fall inside (or within ~50 m of) its own reserve
    for r, grp in pts.groupby("reserve"):
        poly = bnd.loc[bnd.reserve == r].geometry.iloc[0].buffer(0.0005)
        outside = (~grp.geometry.within(poly)).sum()
        assert outside == 0, f"{outside} points of {r} fall outside the reserve"

    ref_dir = os.path.join(a.out, "data", "reference")
    os.makedirs(ref_dir, exist_ok=True)
    cols = ["ref_id", "reserve", "year", "class", "class_name", "name_original",
            "longitude", "latitude", "source_file"]
    pts[cols].to_csv(os.path.join(ref_dir, "original_reference_points.csv"), index=False)

    summary = pts.pivot_table(index=["reserve", "year"], columns="class_name",
                              values="class", aggfunc="size", fill_value=0)
    summary["total"] = summary.sum(axis=1)
    summary.to_csv(os.path.join(ref_dir, "original_reference_points_summary.csv"))
    print(summary)

    up = os.path.join(a.out, "gee_upload")
    os.makedirs(up, exist_ok=True)
    shp = os.path.join(up, "reserves.shp")
    bnd.to_file(shp)
    with zipfile.ZipFile(os.path.join(up, "reserves.zip"), "w", zipfile.ZIP_DEFLATED) as z:
        for ext in ("shp", "shx", "dbf", "prj", "cpg"):
            p = os.path.join(up, f"reserves.{ext}")
            if os.path.exists(p):
                z.write(p, f"reserves.{ext}")
                os.remove(p)
    bnd.drop(columns="geometry").to_csv(os.path.join(a.out, "data", "reserves.csv"), index=False)


if __name__ == "__main__":
    main()
