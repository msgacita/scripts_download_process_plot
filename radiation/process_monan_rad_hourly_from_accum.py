#!/usr/bin/env python3
"""
Quick check + hourly amounts from accumulated precipitation (MONAN)

What this script does (phase 1):
1) Scans your DIAG files, keeps only the existing hours in December 2018.
2) Loads the accumulated precipitation variable `rainnc(Time, nCells)` lazily.
3) Verifies whether `rainnc` is (near) monotonically non-decreasing over time
   per grid cell (allowing tiny numerical jitter).
4) Computes hourly precipitation amounts by differencing the accumulated field.
   Negative diffs (beyond tolerance) are treated as invalid and set to NaN.
5) Saves:
   - hourly precipitation (mm per hour) as `rain_hourly`
   - a per-timestep fraction of cells with negative diffs as `neg_frac`
   - lon/lat (degrees) for each cell
   to a single NetCDF.

USAGE EXAMPLE
-------------
python getting_MONAN_hourly_precip.py \
  --diag-root /pesq/dados/monan/users/lianet.hernandez/global_clm_2018-2019 \
  --init 2018111500 \
  --out-nc /pesq/dados/monan/users/lianet.hernandez/global_clm_2018-2019/derived_data/dec2018_hourly_rain_from_accum.nc

NOTES
-----
- Assumes one file per output time, 1-hourly, with filenames like:
  MONAN_DIAG_G_MOD_GFS_<INIT>_<YYYYMMDDHH>.00.00.x655362L55.nc
- Works on unstructured grid (dimension: nCells).
- `rainnc` is assumed to be accumulated total precipitation in mm since model start.
- Tolerance for numerical jitter can be adjusted via --neg-eps.

After you verify the hourly product looks good, we can run the diurnal
phase–amplitude script on this output.

Author: Adapted from script by Lianet H. Pardo
"""

import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import datetime

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--diag-root", required=True,
                   help="Folder containing DIAG files (one per time).")
    p.add_argument("--year", type=int, default=2018, help="Year of interest (default: 2018)")
    p.add_argument("--month", type=int, default=12, help="Month of interest (default: 12)")
    p.add_argument("--init", required=True, help="Initialization tag, e.g., 2018111500")
    p.add_argument("--var", default="swdnb",
                   help="Accumulated variable name (default: acswdnb). Starts with ac.")
    p.add_argument("--out-dir", required=True, help="Output directory for NetCDF with hourly data")
    p.add_argument("--chunks", default="auto",
                   help="Dask chunking for nCells (e.g., 200000, or 'auto').")
    return p.parse_args()

TS_RE = re.compile(r".*_(\d{10})\.00\.00\..*\.nc$")  # grabs YYYYMMDDHH before .00.00.

def list_monthly_files(diag_root: Path, init_tag: str, year: int, month: int):
    # Pattern is broad; we filter by regex and by date below.
    pattern = f"MONAN_DIAG_G_MOD_GFS_{init_tag}_*.00.00.*.nc"
    files = sorted(diag_root.glob(pattern))
    out = []
    for f in files:
        m = TS_RE.match(f.name)
        if not m:
            continue
        ts = m.group(1)  # YYYYMMDDHH
        # Keep only files in the specified month 
        if ts.startswith(f"{year}{month:02d}"):
            out.append((pd.to_datetime(ts, format="%Y%m%d%H"), f))
    out.sort(key=lambda x: x[0])
    return out

def build_time_index(files):
    if not files:
        raise FileNotFoundError("No DIAG files found for the specified initialization tag.")
    times = [t for t, _ in files]
    paths = [p for _, p in files]
    idx = pd.DatetimeIndex(times)
    return idx, paths

def load_accum(paths, var, chunks, lat_name="latCell", lon_name="lonCell"):
    """
    Open each file individually and concatenate along Time (lazy). Also read
    lat/lon *once* from the first file, convert from radians to degrees, and
    wrap lon to [-180, 180).

    Returns
    -------
    da  : xr.DataArray  [Time, nCells]   (concatenated variable)
    lat : xr.DataArray  [nCells]         (degrees)
    lon : xr.DataArray  [nCells]         (degrees, wrapped)
    """
    # chunks handling
    if chunks == "auto":
        chunks_dict = {"Time": 1}
    else:
        try:
            n = int(chunks)
            chunks_dict = {"nCells": n, "Time": 1}
        except Exception:
            chunks_dict = {"Time": 1}

    das = []
    lat = None
    lon = None

    for i, p in enumerate(paths, 1):
        print(f"Reading file {p}")
        with xr.open_dataset(
            str(p),
            engine="netcdf4",
            chunks=chunks_dict,
            decode_times=True,
        ) as ds_one:
            if var not in ds_one:
                raise KeyError(f"Variable '{var}' not found in {p}")

            # Read lat/lon once (from the first file we open)
            if (lat is None) or (lon is None):
                if (lat_name not in ds_one) or (lon_name not in ds_one):
                    raise KeyError(f"Expected '{lat_name}' and '{lon_name}' in {p}")

                lat_rad = ds_one[lat_name].astype("float64").load()
                lon_rad = ds_one[lon_name].astype("float64").load()

                # Convert radians -> degrees; wrap lon to [-180, 180)
                lat_deg = np.degrees(lat_rad)
                lon_deg = (np.degrees(lon_rad) + 180.0) % 360.0 - 180.0

                # Keep original 1D dim (usually 'nCells')
                lat = xr.DataArray(lat_deg.values, dims=lat_rad.dims, coords=lat_rad.coords, name="lat")
                lon = xr.DataArray(lon_deg.values, dims=lon_rad.dims, coords=lon_rad.coords, name="lon")

            da_one = ds_one[var]

            # sanity check: Time present and length 1
            if "Time" not in da_one.dims:
                raise KeyError(f"'Time' dim not found in {p} for variable '{var}'")

            das.append(da_one)

    # Concatenate lazily along Time
    da = xr.concat(das, dim="Time", compat="equals", join="exact")

    return da, lat, lon

def compute_hourly_from_accum(accum: xr.DataArray, var: str):
    # Sort by time just in case
    accum = accum.sortby("Time")
    # Ensure units info
    if "units" in accum.attrs:
        units = accum.attrs["units"]
    else:
        units = "W/m^{-2}"
    # Hourly diffs
    diff = accum.diff("Time")  # (Time-1, nCells), in mm
    # Clean diffs: values in [-neg_eps, 0) -> 0; values < -neg_eps -> NaN
    diff_clean = diff/3600
    diff_clean.attrs.update(long_name=f"Hourly {var.removeprefix('ac')} from diff({units} accum)", units="W/m^{-2}")

    # Build an aligned Time coordinate for hourly values (the diff is between t[i-1] -> t[i], we assign to t[i])
    hourly = diff_clean
    hourly = hourly.assign_coords(Time=accum["Time"].isel(Time=slice(1, None)))
    # Attach neg_frac as a DataArray on the same Time

    # Simple summary (print to stdout when executed directly)
    total_steps = hourly.sizes["Time"]
    print(total_steps)
    return hourly

def main():
    args = parse_args()
    diag_root = Path(args.diag_root)
    month_abbr = datetime.datetime.strptime(str(args.month), "%m").strftime("%b")

    # 1) Find existing December files
    monthly_files = list_monthly_files(diag_root, args.init, args.year, args.month)
    if not monthly_files:
        raise SystemExit(f"No files found for {args.year}-{args.month:02d} with the expected pattern.")

    # 2) Build time index
    t_idx, paths = build_time_index(monthly_files)
    print(f"[INFO] Found {len(paths)} files (up to {t_idx[-1]}) for {args.year}-{args.month:02d}")

    # 3) Load accumulated precipitation lazily
    accum, lat, lon = load_accum(paths, args.var, args.chunks)
   
    # (optional) ensure lat/lon attrs are CF-ish
    lat = lat.assign_attrs(units="degrees_north", long_name="latitude")
    lon = lon.assign_attrs(units="degrees_east",  long_name="longitude")

    # Attach our parsed Time index as coordinate to be safe
    accum = accum.assign_coords(Time=xr.DataArray(t_idx, dims=("Time",)))
    accum.name = args.var

    # 4) Compute hourly from accumulated with QC
    hourly = compute_hourly_from_accum(accum, args.var)

    # 6) Save to NetCDF
    ds_out = xr.Dataset(
        data_vars={
            f"{args.var.removeprefix('ac')}": hourly.astype("float32"),
        },
        coords={
            "Time": hourly["Time"],
            "lat": (("nCells",), lat.values, lat.attrs),
            "lon": (("nCells",), lon.values, lon.attrs),
        },
    )
    ds_out[args.var.removeprefix("ac")].attrs.update(
        long_name=f"Hourly {args.var.removeprefix('ac')} amount derived from accumulated {args.var}",
        units="W/m^{-2}",
    )

    out_path = Path(args.out_dir) / f"{month_abbr}{year}_hourly_{args.var.removeprefix('ac')}.nc"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds_out.to_netcdf(out_path, engine="h5netcdf", mode="w")
    print(f"[OK] Wrote {out_path}")

    # Calculate monthly hourly mean (diurnal cycle)
    print(f"\n[INFO] Calculating monthly hourly mean (diurnal cycle) for {args.year}-{args.month:02d}...")
    times = pd.DatetimeIndex(hourly["Time"].values)
    hours = times.hour
    hourly_with_hour = hourly.assign_coords(hour=("Time", hours))
    
    # Disable flox to avoid TypeAlias import error in Python 3.9
    with xr.set_options(use_flox=False):
        monthly_hourly_mean = hourly_with_hour.groupby("hour").mean(dim="Time")
    
    print(f"[INFO] Monthly hourly mean shape: {monthly_hourly_mean.shape}")
    print(f"[INFO] Hours: {list(monthly_hourly_mean['hour'].values)}")
    
    # Calculate local time hour for each cell based on longitude
    # Local time = UTC time + longitude/15 (converted to hours)
    print("[INFO] Calculating local time hour for each cell...")
    lon_offset_hours = lon.values / 15.0  # Longitude offset in hours
    
    # For each cell, calculate the local hour for each UTC hour
    # Result shape: (hour, nCells)
    utc_hours = monthly_hourly_mean["hour"].values  # 0-23
    local_hour = np.zeros((len(utc_hours), len(lon_offset_hours)))
    
    for i, utc_h in enumerate(utc_hours):
        local_hour[i, :] = (utc_h + lon_offset_hours) % 24
    
    # Save monthly hourly mean to separate file
    var_name = args.var.removeprefix('ac')
    out_mhourly = out_path.parent / f"{out_path.stem}_MHourly.nc"
    
    ds_mhourly = xr.Dataset(
        data_vars={
            f"{var_name}_MHourly": monthly_hourly_mean.astype("float32"),
            "local_hour": (("hour", "nCells"), local_hour.astype("float32"))
        },
        coords={
            "hour": monthly_hourly_mean["hour"],
            "lat": (("nCells",), lat.values, lat.attrs),
            "lon": (("nCells",), lon.values, lon.attrs),
        },
    )
    ds_mhourly["local_hour"].attrs.update(
        long_name="Local Solar Time Hour",
        units="hours",
        description="Local time hour (0-24) calculated from UTC hour + longitude/15"
    )
    ds_mhourly[f"{var_name}_MHourly"].attrs.update(
        long_name=f"Monthly Hourly Mean {var_name}",
        units="W m^{-2}",
        description=f"Diurnal cycle: mean for each hour (0-23 UTC) averaged over {args.year}-{args.month}"
    )
    
    ds_mhourly.to_netcdf(out_mhourly, engine="h5netcdf", mode="w")
    print(f"[OK] Wrote monthly hourly mean to {out_mhourly}")

    

if __name__ == "__main__":
    main()

