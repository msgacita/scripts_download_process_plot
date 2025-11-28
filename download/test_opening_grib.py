import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import sys
from pathlib import Path
import xarray as xr
import cfgrib

def file_diagnostics(path: Path):
    """Return quick diagnostics about a file to help debug EOF/truncation issues."""
    info = {}
    try:
        stat = path.stat()
        info["size_bytes"] = stat.st_size
    except Exception as e:
        info["stat_error"] = str(e)

    try:
        with path.open("rb") as f:
            head = f.read(256)
        info["head_bytes"] = head[:64]
        info["head_hex"] = head[:32].hex()
    except Exception as e:
        info["read_error"] = str(e)

    return info

def try_cfgrib(path: str):
    """Try opening with xarray+cfgrib. Return dataset or raise exception."""
    ds = xr.open_dataset(path, engine="cfgrib")
    return ds

def try_pygrib(path: str):
    """Try inspecting the GRIB file with pygrib (if installed). Returns list of msgs or raises."""
    try:
        import pygrib
    except Exception:
        pygrib = None
    if pygrib is None:
        raise RuntimeError("pygrib not installed")
    grbs = pygrib.open(path)
    msgs = []
    try:
        for i, g in enumerate(grbs):
            # collect a small summary of the first few messages
            msgs.append({
                "msg_no": i + 1,
                "name": getattr(g, "name", None),
                "shortName": getattr(g, "shortName", None),
                "typeOfLevel": getattr(g, "typeOfLevel", None),
                "level": getattr(g, "level", None),
            })
            if i >= 9:
                break
    finally:
        grbs.close()
    return msgs

def main(argv=None):
    argv = argv or sys.argv[1:]
    var_list = [
    # "divergence",
    # "fraction_of_cloud_cover",
    # "geopotential",
    # "ozone_mass_mixing_ratio",
    # "potential_vorticity",
    # "relative_humidity",
    # "specific_cloud_ice_water_content",
    # "specific_cloud_liquid_water_content",
    # "specific_humidity",
    # "specific_rain_water_content",
    # "specific_snow_water_content",
    # "temperature",
    # "u_component_of_wind",
    # "v_component_of_wind",
    # "vertical_velocity",
    # "vorticity"
    'clear_sky_direct_solar_radiation_at_surface'
    ]
    year='2018'
    month='12'
    for var in var_list:
        print(var)
        file_arg = argv[0] if len(argv) >= 1 else  "/pesq/dados/monan/users/madeleine.gacita/global_data/era5/single_levels/"+year+"/"+month+"/"+var+"_MHour.grib"
        path = Path(file_arg)

        if not path.exists():
            print(f"File does not exist: {path!s}")
            # return 1

        print(f"Diagnostics for: {path}")
        diag = file_diagnostics(path)
        for k, v in diag.items():
            if k == "head_bytes":
                print(f"  {k}: <{len(v)} bytes>")
            else:
                print(f"  {k}: {v}")

        # First try xarray+cfgrib (the common path)
        try:
            ds = try_cfgrib(str(path))
            print("Successfully opened with cfgrib/xarray")
            print(ds)
            print(ds['cdir'].shape)
            print(ds['cdir'].time)
            print(ds['cdir'].step)
            ds.close()
            # return 0
        except Exception as e:
            msg = f"{type(e).__name__}: {e}"
            print(f"cfgrib/xarray failed: {msg}")

            # # Common cfgrib/ecCodes failure: stale or corrupt index (.idx) file
            # # Example index name: 2019/01/vorticity.grib.5b7b6.idx
            # if "Can't read index file" in msg or ".idx" in msg:
            #     print("Detected an .idx-related error from cfgrib. Searching for matching .idx files to remove and retry...")
            #     try:
            #         removed = []
            #         pattern = path.name + ".*.idx"
            #         for idx in path.parent.glob(pattern):
            #             try:
            #                 idx.unlink()
            #                 removed.append(str(idx))
            #             except Exception as ex:
            #                 print(f"  Could not remove index file {idx}: {ex}")

            #         if removed:
            #             print("Removed index files:")
            #             for r in removed:
            #                 print(f"  - {r}")
            #             print("Retrying cfgrib open after removing .idx files...")
            #             try:
            #                 ds = try_cfgrib(str(path))
            #                 print("Successfully opened with cfgrib/xarray after removing .idx files")
            #                 print(ds)
            #                 ds.close()
            #                 return 0
            #             except Exception as e2:
            #                 print(f"Retry after removing .idx files still failed: {type(e2).__name__}: {e2}")
            #         else:
            #             print(f"No .idx files matching pattern {pattern} were found next to {path}")
            #     except Exception as ex:
            #         print(f"Error while attempting to remove .idx files: {ex}")

            #     # As an additional fallback, try opening cfgrib with a transient/in-memory index (backend option)
            #     try:
            #         print("Attempting to open with cfgrib using an in-memory/temporary index (backend_kwargs={'indexpath': ''})...")
            #         ds = xr.open_dataset(str(path), engine="cfgrib", backend_kwargs={"indexpath": ""})
            #         print("Successfully opened with cfgrib using temporary index")
            #         print(ds)
            #         ds.close()
            #         return 0
            #     except Exception as e3:
            #         print(f"Opening with temporary index also failed: {type(e3).__name__}: {e3}")

        # # If cfgrib failed, try pygrib as a lighter-weight reader to check file integrity
        # try:
        #     msgs = try_pygrib(str(path))
        #     print(f"pygrib was able to read {len(msgs)} messages (showing up to 10):")
        #     for m in msgs:
        #         print(f"  - msg {m['msg_no']}: {m['shortName']} ({m['name']}) level={m['level']}")
        #     return 0
        # except Exception as e:
        #     print(f"pygrib inspection failed: {type(e).__name__}: {e}")

if __name__ == "__main__":
    
    main()
