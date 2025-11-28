import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path
import xarray as xr
import cfgrib
import argparse
from copy import deepcopy
import matplotlib.pyplot as plt

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--year", required=True,
                   help="Year of the ERA5 files.")
    p.add_argument("--month", required=True,
                   help="Month of the ERA5 files.")
    p.add_argument("--var", required=True, default='ssrd',
                   help="Variable name in the ERA5 files.")
    p.add_argument("--era5-root", required=False, default='/pesq/dados/monan/users/madeleine.gacita/global_data/era5/single_levels/',
                   help="Folder containing ERA5 single level files (one per variable and month).")
    args = p.parse_args()
    
    # Validate variable choice
    if args.var not in ['ssrd', 'ssrdc', 'strd', 'ssr', 'str', 'sshf', 'slhf','ttr','ttrc','tisr']:
        print(f"Error: Variable {args.var} is not supported. Choose from 'ssrd', 'strd', 'ssr', 'str', 'sshf', 'slhf','ttr','ttrc','tisr'.")
        sys.exit(1)
    
    return args

def flatten_dataset(ds, var, time_selection=[]):
    """ ECMWF code, see more at https://confluence.ecmwf.int/pages/viewpage.action?pageId=536218894
    Rearange the data so we flatten time samples.
    E.g. from
    {'time': 2, 'step': 24, 'latitude': 1801, 'longitude': 3600}
    to
    {'time': 48, 'latitude': 1801, 'longitude': 3600}
    """
    total_samples = ds.valid_time.shape[0]*ds.valid_time.shape[1]
     
    reshaped_var = ds[var].values.reshape((total_samples, ds[var].shape[2], ds[var].shape[3]))
 
    ds_new = xr.Dataset(
        data_vars={var: (["time", "latitude", "longitude"], reshaped_var)},
        coords={
            "time": ("time", ds.valid_time.values.reshape((total_samples))),
            "latitude": ("latitude", ds.latitude.values),
            "longitude": ("longitude", ds.longitude.values),
        },
    )
 
    ds_new.attrs = ds.attrs
    ds_new[var].attrs = ds[var].attrs
    ds_new.longitude.attrs = ds.longitude.attrs
    ds_new.latitude.attrs = ds.latitude.attrs
    ds_new.time.attrs = ds.valid_time.attrs
    if time_selection is not None:
        ds_new = ds_new.isel(time=slice(*time_selection))
 
    # Remove additional variables, same as in ARCO hourly ERA5 on single levels
    for drop_var in ["number", "surface"]:
        if drop_var in ds_new:
            ds_new = ds_new.drop_vars([drop_var])
 
    return ds_new

def apply_lon_lat_conventions(ds):
    # Renames
    if "lon" in ds.dims:
        ds = ds.rename({"lon": "longitude"})
    if "lat" in ds.dims:
        ds = ds.rename({"lat": "latitude"})
 
    # Flip latitudes (ensure they are monotonic increasing)
    if "latitude" in ds.dims:
        lats = ds["latitude"]
        if len(lats) > 1 and lats[0] > lats[-1]:
            ds = ds.reindex(latitude=ds.latitude[::-1])
 
    # Convert longitude to [-180, 180[
    if "longitude" in ds.dims and ds["longitude"].max() > 180:
        lons = ds["longitude"]
        lons_attrs = lons.attrs
        new_lons = np.concatenate([lons[lons >= 180], lons[lons < 180]])
        ds = ds.reindex(longitude=new_lons)
        ds = ds.assign_coords(longitude=(((ds["longitude"] + 180) % 360) - 180))
        ds["longitude"].attrs = lons_attrs
    return ds
 


def deaccumulate_monthly(ds0, variable):
    '''  ECMWF code, see more at https://confluence.ecmwf.int/pages/viewpage.action?pageId=536218894
     De-accumulate the variable of one Dataset
    Args:
        ds0: Date 0 (e.g. 1970-01-02) dataset.
        variable: The dataset variable
     
    Returns:
        A de-accumulated Dataset.'''

    ds0_flat = flatten_dataset(ds0, variable, None)
    ds_deaccumulated = deepcopy(ds0_flat)
    print(ds0_flat.time)
    input()
     
    # Array indexes
    # ds-1, date-1, step23
    i_ds_1_d_1_s23 = 46
    # ds0, date-1, step24
    i_ds0_d_1_s24 = 23
    # ds0, date0, step1
    i_ds0_d0_s1 = 24
     
    # # 00:00    <date-1, step24> - <date-1, step23>
    ds_deaccumulated[variable][i_ds0_d_1_s24] -= ds_1_flat[variable][i_ds_1_d_1_s23]
     
    # Nothing to do for 01:00    <date0, step1>,
     
    # From
    #    02:00    <date0, step2> - <date0, step1>,
    # to
    #    23:00    <date0, step23> - <date0, step22>
    for i, t in enumerate(range(ds_deaccumulated.time.shape[0])):
        print(t)
        if i > i_ds0_d0_s1:
            ds_deaccumulated[variable][t] -= ds0_flat[variable][t - 1]
     
    # Slice to the times of interest (24 samples)
    ds_deaccumulated = ds_deaccumulated.isel(time=slice(i_ds0_d_1_s24, i_ds_1_d_1_s23 + 1))
    print(ds_deaccumulated)
     
    return ds_deaccumulated
 



def accumulate_monthly(ds0_deaccumulated, variable, ds_1_accumulated=None, ds_1_deaccumulated=None):
    """
     ECMWF code, see more at https://confluence.ecmwf.int/pages/viewpage.action?pageId=536218894
    Accumulate a variable over each time sample (revert de de-accumulation to obtain the original values).
     
    Args:
        ds0_deaccumulated: Date 0 (e.g. 1970-01-02) de-accumulated
        variable: Dataset variable
        ds_1_accumulated: Date -1 (e.g. 1970-01-01) accumulated.  Expected time and steps dimensions.
                          Needed either this or ds_1_deaccumulated
        ds_1_deaccumulated: Date -1 (e.g. 1970-01-01) de-accumulated. Expected time dimension.
                            Needed either this or ds_1_accumulated
         
    Returns:
        Date 0 xarray.Dataset
    """
    ds0_accumulated = deepcopy(ds0_deaccumulated)
     
    if ds_1_accumulated:
        ds0_accumulated[variable][0] += ds_1_accumulated[variable].values[1][22]
    elif ds_1_deaccumulated:
        for i in range(1, 24):
            ds0_accumulated[variable][0] += ds_1_deaccumulated[variable].values[i]
 
    for step in range(2, 24):
        tmp = deepcopy(ds0_deaccumulated[variable].values[1])
        for acc_step in range(2, step+1):
            tmp += ds0_deaccumulated[variable].values[acc_step]
        ds0_accumulated[variable][step] = tmp
     
    return ds0_accumulated
 


def main():
    args = parse_args()
    era5_root = Path(args.era5_root)
    var_dict = {
        'ssrd': 'surface_solar_radiation_downwards',
        'ssrdc': 'surface_solar_radiation_downward_clear_sky',
        'strd': 'surface_thermal_radiation_downwards',
        'ssr': 'surface_net_solar_radiation',
        'str': 'surface_net_thermal_radiation',
        'sshf': 'surface_sensible_heat_flux',
        'slhf': 'surface_latent_heat_flux', 
        'ttr': 'top_net_thermal_radiation',
        'ttrc': 'top_net_thermal_radiation_clear_sky',
        'tisr': 'toa_incident_solar_radiation'
    }
    # Example: Load variable for a given month
    var_file = era5_root / f"{args.year}/{args.month}/{var_dict[args.var]}_MHour.nc"
    print(var_file)
    if not var_file.exists():
        print(f"File not found: {var_file}")
        sys.exit(1)
    print(f"Loading {args.var} data from {var_file}...")
    
    ds0 = xr.open_dataset(var_file)
    ds = apply_lon_lat_conventions(ds0)

    print(ds)
    print(ds[args.var].shape)
    print(ds['valid_time'])

    print(ds[args.var].attrs)
    flux_var_name = f'{args.var}_flux'
    ds[flux_var_name] = ds[args.var]/3600
    ds[flux_var_name].attrs = {'units': 'W m**-2'}
    ds[flux_var_name].attrs['long_name'] = var_dict[args.var]
    ds[flux_var_name].attrs['description'] = f'Instantaneous flux calculated from accumulated {args.var.upper()}'
    
    ax=plt.figure(figsize=(8, 6))
    ds[flux_var_name][14].plot(cmap='RdYlBu')
    plt.savefig(f'{args.var}_monthly_mean_15UTC_{args.year}{args.month}.png', dpi=150, bbox_inches='tight')
    plt.close()
    # Calculate mean over the 24 hours (valid_time dimension)
    print(f"\nCalculating mean over 24 hours...")
    flux_mean_24h = ds[flux_var_name].mean(dim='valid_time')
    print(f"24-hour mean flux shape: {flux_mean_24h.shape}")
    print(f"24-hour mean flux value (global): {float(flux_mean_24h.mean()):.2f} W/m²")
    
    # Add mean as a new variable to the dataset
    mean_var_name = f'{args.var}_flux_24h_mean'
    ds[mean_var_name] = flux_mean_24h
    ds[mean_var_name].attrs = {
        'units': 'W m**-2',
        'long_name': f'{var_dict[args.var]} - 24h Mean',
        'description': f'24-hour mean of instantaneous flux from {args.var.upper()}'
    }
    plt.close()
    # Plot the 24-hour mean
    ax=plt.figure(figsize=(8, 6))
    ds[mean_var_name].plot(cmap='RdYlBu', ax=ax)
    plt.title(f'{args.var.upper()} 24-hour Mean Flux - {args.year}/{args.month}')
    plt.savefig(f'{args.var}_monthly_mean_15UTC_{args.year}{args.month}.png', dpi=150, bbox_inches='tight')

    # if ((ds0[args.var].GRIB_stepType=='accum') & (ds0[args.var].GRIB_units=='J m**-2')):
    #     # print("Variable {args.var} is accumulated energy ")
    #     # print(ds0.time.values)

        # ds_mhourly_mean = deaccumulate_monthly(ds0, args.var)
        # print(ds_mhourly_mean)
    #     # # ds_accumulated = accumulate(ds_deaccumulated, args.var)


    # # if ((ds[args.var].GRIB_stepType=='accum') & (ds[args.var].GRIB_units=='J m**-2')):
    # #     # Variable is accumulated energy
    #     # Converting to instantaneous flux (W/m²) from differences and timestep intervals
    #     print(f"\nCalculating instantaneous flux from accumulated {args.var.upper()}...")
    #     print(f"Original dataset shape: {ds0[args.var].shape}")
    #     print(f"Dimensions: {ds0[args.var].dims}")
        
    #     # Check if step is a dimension or scalar coordinate
    #     if 'step' in ds0.dims:
    #         print(f"Coordinates: time={len(ds0['time'])}, step={len(ds0['step'])}")
    #     elif 'step' in ds0.coords:
    #         print(f"Coordinates: time={len(ds0['time'])}, step={ds0['step'].values} (scalar)")
    #     else:
    #         print(f"Coordinates: time={len(ds0['time'])}")

    #     # Get step values in seconds for all steps
    #     all_step_seconds = np.array([pd.Timedelta(s).total_seconds() for s in ds0['step'].values])

    #     # Initialize flux array with same shape as original
    #     flux_values = np.zeros_like(ds0[args.var].values)

    #     # First step: accumulated value / step time
    #     flux_values[:, 0, :, :] = ds0[args.var].values[:, 0, :, :] / all_step_seconds[0]

    #     # Remaining steps: difference / time interval
    #     for i in range(1, len(all_step_seconds)):
    #         time_interval = all_step_seconds[i] - all_step_seconds[i-1]
    #         flux_values[:, i, :, :] = (ds0[args.var].values[:, i, :, :] - ds0[args.var].values[:, i-1, :, :]) / time_interval

    #     # Create DataArray with proper coordinates
    #     flux_var_name = f'{args.var}_flux'
    #     ds0[flux_var_name] = xr.DataArray(
    #         flux_values,
    #         dims=ds0[args.var].dims,
    #         coords=ds0[args.var].coords,
    #         attrs={'long_name': f'{args.var.upper()} Flux',
    #                 'units': 'W m-2',
    #                 'description': f'Instantaneous flux calculated from accumulated {args.var.upper()}'}
    #     )

    #     print(f"\nCreated {flux_var_name} variable:")
    #     print(f"  Shape: {ds0[flux_var_name].shape}")

    #     # Verify no NaN values from the calculation
    #     nan_count = np.isnan(ds0[flux_var_name].values).sum()
    #     print(f"  NaN values: {nan_count} ({100*nan_count/ds0[flux_var_name].size:.2f}%)")

    #     # Convert longitude from 0-360 to -180 to 180 format if needed
    #     if 'longitude' in ds0.coords:
    #         lon_values = ds0['longitude'].values
    #         if lon_values.max() > 180:
    #             print("\nConverting longitude from 0-360° to -180-180° format...")
    #             lon_converted = np.where(lon_values > 180, lon_values - 360, lon_values)
    #             ds0['longitude'] = lon_converted
    #             ds0 = ds0.sortby('longitude')
    #             print("Longitude converted and sorted")

    #     # Save full dataset to NetCDF
    #     output_file = era5_root / f"{args.year}/{args.month}/{var_dict[args.var]}_flux.nc"
    #     print(f"\nSaving full dataset to: {output_file}")
    #     ds0.to_netcdf(output_file, engine='netcdf4')
    #     print(f"✅ Successfully saved {flux_var_name}")

    #     # Print sample statistics
    #     print(f"\nSample flux statistics:")
    #     for i in range(min(5, len(ds0['step']))):
    #         step_hours = all_step_seconds[i] / 3600
    #         flux_data = ds0[flux_var_name].isel(step=i)
    #         print(f"  Step {int(step_hours):02d}h: mean={float(flux_data.mean()):.2f} W/m², "
    #                 f"min={float(flux_data.min()):.2f}, max={float(flux_data.max()):.2f}")

    #     # ====================================================================
    #     # CALCULATE HOURLY MONTHLY MEAN (preserves all 24 hours)
    #     # ====================================================================
    #     print("\n" + "="*60)
    #     print("Calculating hourly monthly mean climatology...")
    #     print("="*60)

    #     # Step 1: Calculate actual valid time for each (time, step) combination
    #     print("\nCalculating valid times and extracting hour-of-day...")

    #     time_values = ds0['time'].values
    #     step_values = ds0['step'].values

    #     # Create array of valid times and corresponding hours
    #     valid_times = []
    #     for t in time_values:
    #         for s in step_values:
    #             valid_time = pd.Timestamp(t) + pd.Timedelta(s)
    #             valid_times.append(valid_time)

    #     valid_times = np.array(valid_times)
    #     hours = pd.DatetimeIndex(valid_times).hour

    #     print(f"  Total time points: {len(valid_times)}")
    #     print(f"  Unique hours: {sorted(set(hours))}")

    #     # Step 2: Stack time and step into single dimension
    #     print("\nStacking time and step dimensions...")
    #     ds_stacked = ds0[[flux_var_name]].stack(valid_time=('time', 'step'))

    #     # Step 3: Assign hour coordinate
    #     ds_stacked['hour'] = ('valid_time', hours)
    #     ds_stacked = ds_stacked.set_coords('hour')

    #     # Step 4: Group by hour and calculate mean
    #     print("Grouping by hour-of-day and calculating means...")
    #     with xr.set_options(use_flox=False):
    #         ds_hourly_mean = ds_stacked.groupby('hour').mean(dim='valid_time')

    #     # Update attributes to reflect monthly mean
    #     ds_hourly_mean[flux_var_name].attrs['long_name'] = f'{args.var.upper()} Flux - Monthly Hourly Mean'
    #     ds_hourly_mean[flux_var_name].attrs['units'] = 'W m-2'
    #     ds_hourly_mean[flux_var_name].attrs['description'] = (
    #         f'Hourly monthly means of instantaneous flux calculated from accumulated {args.var.upper()}. '
    #         f'Averaged over all days in {args.year}-{args.month} for each hour of day (UTC).'
    #     )

    #     print(f"\n✅ Hourly climatology complete!")
    #     print(f"  Hours: {len(ds_hourly_mean['hour'])} (0-23)")
    #     print(f"  Variables: {list(ds_hourly_mean.data_vars)}")
    #     print(f"  Shape: {ds_hourly_mean[flux_var_name].shape}")

    #     # Save hourly monthly mean
    #     output_hourly = era5_root / f"{args.year}/{args.month}/{var_dict[args.var]}_flux_MHour.nc"
    #     print(f"\nSaving hourly monthly mean to: {output_hourly}")
    #     ds_hourly_mean.to_netcdf(output_hourly, engine='netcdf4')
    #     print(f"✅ Successfully saved hourly monthly mean")

#     # Print hourly statistics
    # for t in ds0.valid_time:
    #     hour_data = ds0[args.var].sel(valid_time=t)/3600
    #     print(f"  Hour {t.dt.hour.values:02d} UTC: mean={float(hour_data.mean()):.2f} J/m², "
    #             f"min={float(hour_data.min()):.2f}, max={float(hour_data.max()):.2f}")

    #     print(f"\nHourly monthly mean statistics for {args.var}:")
    #     for t in ds0.time:
    #         hour_data = ds_mhourly_mean[args.var].sel(time=t)
    #         print(f"  Hour {t.dt.hour.values:02d} UTC: mean={float(hour_data.mean()):.2f} W/m², "
    #                 f"min={float(hour_data.min()):.2f}, max={float(hour_data.max()):.2f}")

    #     ds_mhourly_mean.close()
    #     ds0.close()
    # #     print("\n" + "="*60)
    # #     print("✅ All processing complete!")
    # #     print("="*60)
    # #     print(f"\nOutput files:")
    # #     print(f"  1. Full dataset:     {output_file}")
    # #     print(f"  2. Hourly mean:      {output_hourly}")
    # if args.var not in ds0.variables:
    #     print(f"Variable {args.var} not found in dataset.")
    #     sys.exit(1)
    # else:
    #     print(f"\nVariable {args.var} is not accumulated energy. No flux calculation performed.")
    #     print("Exiting without changes.")
    #     sys.exit(1)
if __name__ == "__main__":
    main()
