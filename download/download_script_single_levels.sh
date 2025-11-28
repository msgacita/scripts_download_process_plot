#!/bin/bash


# Check for month and year arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <month> <year>"
    exit 1
fi

month="$1"
year="$2"
mkdir -p "/pesq/dados/monan/users/madeleine.gacita/global_data/era5/single_levels/$year/$month"


#List of meteorological variables
variables=(
        "10m_u_component_of_wind"
        "10m_v_component_of_wind"
        "2m_dewpoint_temperature"
        "2m_temperature"
        "mean_sea_level_pressure"
        "sea_surface_temperature"
        "surface_pressure"
        "total_precipitation"
        "skin_temperature"
        "100m_u_component_of_wind"
        "100m_v_component_of_wind"
        "clear_sky_direct_solar_radiation_at_surface"
        "surface_latent_heat_flux"
        "surface_net_solar_radiation"
        "surface_net_solar_radiation_clear_sky"
        "surface_net_thermal_radiation"
        "surface_net_thermal_radiation_clear_sky"
        "surface_sensible_heat_flux"
        "surface_solar_radiation_downward_clear_sky"
        "surface_solar_radiation_downwards"
        "surface_thermal_radiation_downward_clear_sky"
        "surface_thermal_radiation_downwards"
        "toa_incident_solar_radiation"
        "top_net_solar_radiation"
        "top_net_solar_radiation_clear_sky"
        "top_net_thermal_radiation"
        "top_net_thermal_radiation_clear_sky"
        "total_sky_direct_solar_radiation_at_surface"
        "cloud_base_height"
        "high_cloud_cover"
        "low_cloud_cover"
        "medium_cloud_cover"
        "total_cloud_cover"
        "total_column_cloud_ice_water"
        "total_column_cloud_liquid_water"
        "convective_precipitation"
        "convective_rain_rate"
        "large_scale_rain_rate"
        "large_scale_precipitation"
        "large_scale_precipitation_fraction"
        "precipitation_type"
        "total_column_rain_water"
        "convective_snowfall"
        "snowfall"
        "total_column_snow_water"
        "high_vegetation_cover"
        "leaf_area_index_high_vegetation"
        "leaf_area_index_low_vegetation"
        "low_vegetation_cover"
        "type_of_high_vegetation"
        "type_of_low_vegetation"
        "boundary_layer_height"
        "convective_available_potential_energy"
        "convective_inhibition"
        "forecast_albedo"
        "geopotential"
        "land_sea_mask"
        "sea_ice_cover"
        "total_column_water"
        "total_column_water_vapour"
)
# variables=(
#     "fraction_of_cloud_cover"
# )

template="request_template_era5_single_levels.py"
# Path to the Python script
script="request_code_single_levels.py"

for var in "${variables[@]}"; do
    echo "Requesting $var for $year-$month"

    # Replace placeholders in template
    sed "s/{{VAR}}/$var/g; s/{{MONTH}}/$month/g; s/{{YEAR}}/$year/g" "$template" > "$script"



    # Run the modified Python script
    python3 "$script"

    echo "--------------------------------------"
done

rm "$script"
