#!/bin/bash


# Check for month and year arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <month> <year>"
    exit 1
fi

month="$1"
year="$2"
mkdir -p "/pesq/dados/monan/users/madeleine.gacita/global_data/era5_land/$year/$month"


#List of meteorological variables
variables=(
    "2m_dewpoint_temperature"
    "2m_temperature"
    "skin_temperature"
    "forecast_albedo"
    "surface_latent_heat_flux"
    "surface_net_solar_radiation"
    "surface_net_thermal_radiation"
    "surface_sensible_heat_flux"
    "surface_solar_radiation_downwards"
    "surface_thermal_radiation_downwards"
    "10m_u_component_of_wind"
    "10m_v_component_of_wind"
    "surface_pressure"
    "total_precipitation"
    "leaf_area_index_high_vegetation"
    "leaf_area_index_low_vegetation"
    "high_vegetation_cover"
    "glacier_mask"
    "lake_cover"
    "low_vegetation_cover"
    "lake_total_depth"
    "land_sea_mask"
    "soil_type"
    "type_of_high_vegetation"
    "type_of_low_vegetation"
)
# variables=(
#     "fraction_of_cloud_cover"
# )

template="request template_land.py"
# Path to the Python script
script="request_code_land.py"

for var in "${variables[@]}"; do
    echo "Requesting $var for $year-$month"

    # Replace placeholders in template
    sed "s/{{VAR}}/$var/g; s/{{MONTH}}/$month/g; s/{{YEAR}}/$year/g" "$template" > "$script"



    # Run the modified Python script
    python3 "$script"

    echo "--------------------------------------"
done

rm "$script"
