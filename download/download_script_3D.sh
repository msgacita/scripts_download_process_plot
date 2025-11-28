#!/bin/bash


# Check for month and year arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <month> <year>"
    exit 1
fi

month="$1"
year="$2"


# Create directory for the month if it doesn't exist
mkdir -p "/pesq/dados/monan/users/madeleine.gacita/global_data/era5/pressure_levels/$year/$month"
#List of meteorological variables
variables=(
    "divergence"
    "fraction_of_cloud_cover"
    "geopotential"
    "ozone_mass_mixing_ratio"
    "potential_vorticity"
    "relative_humidity"
    "specific_cloud_ice_water_content"
    "specific_cloud_liquid_water_content"
    "specific_humidity"
    "specific_rain_water_content"
    "specific_snow_water_content"
    "temperature"
    "u_component_of_wind"
    "v_component_of_wind"
    "vertical_velocity"
    "vorticity"
)
# variables=(
#     "fraction_of_cloud_cover"
# )

template="request_template_era5_pressure_levels.py"
# Path to the Python script
script="request_code_pressure_levels.py"

for var in "${variables[@]}"; do
   echo "Requesting $var for $year-$month"

    # Replace placeholders in template
    sed "s/{{VAR}}/$var/g; s/{{MONTH}}/$month/g; s/{{YEAR}}/$year/g" "$template" > "$script"


    # Run the modified Python script
    python3 "$script"

    echo "--------------------------------------"
done


