import cdsapi
import calendar

# Input parameters
month = "{{MONTH}}"
year = "{{YEAR}}"


# Calculate number of days in the month
num_days = calendar.monthrange(int(year), int(month))[1]
days = [f"{day:02d}" for day in range(1, num_days + 1)]

dataset = "reanalysis-era5-pressure-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": ["{{VAR}}"], # This line will be modified
    "year": ["{{YEAR}}"],
    "month": ["{{MONTH}}"],
    "day": days,
    "time": [
        "00:00", "03:00", "06:00", 
        "09:00", "12:00", "15:00", 
        "18:00", "21:00"
    ],
    "pressure_level": [
        "3", "10", "20",
        "30", "50", "70", 
        "100", "150", "200",
        "250", "300", "400", 
        "500",  "700", "775", 
        "850", "925", "1000"
    ],
    "data_format": "grib",
    "download_format": "unarchived"
}

client = cdsapi.Client()
#client.retrieve(dataset, request).download()

# Download the data directly to the target file
target_file = f"/pesq/dados/monan/users/madeleine.gacita/global_data/era5/pressure_levels/{{YEAR}}/{{MONTH}}/{{VAR}}.grib"
client.retrieve(dataset, request, target_file)
