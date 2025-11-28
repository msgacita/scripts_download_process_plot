import cdsapi
import calendar

# Input parameters
month = "{{MONTH}}"
year = "{{YEAR}}"


# Calculate number of days in the month
num_days = calendar.monthrange(int(year), int(month))[1]
days = [f"{day:02d}" for day in range(1, num_days + 1)]


dataset = "reanalysis-era5-land"
request = {
    "variable": ["{{VAR}}"], # This line will be modified
    "year": ["{{YEAR}}"],
    "month": ["{{MONTH}}"],
    "day": days,
    "time": [
        "00:00", "01:00", "02:00",
        "03:00", "04:00", "05:00",
        "06:00", "07:00", "08:00",
        "09:00", "10:00", "11:00",
        "12:00", "13:00", "14:00",
        "15:00", "16:00", "17:00",
        "18:00", "19:00", "20:00",
        "21:00", "22:00", "23:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived"
}

client = cdsapi.Client()
#client.retrieve(dataset, request).download()

# Download the data directly to the target file
target_file = f"/pesq/dados/monan/users/madeleine.gacita/global_data/era5/land/{{YEAR}}/{{MONTH}}/{{VAR}}.grib"
client.retrieve(dataset, request, target_file)