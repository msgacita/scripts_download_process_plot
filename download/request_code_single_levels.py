import cdsapi
import calendar

# Input parameters
month = "11"
year = "2018"

dataset = "reanalysis-era5-single-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis_by_hour_of_day"],
    "variable": ["2m_dewpoint_temperature"],
    "year": ["2018"],
    "month": ["11"],
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
target_file = f"/pesq/dados/monan/users/madeleine.gacita/global_data/era5/single_levels/2018/11/2m_dewpoint_temperature_MHour.nc"
client.retrieve(dataset, request, target_file)