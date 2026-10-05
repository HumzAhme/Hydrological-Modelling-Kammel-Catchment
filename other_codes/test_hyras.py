import os
import pandas as pd
import xarray as xr
import rioxarray
import geopandas as gpd

# ── Paths ─────────────────────────────────────────────────────────────────────
dest_dir = r"C:\Uni\Hydrological Modeling\kammel_project"

dir_tmean = os.path.join(dest_dir, "mean_temp")
dir_tmax  = os.path.join(dest_dir, "max_temp")
dir_tmin  = os.path.join(dest_dir, "min_temp")

catchment_path = os.path.join(dest_dir, "catchment.gpkg")

out_csv = os.path.join(dest_dir, "Kammel_temperature_DWD_1983_2020_partial.csv")
out_pkl = os.path.join(dest_dir, "Kammel_temperature_DWD_1983_2020_partial.pkl")

years = range(1983, 2021)

laea_crs = (
    "+proj=laea +lat_0=52 +lon_0=10 "
    "+x_0=4321000 +y_0=3210000 "
    "+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 "
    "+units=m +no_defs +type=crs"
)

# ── Load catchment ────────────────────────────────────────────────────────────
kammel = gpd.read_file(catchment_path)
kammel = kammel.to_crs(laea_crs)


def open_dataset_with_debug(filepath):
    print("  exists:", os.path.exists(filepath))

    if not os.path.exists(filepath):
        raise FileNotFoundError(filepath)

    print("  size MB:", round(os.path.getsize(filepath) / 1024 / 1024, 2))

    try:
        ds = xr.open_dataset(filepath, engine="netcdf4")
        print("  opened with netcdf4")
        return ds
    except Exception as e:
        print("  netcdf4 failed:")
        print(" ", repr(e))

    try:
        ds = xr.open_dataset(filepath, engine="h5netcdf")
        print("  opened with h5netcdf")
        return ds
    except Exception as e:
        print("  h5netcdf failed:")
        print(" ", repr(e))
        raise


def process_hyras_temp(filepath, varname):
    print(f"Processing {os.path.basename(filepath)}")

    ds = open_dataset_with_debug(filepath)

    da = ds[varname]

    da = da.rio.write_crs(laea_crs)
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")

    da_clip = da.rio.clip(kammel.geometry, kammel.crs, drop=True)

    ts = da_clip.mean(dim=("x", "y"), skipna=True)

    df = ts.to_dataframe(name=varname).reset_index()
    df = df[["time", varname]]
    df["time"] = pd.to_datetime(df["time"]).dt.date

    ds.close()
    return df


# ── Load existing progress if available ───────────────────────────────────────
all_years = []
processed_years = set()

if os.path.exists(out_csv):
    print("Existing partial file found. Loading progress...")
    old_df = pd.read_csv(out_csv)
    old_df["date"] = pd.to_datetime(old_df["date"])

    old_df = old_df.drop_duplicates(subset="date")
    old_df = old_df.sort_values("date")

    processed_years = set(old_df["date"].dt.year.unique())
    all_years.append(old_df)

    print("Already processed years:", sorted(processed_years))
    print("Last saved date:", old_df["date"].max().date())
else:
    print("No partial file found. Starting from scratch.")


# ── Main loop ─────────────────────────────────────────────────────────────────
for yr in years:

    if yr in processed_years:
        print(f"Skipping {yr}, already processed.")
        continue

    print(f"\n=== Year {yr} ===")

    try:
        f_tmean = os.path.join(dir_tmean, f"tas_hyras_1_{yr}_v6-0_de.nc")
        f_tmax  = os.path.join(dir_tmax,  f"tasmax_hyras_1_{yr}_v6-0_de.nc")
        f_tmin  = os.path.join(dir_tmin,  f"tasmin_hyras_1_{yr}_v6-0_de.nc")

        tmean = process_hyras_temp(f_tmean, "tas")
        tmax  = process_hyras_temp(f_tmax, "tasmax")
        tmin  = process_hyras_temp(f_tmin, "tasmin")

        df_year = (
            tmean
            .merge(tmax, on="time")
            .merge(tmin, on="time")
            .rename(columns={
                "time": "date",
                "tas": "Tmean_degC",
                "tasmax": "Tmax_degC",
                "tasmin": "Tmin_degC"
            })
        )

        df_year["date"] = pd.to_datetime(df_year["date"])

        all_years.append(df_year)

        temp_df = pd.concat(all_years, ignore_index=True)
        temp_df = temp_df.drop_duplicates(subset="date")
        temp_df = temp_df.sort_values("date")

        temp_df.to_csv(out_csv, index=False)
        temp_df.to_pickle(out_pkl)

        print(f"Saved progress up to {yr}")
        print("Current last date:", temp_df["date"].max().date())

    except Exception as e:
        print(f"\nERROR in year {yr}")
        print(repr(e))
        print("\nProgress before this year is saved.")
        break