"""
Compute catchment-averaged daily potential evapotranspiration (PET) using the
Hargreaves-Samani (1985) formula.

PET = 0.0023 * Ra * (Tmean + 17.8) * (Tmax - Tmin)^0.5      [mm/day]

Ra (extraterrestrial radiation) is calculated from latitude and day-of-year
following the FAO-56 procedure (Allen et al., 1998).

Input : CSV with columns date, Tmean_degC, Tmax_degC, Tmin_degC
Output: CSV with an added PET_mm column
"""

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------
INPUT_CSV = "Kammel_temperature_DWD_1983_2020_partial.csv"
OUTPUT_CSV = "Kammel_PET_Hargreaves_1983_2020.csv"
LATITUDE_DEG = 48.46342  # catchment centroid latitude, °N

GSC = 0.0820  # solar constant, MJ m-2 min-1


def extraterrestrial_radiation(doy: np.ndarray, lat_deg: float) -> np.ndarray:
    """
    FAO-56 extraterrestrial radiation Ra [MJ m-2 day-1] for given day-of-year.
    """
    lat_rad = np.radians(lat_deg)

    # Inverse relative distance Earth-Sun
    dr = 1 + 0.033 * np.cos(2 * np.pi / 365 * doy)

    # Solar declination
    decl = 0.409 * np.sin(2 * np.pi / 365 * doy - 1.39)

    # Sunset hour angle
    ws = np.arccos(np.clip(-np.tan(lat_rad) * np.tan(decl), -1, 1))

    # Extraterrestrial radiation
    ra = (24 * 60 / np.pi) * GSC * dr * (
        ws * np.sin(lat_rad) * np.sin(decl)
        + np.cos(lat_rad) * np.cos(decl) * np.sin(ws)
    )
    return ra  # MJ m-2 day-1


def main():
    df = pd.read_csv(INPUT_CSV, parse_dates=["date"])

    doy = df["date"].dt.dayofyear.to_numpy()
    ra_mj = extraterrestrial_radiation(doy, LATITUDE_DEG)

    # Convert Ra from MJ m-2 day-1 to mm/day equivalent
    ra_mm = ra_mj * 0.408

    df["Ra_MJ_m2_day"] = ra_mj
    df["Ra_mm_day"] = ra_mm

    df["PET_mm"] = (
        0.0023
        * ra_mm
        * (df["Tmean_degC"] + 17.8)
        * np.sqrt(np.clip(df["Tmax_degC"] - df["Tmin_degC"], 0, None))
    )

    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Saved {len(df)} rows to {OUTPUT_CSV}")
    print(df[["date", "Tmean_degC", "Tmax_degC", "Tmin_degC", "PET_mm"]].head())


if __name__ == "__main__":
    main()
