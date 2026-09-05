from pathlib import Path

import numpy as np
import pandas as pd


POSITIONS = Path("outputs/established_dispersal_centrality/event_member_positions_established_segments.parquet")
EVENTS = Path("outputs/established_dispersal_centrality/established_events.csv")
OUT = Path("outputs/established_dispersal_radius_samples")
RADII = (1, 2, 5)


def haversine_m(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, (lon1, lat1, lon2, lat2))
    dlon, dlat = lon2 - lon1, lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * 6_371_000 * np.arcsin(np.sqrt(a))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pos = pd.read_parquet(POSITIONS)
    events = pd.read_csv(EVENTS, usecols=["event_id", "established_outcome"])
    rows = []
    for (event_id, bin_2min), z in pos.groupby(["event_id", "bin_2min"], observed=True, sort=False):
        focal_id = str(z.animal_id_focal.iloc[0])
        focal = z[z.animal_id.astype(str).eq(focal_id)]
        others = z[~z.animal_id.astype(str).eq(focal_id)]
        if focal.empty or others.empty:
            continue
        f = focal.iloc[0]
        distances = haversine_m(f.lon, f.lat, others.lon.to_numpy(), others.lat.to_numpy())
        nearest = float(np.nanmin(distances))
        rows.append((event_id, focal_id, bin_2min, int(z.week.iloc[0]), nearest))
    bins = pd.DataFrame(rows, columns=["event_id", "animal_id", "bin_2min", "week", "nearest_recipient_m"])
    bins = bins.merge(events, on="event_id", how="left")
    for radius in RADII:
        bins[f"contact_{radius}m"] = bins.nearest_recipient_m.le(radius)
    bins.to_parquet(OUT / "recipient_contact_2min_bins_1m_2m_5m.parquet", index=False)

    agg = {"event_id": "nunique", "bin_2min": "size", "week": "nunique"}
    for radius in RADII:
        agg[f"contact_{radius}m"] = "sum"
    summary = (bins.groupby(["animal_id", "established_outcome"], observed=True)
               .agg(**{
                   "recipient_segments": ("event_id", "nunique"),
                   "two_minute_values": ("bin_2min", "size"),
                   "event_week_values": ("week", "size"),
                   **{f"contact_values_{r}m": (f"contact_{r}m", "sum") for r in RADII},
               }).reset_index())
    # Count unique event-week estimates, not merely distinct week numbers.
    ew = (bins[["animal_id", "established_outcome", "event_id", "week"]].drop_duplicates()
          .groupby(["animal_id", "established_outcome"], observed=True).size().rename("event_week_values").reset_index())
    summary = summary.drop(columns="event_week_values").merge(ew, on=["animal_id", "established_outcome"])
    for radius in RADII:
        summary[f"proportion_within_{radius}m"] = summary[f"contact_values_{radius}m"] / summary.two_minute_values
    summary.to_csv(OUT / "values_per_individual_1m_2m_5m.csv", index=False)
    print(summary.to_string(index=False))
    print(f"\nIndividuals={summary.animal_id.nunique()}, individual-outcome rows={len(summary)}, bins={len(bins):,}")


if __name__ == "__main__":
    main()
