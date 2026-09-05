from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


SOURCE = Path("outputs/established_dispersal_radius_samples/values_per_individual_1m_2m_5m.csv")
OUT = Path("outputs/established_dispersal_radius_samples/high_support_radius_integration.png")


def main():
    data = pd.read_csv(SOURCE)
    data = (data.groupby("animal_id", as_index=False)
            .agg(two_minute_values=("two_minute_values", "sum"),
                 contact_values_1m=("contact_values_1m", "sum"),
                 contact_values_2m=("contact_values_2m", "sum"),
                 contact_values_5m=("contact_values_5m", "sum")))
    data = data[data.contact_values_1m > 50].copy()
    radii = [1, 2, 5]
    fig, ax = plt.subplots(figsize=(9, 6))
    for row in data.itertuples(index=False):
        values = [getattr(row, f"contact_values_{r}m") / row.two_minute_values for r in radii]
        ax.plot(radii, values, marker="o", linewidth=2, label=f"{row.animal_id} (n={row.two_minute_values:,})")
    ax.set_xticks(radii, ["1 m", "2 m", "5 m"])
    ax.set_xlabel("Distance threshold")
    ax.set_ylabel("Proportion of eligible 2-minute bins")
    ax.set_title("Recipient-group integration: individuals with >50 contacts at 1 m")
    ax.grid(alpha=.25)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT, dpi=220)
    plt.close(fig)
    print(data[["animal_id", "two_minute_values", "contact_values_1m"]].to_string(index=False))
    print(OUT.resolve())


if __name__ == "__main__":
    main()
