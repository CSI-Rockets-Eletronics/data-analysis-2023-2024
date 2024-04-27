# %%
import os
from datetime import datetime, timedelta
from ipywidgets import interactive, FloatSlider, Layout
from pytz import timezone
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import hashlib

# %%
# Dont edit this cell!!!
# Run this cell once to fetch all the data into memory.

# set the server to fetch from
# base_url = "https://csiwiki.me.columbia.edu/rocketsdata2"
base_url = "http://fs-pi.local:3000"

# set time range for fetching the full data
start = datetime(year=2024, month=4, day=26, hour=16, minute=25)
window = timedelta(minutes=10)


def get_csv_with_cache(url):
    url_hash = hashlib.sha256(url.encode()).hexdigest()

    os.makedirs("cache", exist_ok=True)
    cache_file = os.path.join("cache", url_hash + ".csv")

    if os.path.exists(cache_file):
        return pd.read_csv(cache_file)

    df = pd.read_csv(url)
    df.to_csv(cache_file, index=False)
    return df


def to_ts_range(start, window):
    start_ts = int(timezone("EST").localize(start).timestamp() * 1e6)
    end_ts = start_ts + int(window.total_seconds() * 1e6)
    return start_ts, end_ts


def fetch(device, base_url=base_url):
    # limit fetch window to avoid timeouts
    max_window_ts = int(timedelta(minutes=20).total_seconds() * 1e6)

    start_ts, end_ts = to_ts_range(start, window)
    cur_start_ts = start_ts
    df = None

    while cur_start_ts <= end_ts:
        cur_end_ts = min(cur_start_ts + max_window_ts, end_ts)
        url = f"{base_url}/export/0/latest/{device}/records?&startTs={cur_start_ts}&endTs={cur_end_ts}"
        cur_df = get_csv_with_cache(url)
        df = cur_df if df is None else pd.concat([df, cur_df])
        cur_start_ts += max_window_ts + 1

    df = df.sort_values(by=["ts"])
    df = df.reset_index(drop=True)
    return df


# fetch full data
fs = fetch("FiringStation", base_url)
sci = fetch("Scientific", base_url)
roc = fetch("RocketScientific", base_url)
lc1 = fetch("LoadCell1", base_url)
lc2 = fetch("LoadCell2", base_url)

# %%
# flip load cell data
lc1["thrust"] = -lc1["data"]
lc2["thrust"] = -lc2["data"]

# compute moving median of load cell data
median_window = 10
lc1["thrust_med"] = lc1["thrust"].rolling(window=30).median()
lc2["thrust_med"] = lc2["thrust"].rolling(window=30).median()

# calculate prefix sum of total impulse to speed up window calculations
lc1["thrust_int"] = (lc1["thrust"] * lc1["ts"].diff()).cumsum()
lc2["thrust_int"] = (lc2["thrust"] * lc2["ts"].diff()).cumsum()

# sum load cell data
lc_sum = pd.merge_asof(lc1, lc2, on="ts", suffixes=("_lc1", "_lc2"))
lc_sum["thrust"] = lc_sum["thrust_lc1"] + lc_sum["thrust_lc2"]
lc_sum["thrust_med"] = lc_sum["thrust"].rolling(window=30).median()

# convert mpsi to psi
sci["st1_psi"] = sci["st1"] / 1000
sci["st2_psi"] = sci["st2"] / 1000
roc["bt1_psi"] = roc["bt1"] / 1000
roc["bt2_psi"] = roc["bt2"] / 1000

# clean CC transducer data
roc["bt1_psi_med"] = roc["bt1_psi"].rolling(window=30).median()
roc["bt2_psi_med"] = roc["bt2_psi"].rolling(window=30).median()

# clean thermocouple data
fs["thermo1C_med"] = fs["thermo1C"].rolling(window=10).median()
fs["thermo2C_med"] = fs["thermo2C"].rolling(window=10).median()

# create new ts column that's the relative ts.1, but shifted to absolute ts
sci["ts_precise"] = sci["ts.1"] + int(sci["ts"].mean() - sci["ts.1"].mean())
roc["ts_precise"] = roc["ts.1"] + int(roc["ts"].mean() - roc["ts.1"].mean())


# %%
# set this to be larger than the frequencies of the individuals dataframes,
# but not to large to conserve compute
COMMON_HZ = 1000

# created merged df with dense, uniformly spaced timestamps
merged = pd.DataFrame({"ts": range(*to_ts_range(start, window), int(1e6 / COMMON_HZ))})

merged = pd.merge_asof(merged, fs[["ts", "thermo1C", "thermo2C"]], on="ts")
merged = pd.merge_asof(
    merged,
    sci[["ts_precise", "st1_psi", "st2_psi"]].rename(columns={"ts_precise": "ts"}),
    on="ts",
)
merged = pd.merge_asof(
    merged,
    roc[["ts_precise", "bt1_psi", "bt2_psi"]].rename(columns={"ts_precise": "ts"}),
    on="ts",
)
merged = pd.merge_asof(
    merged, lc1[["ts", "thrust"]].rename(columns={"thrust": "thrust_lc1"}), on="ts"
)
merged = pd.merge_asof(
    merged, lc2[["ts", "thrust"]].rename(columns={"thrust": "thrust_lc2"}), on="ts"
)

# %%
state_names = {
    0: "fire",
    1: "fill",
    2: "purge",
    3: "abort",
    4: "standby",
    5: "keep",
    6: "pulseFillA",
    7: "pulseFillB",
    8: "pulseFillC",
    25: "pulseVentA",
    26: "pulseVentB",
    27: "pulseVentC",
    30: "pulsePurgeA",
    31: "pulsePurgeB",
    32: "pulsePurgeC",
    20: "fireManualIgniter",
    21: "fireManualValve",
    40: "custom",
}


# %%
def create_windows(df_plot, fs, xlim):
    """
    Function to create the rectangles based on Firing Station data given a plot
    to add them to.
    """

    cmap = matplotlib.colormaps["Spectral"]

    ylim1 = df_plot.get_ylim()
    # print(ylim1)
    rects = []
    last_ts = fs["ts"][0]  # not zero

    for i in range(len(fs)):
        if fs["ts"][i] < xlim[0]:
            last_ts = fs["ts"][i]
            continue
        if i == 0:
            continue
        if fs["ts"][i] > xlim[1]:
            rects.append(
                matplotlib.patches.Rectangle(
                    (last_ts, ylim1[0]),
                    fs["ts"][i - 1] - last_ts,
                    ylim1[1] - ylim1[0],
                    color=cmap(fs["stateByte"][i - 1] / 40)[:3] + (0.2,),
                )
            )
            df_plot.text(
                last_ts,
                ylim1[0] + (ylim1[1] - ylim1[0]) / 25,
                state_names[fs["stateByte"][i - 1]],
                fontsize=7,
                rotation="vertical",
            )
            # df_plot.text(last_ts, ylim1[0], state_names[fs_sorted["stateByte"][i-1]], fontsize=7, rotation="vertical")
            df_plot.plot(last_ts, ylim1[0], "ro")
            return rects
        if fs["stateByte"][i] == fs["stateByte"][i - 1]:
            continue
        rects.append(
            matplotlib.patches.Rectangle(
                (last_ts, ylim1[0]),
                fs["ts"][i - 1] - last_ts,
                ylim1[1] - ylim1[0],
                color=cmap(fs["stateByte"][i - 1] / 40)[:3] + (0.2,),
            )
        )
        df_plot.text(
            last_ts,
            ylim1[0] + (ylim1[1] - ylim1[0]) / 25,
            state_names[fs["stateByte"][i - 1]],
            fontsize=7,
            rotation="vertical",
        )
        # df_plot.text(last_ts, ylim1[0], state_names[fs_sorted["stateByte"][i-1]], fontsize=7, rotation="vertical")

        # print(last_ts)
        df_plot.plot(last_ts, ylim1[0], "ro")
        # df_plot.text(xlim[0], 5, "hi")
        last_ts = fs["ts"][i]
    return rects


# %%
SHOW_WINDOWS = False

# plot data
full_xlim = to_ts_range(start, window)


def slider(value, min=full_xlim[0] / 1e6, max=full_xlim[1] / 1e6, step=1):
    return FloatSlider(
        value,
        min=min,
        max=max,
        step=step,
        continuous_update=False,
        layout=Layout(width="100%"),
    )


def plot_with_windows(df, y, title, xlim):
    # trim df so ylim is only for data in range
    df = df[(df["ts"] >= xlim[0]) & (df["ts"] <= xlim[1])]
    x_col = "ts.1" if "ts.1" in df else "ts"
    sci_plot = df.plot(x_col, y, title=title, figsize=(15, 3), style=".")

    if SHOW_WINDOWS:
        # https://www.geeksforgeeks.org/matplotlib-axes-axes-add_patch-in-python/
        for rect in create_windows(sci_plot, fs, xlim):
            sci_plot.add_patch(rect)

    plt.show()


# %%
def test(x_min, x_max):
    lc1_range = lc1[(lc1["ts"] >= x_min) & (lc1["ts"] <= x_max)]
    lc2_range = lc2[(lc2["ts"] >= x_min) & (lc2["ts"] <= x_max)]

    lc1_impulse = lc1_range["thrust_int"].iloc[-1] - lc1_range["thrust_int"].iloc[0]
    lc2_impulse = lc2_range["thrust_int"].iloc[-1] - lc2_range["thrust_int"].iloc[0]

    total_impulse = (lc1_impulse + lc2_impulse) / 1e6

    print("Total impulse (pound*sec): ", total_impulse)

    total_impulse = (lc1_impulse * 2) / 1e6

    print("Total impulse using 2*(load cell 1) (pound*sec): ", total_impulse)


def calculate_mass_flow(df, x_min, x_max):
    """
    Integrates thrust using trapezoidal. Returns thrust in pound*sec.
    """

    df_range = df[(df["ts"] >= x_min) & (df["ts"] <= x_max)]
    mass_flow = 0

    for i in range(1, len(df_range)):
        dt = df_range["ts"].iloc[i] - df_range["ts"].iloc[i - 1]
        dt /= 1e6  # convert to seconds
        avg_thrust = (
            (df_range["thrust"].iloc[i] + df_range["thrust"].iloc[i - 1]) / 2
        ) * 4.4482
        mass_flow += avg_thrust / 9.80665 / dt

    return mass_flow


def test2(x_min, x_max):
    lc1_range = lc1[(lc1["ts"] >= x_min) & (lc1["ts"] <= x_max)]
    lc2_range = lc2[(lc2["ts"] >= x_min) & (lc2["ts"] <= x_max)]

    lc1_mdot = calculate_mass_flow(lc1_range, x_min, x_max)
    lc2_mdot = calculate_mass_flow(lc2_range, x_min, x_max)
    total_liquid_mass_flow = lc1_mdot + lc2_mdot

    print("Total liquid mass flow (kg/sec): ", total_liquid_mass_flow)


def plot_fn(x_min_sec, x_max_sec, x_min_fine_sec, x_max_fine_sec):
    x_min = x_min_sec * 1e6 + x_min_fine_sec * 1e6
    x_max = x_max_sec * 1e6 + x_max_fine_sec * 1e6

    print("Min ts (sec): ", x_min / 1e6)
    print("Max ts (sec): ", x_max / 1e6)

    test(x_min, x_max)
    # test2(x_min, x_max)

    plot_with_windows(sci, ["st1_psi", "st2_psi"], "OX Transducers", [x_min, x_max])
    plot_with_windows(roc, ["bt1_psi", "bt2_psi"], "CC Transducers", [x_min, x_max])
    plot_with_windows(lc1, "thrust_med", "Load Cell 1", [x_min, x_max])
    plot_with_windows(lc2, "thrust_med", "Load Cell 2", [x_min, x_max])
    plot_with_windows(lc_sum, "thrust_med", "Load Cell Sum", [x_min, x_max])
    plot_with_windows(fs, ["thermo1C", "thermo2C"], "Thermocouples", [x_min, x_max])


interactive(
    plot_fn,
    x_min_sec=slider(0),
    x_max_sec=slider(1e12),
    x_min_fine_sec=slider(0, min=-1, max=1, step=0.01),
    x_max_fine_sec=slider(0, min=-1, max=1, step=0.01),
)
