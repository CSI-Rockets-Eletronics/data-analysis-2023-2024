# %%
import os
from datetime import timedelta
from ipywidgets import interactive, FloatSlider, Layout
import pandas as pd
import matplotlib.pyplot as plt
import hashlib

# %%
# Dont edit this cell!!!
# Run this cell once to fetch all the data into memory.

# set the server to fetch from
# base_url = "https://csiwiki.me.columbia.edu/rocketsdata2"
base_url = "https://csi-fs-pi-data-server.ngrok.io"
# base_url = "http://fs-pi.local:3000"

cache_dir = "snapshots/far-launch"

start_ts = 1717880525919739  # in unix microseconds
window = 20 * 60 * 1e6  # 20 minutes

end_ts = start_ts + window


def get_csv_with_cache(url):
    url_hash = hashlib.sha256(url.encode()).hexdigest()

    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, url_hash + ".csv")

    if os.path.exists(cache_file):
        return pd.read_csv(cache_file)

    df = pd.read_csv(url)
    df.to_csv(cache_file, index=False)
    return df


def fetch(device, base_url=base_url):
    # limit fetch window to avoid timeouts
    max_window_ts = int(timedelta(minutes=20).total_seconds() * 1e6)

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
lc = fetch("LoadCell", base_url)

# %%
# flip load cell data
lc["thrust"] = -lc["data"]
# compute moving median of load cell data
median_window = 10
lc["thrust_med"] = lc["thrust"].rolling(window=30).median()
# convert mpsi to psi
sci["t1_psi"] = sci["t1"] / 1000


# %%
# create interactive scatter plot with sliders for start and end ts
def plot_thrust(start_ts, end_ts):
    # filter data based on start and end ts
    filtered_lc = lc[(lc["ts"] >= start_ts) & (lc["ts"] <= end_ts)]

    # plot scatter plot of thrust_med with tiny . markers
    plt.scatter(
        filtered_lc["ts"], filtered_lc["thrust"], marker=".", s=1, label="Thrust"
    )
    plt.scatter(
        filtered_lc["ts"],
        filtered_lc["thrust_med"],
        marker=".",
        s=1,
        label="Moving Median",
    )
    plt.xlabel("Timestamp")
    plt.ylabel("Thrust (lbs)")
    plt.title("Moving Median of Load Cell Data")
    plt.legend()
    plt.show()


# create sliders for start and end ts
start_slider = FloatSlider(
    min=start_ts,
    max=end_ts,
    step=1,
    value=start_ts,
    description="Start TS:",
    layout=Layout(width="100%"),
)
end_slider = FloatSlider(
    min=start_ts,
    max=end_ts,
    step=1,
    value=end_ts,
    description="End TS:",
    layout=Layout(width="100%"),
)

# create interactive plot
interactive_plot = interactive(plot_thrust, start_ts=start_slider, end_ts=end_slider)
interactive_plot

# %%


# create interactive scatter plot with sliders for start and end ts
def plot_t1(start_ts, end_ts):
    # filter data based on start and end ts
    filtered_sci = sci[(sci["ts"] >= start_ts) & (sci["ts"] <= end_ts)]

    # plot scatter plot of ox tank transducer data
    plt.scatter(filtered_sci["ts"], filtered_sci["t1_psi"], marker=".", s=1)
    plt.xlabel("Timestamp")
    plt.ylabel("Pressure (psi)")
    plt.title("Interactive Plot of Ox Tank Transducer")
    plt.show()


# create interactive plot
interactive_plot = interactive(plot_t1, start_ts=start_slider, end_ts=end_slider)
interactive_plot

# %%
