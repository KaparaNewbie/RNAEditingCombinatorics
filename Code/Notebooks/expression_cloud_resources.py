# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import re
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np
import plotly.io as pio
import plotly.express as px
from icecream import ic

pio.templates.default = "plotly_white"


# %%
def parse_julia_log_file(
    log_file: Path | str,
    run_name: str
) -> pd.DataFrame:

    # keys we care about inside each block
    FIELDS = ["samplename", "distinctfile", "allprotsfile", "readsfile"]
    
    header_re = re.compile(r"┌ Info: (\d{2}\.\d{2}\.\d{4}) - (\d{2}:\d{2}:\d{2})\s+([^\n]+)")
    kv_re = re.compile(r"│\s+(\w+)\s*=\s*(.+)")

    rows = []
    current = None

    log_file = Path(log_file) # in case log file is a str

    for line in log_file.read_text().splitlines():
        m = header_re.match(line)
        if m:
            if current:
                rows.append(current)
            date, time, stage = m.groups()
            current = {
                "timestamp": datetime.strptime(f"{date} {time}", "%d.%m.%Y %H:%M:%S"),
                "stage": stage.strip(),
            }
            continue

        if current:
            m = kv_re.match(line)
            if m:
                k, v = m.groups()
                if k in FIELDS:
                    current[k] = v.strip().strip('"')

    # last block
    if current:
        rows.append(current)

    df = pd.DataFrame(rows).sort_values("timestamp").reset_index(drop=True)
    
    df.insert(0, "runname", run_name)
    
    return df


# %%
runs_names = [
    "RUSC2_40", 
    "RUSC2_60",
    "RUSC2_80",
    "TWK7_60",
    "ANR17_60",
    "ANR17_80"
]

pidstat_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/RUSC2_40/julia.RUSC2..pidstat",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/RUSC2_60/julia.RUSC2..pidstat",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/RUSC2_80/julia.RUSC2..pidstat",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/TWK7_60/julia.TWK7..pidstat",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/ANR17_60/julia.ANR17..pidstat",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/ANR17_80/julia.ANR17_80T..pidstat"
]

julia_log_files = [
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/RUSC2_40/expressionlevels.EntropyConsidered.28.01.2026.log",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/RUSC2_60/expressionlevels.EntropyConsidered.28.01.2026.log",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/RUSC2_80/expressionlevels.EntropyConsidered.28.01.2026.log",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/TWK7_60/expressionlevels.EntropyConsidered.06.02.2026.log",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/ANR17_60/expressionlevels.EntropyConsidered.05.02.2026.log",
    "/private6/projects/Combinatorics/D.pealeii/MpileupAndTranscripts/Illumina/TestFixedExpression/CloudLogs/ANR17_80/expressionlevels.EntropyConsidered.09.02.2026.log"
]

# %%
main_stages = [
    "run_sample",
    "prepare_allprotsdf!",
    "prepare_allprotsdf!",
    "prepare_readsdf",
    "prepare_distinctdf",
    "prepare_allprotsdf! - finished",
    "distances",
    "additional_assignments",
    "distances - main function",
    "finddistinctAAsets",
    "choosesolutions",
    "prepare_solution_data_for_reassignment",
    "threaded_one_solution_additional_assignment_considering_available_reads",
    "additional_assignments - solution start",
    "threaded_one_solution_additional_assignment_considering_available_reads - assignment plan",
    "merge_chosendf_results!",
    "finalize_solution_data",
    "finalize_solution_data",
    "main - finished sucessfully!"
]

# %%
julia_log_dfs = [
    parse_julia_log_file(log_file, run_name)
    for run_name, log_file in zip(runs_names, julia_log_files)
]
concat_julia_log_df = pd.concat(julia_log_dfs, ignore_index=True)
concat_julia_log_df = concat_julia_log_df.loc[
    concat_julia_log_df["stage"].isin(main_stages),
    ["runname", "timestamp", "stage"]
]
concat_julia_log_df

# %%
starting_dates_df = concat_julia_log_df.groupby("runname")["timestamp"].min().dt.date
starting_dates_df

# %%
# df = pd.read_csv(pidstat_file, sep='\s+', skiprows=2)

# # the titles show every time a sample is taken -
# # we want to remove these rows
# df = df.loc[
#     ~df["#"].apply(lambda x: x.startswith("#"))
# ].reset_index(drop=True)


# # the 1st column (#) is time in 12-hour format, 
# # the 2nd column (Time) is AM/PM indicator
# # we want to combine these two columns into a single datetime column
# df = df.rename(
#     columns={
#         "#": "Time_12",
#         "Time": "AM/PM"
#     }
# )

# df.insert(
#     0,
#     "Time",
#     df["Time_12"] + " " + df["AM/PM"]
# )

# df = df.drop(
#     columns=["Time_12", "AM/PM"]
# )

# # now we can change the Time column to datetime
# df["Time"] = pd.to_datetime(df["Time"], format="%I:%M:%S %p").dt.time

# df["%CPU"] = pd.to_numeric(df["%CPU"], errors='coerce')
# df["%MEM"] = pd.to_numeric(df["%MEM"], errors='coerce')
# df["RSS"] = pd.to_numeric(df["RSS"], errors='coerce')

# df.insert(
#     df.columns.get_loc("RSS") + 1,
#     "RSS_MB",
#     df["RSS"] / 1024
# )

# df_tg = df.loc[
#     df["TGID"].ne("-")
# ].reset_index(drop=True)
# df_tg

# df

# %%
def parse_pidstat_file(
    pidstat_file: str | Path,
    run_name: str,
    start_date: datetime.date
) -> pd.DataFrame:

    df = pd.read_csv(pidstat_file, sep='\s+', skiprows=2)

    # the titles show every time a sample is taken -
    # we want to remove these rows
    df = df.loc[
        ~df["#"].apply(lambda x: x.startswith("#"))
    ].reset_index(drop=True)


    df = df.rename(
        columns={
            "Time": "TempBadCol1"
        }
    )
    df = df.rename(
        columns={
            "#": "Time"
        }
    )

    for i in range(1, len(df.columns[1:])):
        j = i + 1
        col_i = df.columns[i]
        col_j = df.columns[j]
        # ic(i, col_i, col_j)
        new_col_j = f"TempBadCol{j}"
        df = df.rename(
            columns={
                col_j: new_col_j
            }
        )
        df = df.rename(
            columns={
                col_i: col_j
            }
        )
        
    last_bad_col = df.columns[-1]
    df = df.drop(columns=last_bad_col)

    df["%CPU"] = pd.to_numeric(df["%CPU"], errors='coerce')
    df["%MEM"] = pd.to_numeric(df["%MEM"], errors='coerce')
    df["RSS"] = pd.to_numeric(df["RSS"], errors='coerce')

    df.insert(
        df.columns.get_loc("RSS") + 1,
        "RSS_MB",
        df["RSS"] / 1024
    )
    
    df.insert(0, "runname", run_name)
    
    # Create a real timestamp from pidstat_df["Time"] using start_date + rollover past midnight
    # Assumes pidstat_df is already filtered to a single run (it is in your notebook).

    times = pd.to_datetime(df["Time"], format="%H:%M:%S", errors="coerce")
    if times.isna().any():
        bad = df.loc[times.isna(), "Time"].head()
        raise ValueError(f"Unparseable Time values (showing up to 5):\n{bad}")

    # detect day rollovers (time going "backwards")
    secs = times.dt.hour * 3600 + times.dt.minute * 60 + times.dt.second
    day_offset = (secs.diff().fillna(0) < 0).cumsum()

    df["timestamp"] = (
        pd.to_datetime(start_date)  # midnight of start_date
        + pd.to_timedelta(day_offset, unit="D")
        + pd.to_timedelta(times.dt.strftime("%H:%M:%S"))
    )
    
    df["Time"] = df["timestamp"]
    del df["timestamp"]

    return df

# %%
pidstat_dfs = [
    parse_pidstat_file(pidstat_file, run_name, starting_dates_df.loc[run_name])
    for run_name, pidstat_file in zip(runs_names, pidstat_files)
]
concat_pidstat_df = pd.concat(pidstat_dfs, ignore_index=True)
concat_pidstat_df

# %%
for run_name in runs_names:
    avg_cpu = np.round(
        concat_pidstat_df.loc[
            concat_pidstat_df["runname"] == run_name,
            "%CPU"
        ].mean() / 100,
        2
    )
    ic(run_name, avg_cpu)

# %%
concat_julia_log_df

# %%
run_name = runs_names[-1]
julia_log_df = concat_julia_log_df.loc[
        concat_julia_log_df["runname"] == run_name
    ].sort_values("timestamp").reset_index(drop=True)
julia_log_df

# %%
main_main_stages = [
    "distances",
    "additional_assignments",
]

# %%
for run_name in runs_names:
    fig = px.area(
        concat_pidstat_df.loc[
            concat_pidstat_df["runname"] == run_name
        ], 
        x="Time", 
        y="%CPU",
    )
    julia_log_df = concat_julia_log_df.loc[
        concat_julia_log_df["runname"] == run_name
    ]
    # for stage in main_stages:
    for stage in main_main_stages:
        if stage in julia_log_df["stage"].values:
            timestamp = (
                julia_log_df.loc[julia_log_df["stage"] == stage, "timestamp"]
                .iloc[0]
                .to_pydatetime()
            )
            fig.add_vline(
                x=timestamp,
                line_width=4,
                line_color="red",
            )            
    fig.update_layout(
        title="CPU usage over time - " + run_name,
        xaxis_title="Time",
        yaxis_title="% CPU",
        width=1200,
        # height=500
        height=280
    )
    fig.show()

# %%
for run_name in runs_names:
    fig = px.area(
        concat_pidstat_df.loc[
            concat_pidstat_df["runname"] == run_name
        ], 
        x="Time", 
        y="%MEM",
    )
    julia_log_df = concat_julia_log_df.loc[
        concat_julia_log_df["runname"] == run_name
    ]
    # for stage in main_stages:
    for stage in main_main_stages:
        if stage in julia_log_df["stage"].values:
            timestamp = (
                julia_log_df.loc[julia_log_df["stage"] == stage, "timestamp"]
                .iloc[0]
                .to_pydatetime()
            )
            fig.add_vline(
                x=timestamp,
                line_width=4,
                line_color="red",
            )            
    fig.update_layout(
        title="Memory usage over time - " + run_name,
        xaxis_title="Time",
        yaxis_title="% Memory",
        width=1200,
        # height=500
        height=280
    )
    fig.show()

# %%
for run_name in runs_names:
    fig = px.area(
        concat_pidstat_df.loc[
            concat_pidstat_df["runname"] == run_name
        ], 
        x="Time", 
        y="RSS_MB",
    )
    julia_log_df = concat_julia_log_df.loc[
        concat_julia_log_df["runname"] == run_name
    ]
    # for stage in main_stages:
    for stage in main_main_stages:
        if stage in julia_log_df["stage"].values:
            timestamp = (
                julia_log_df.loc[julia_log_df["stage"] == stage, "timestamp"]
                .iloc[0]
                .to_pydatetime()
            )
            fig.add_vline(
                x=timestamp,
                line_width=4,
                line_color="red",
            )            
    fig.update_layout(
        title="Memory usage over time - " + run_name,
        xaxis_title="Time",
        yaxis_title="RSS (MB)",
        width=1200,
        # height=500
        height=400
    )
    fig.show()
