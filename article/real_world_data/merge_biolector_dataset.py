import pandas as pd
import pathlib

PROJECT_DIR = pathlib.Path(__file__).parents[2]

# load online data
online_dat_filepath = (
    PROJECT_DIR
    / "article"
    / "real_world_data"
    / "COU_S2_004_RL_Repeat_processed.csv"
)
online_dat = pd.read_csv(online_dat_filepath, index_col=0)

online_dat["Biolector row"] = online_dat["Biolector well"].str.extract(
    r"([A-Z])"
)
online_dat["Biolector column"] = online_dat["Biolector well"].str.extract(
    r"(\d.*)"
)

# Rename columns
online_dat = online_dat.rename(columns={"Biomass [LS Gain=3]": "Biomass concentration [light scatter]"})

# Remove measurements after feeding has stopped
online_dat = online_dat[(online_dat["Feeding time"] > 0) & (online_dat["Feeding"] == True)]

# loading sampling data
sample_dat_filepath = (
    PROJECT_DIR
    / "article"
    / "real_world_data"
    / "Ct000259.LG_Pipetting_processed.csv"
)
sample_dat = pd.read_csv(sample_dat_filepath, index_col=0)

# merge datasets
df = online_dat.merge(
    sample_dat, on=["Biolector well", "Cycle"], how="left", validate="1:1"
)
df["Sample volume"] = df["Sample volume"].fillna(0)

# Calculate before sampling volume
df['Volume'] = df['Volume'] + df['Sample volume']

# Add "metadata"
## Here we will assign placeholder strain names because the real strain names are not important
strain_biolectorcolumn_map = {
    "01": "Yeast_strain_1",
    "02": "Yeast_strain_1",
    "03": "Yeast_strain_2",
    "04": "Yeast_strain_2",
    "05": "Yeast_strain_3",
    "06": "Yeast_strain_3",
    "07": "Yeast_strain_4",
    "08": "Yeast_strain_4",
}
df["Strain"] = df["Biolector column"].apply(
    lambda x: strain_biolectorcolumn_map[x]
)

# Keep only required columns
df = df[
    [
        "Strain",
        "Biolector well",
        "Biolector column",
        "Biolector row",
        "Cycle",
        "Time",
        "Feeding time",
        "Volume",
        "Accum. feed [uL]",
        "Biomass concentration [light scatter]",
        "DO%",
        "Sample volume",
    ]
]

# Save
df.to_csv(
    PROJECT_DIR
    / "article"
    / "real_world_data"
    / "biolector_yeast_fedbatch.csv",
    index=False,
)
