import os
import sys

import numpy as np
import pandas as pd

OUTPUT_FILE = os.path.join(
    "tests", "test_data", "fed-batch3_measurements_only.csv"
)

fedbatch_file = os.path.join("tests", "test_data", "fed-batch3.csv")
fedbatch_df = pd.read_csv(fedbatch_file, index_col=0)

measurement_times = [10, 30, 40, 50]
exact_measurements = (
    fedbatch_df.query("timestamp.isin(@measurement_times)")
    .drop_duplicates(
        subset="timestamp", keep="last"
    )  # ODE solver save both data before and after sampling event
    .fillna({"sample_volume": 0})
)
exact_measurements.to_csv(OUTPUT_FILE)
