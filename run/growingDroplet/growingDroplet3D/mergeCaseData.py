# This Python snippet uses Pandas to combine contact angle data
# and wetted area data into a single data table and write it
# as a proper CSV file.
# NOTE: this file is not intended to be used on its own,
# but rather to be called from 'combineCaseData.sh'.
import pandas as pd

wadf = pd.read_csv("postProcessing/wettedArea/0/surfaceFieldValue.dat",
                    delim_whitespace=True, comment='#', header=None,
                    names=["time", "wetted_area"])
cadf = pd.read_csv("contactAngleData.csv")

# Remove the dummy column that is only used as key for 'grep'
cadf = cadf.drop("dummyCol", axis=1)

# Combine data tables
cdf = wadf.join(cadf)

# Drop first row: data is at least partially invalid
cdf = cdf.drop(0)

cdf.to_csv("growing_droplet_data.csv", index=False)

