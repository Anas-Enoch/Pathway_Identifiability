import pandas as pd

# read HTML tables
tables = pd.read_html("ST000356-metabolite data.php.html")

# usually the first large table contains the metabolite matrix
df = tables[0]

print(df.head())

# save as csv
df.to_csv("data/processed_metabolite_matrix_ST000356.csv", index=False)
