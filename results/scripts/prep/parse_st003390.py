from pathlib import Path

input_file = Path("results/data/raw/ST003390_AN005559.txt")
output_file = Path("results/data/processed_metabolite_matrix_ST003390.csv")

lines = input_file.read_text(encoding="utf-8", errors="ignore").splitlines()

start = None
end = None

for i, line in enumerate(lines):
    if line.strip() == "MS_METABOLITE_DATA_START":
        start = i + 1
    elif line.strip() == "MS_METABOLITE_DATA_END":
        end = i
        break

if start is None or end is None:
    raise ValueError("Could not find MS_METABOLITE_DATA_START / END block")

data_block = lines[start:end]

output_file.parent.mkdir(parents=True, exist_ok=True)
output_file.write_text("\n".join(data_block), encoding="utf-8")

print(f"Saved: {output_file}")
