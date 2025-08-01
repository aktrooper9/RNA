import pandas as pd
import re

# --- Load lines ---
with open("gse115828_series_matrix.txt", "r", encoding="utf-8") as f:
    lines = [line.strip() for line in f if line.startswith("!")]

# --- Extract sample IDs ---
sample_line = next(line for line in lines if line.startswith("!Series_sample_id"))
sample_ids = sample_line.split('\t')[1].strip('"').split()
num_samples = len(sample_ids)

# --- Extract characteristics lines ---
char_lines = [line for line in lines if line.startswith("!Sample_characteristics_ch1")]

# --- Build metadata dictionary ---
meta_dict = {"sample_id": sample_ids}

for i, line in enumerate(char_lines):
    tokens = line.split('\t')[1:]  # remove field label
    clean_tokens = []
    for token in tokens:
        try:
            token = token.strip('"')
            key, value = re.split(r":\s*", token, maxsplit=1)
            clean_tokens.append((key, value))
        except Exception:
            # fallback for broken entries
            clean_tokens.append((f"Field_{i}", None))

    keys = [k for k, _ in clean_tokens]
    values = [v for _, v in clean_tokens]

    # Use the most frequent key name as column name
    key_name = max(set(keys), key=keys.count)
    meta_dict[key_name] = values

# --- Create DataFrame ---
sample_metadata = pd.DataFrame(meta_dict)

# --- Display first few rows ---
print(sample_metadata.head())

sample_metadata.to_csv("sample_metadata_gse115828.csv", index=False)

