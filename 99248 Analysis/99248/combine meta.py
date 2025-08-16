files = [

    "GSE99248_series_matrix 2.txt"
]

geo_accessions = []
studies = []
conditions = []
log = []

for fname in files:
    study = fname.split("_")[0]
    with open(fname, 'r') as f:
        lines = f.readlines()
        geo_line = None
        char_lines = []
        for line in lines:
            if line.startswith("!Sample_geo_accession") or line.startswith("!sample_geo_accession"):
                geo_line = [x.replace('"', '').strip() for x in line.strip().split('\t')[1:]]
            if line.startswith("!Sample_characteristics_ch"):
                char_lines.append([x.replace('"', '').strip() for x in line.strip().split('\t')[1:]])
        if geo_line is None or not char_lines:
            log.append(f"{fname}: Missing geo_accession or characteristics lines")
            continue
        for i, geo in enumerate(geo_line):
            cond = None
            for char_line in char_lines:
                if i >= len(char_line):
                    continue
                char = char_line[i]
                if ':' in char:
                    key, val = char.split(':', 1)
                    key = key.strip().lower()
                    val = val.strip()
                else:
                    key = ''
                    val = char.strip()
                # Normalize condition
                if key == 'donor id':
                    if val.lower().startswith('normal'):
                        cond = '1'
                        break
                    elif val.lower().startswith('amd'):
                        cond = '4'
                        break
                elif key == 'mgs_level':
                    if val in {'1', '2', '3', '4'}:
                        cond = val
                        break
                elif key == 'amd_status':
                    if val.lower() == 'control':
                        cond = '1'
                        break
                    elif val.upper() == 'AMD':
                        cond = '4'
                        break
            if cond is not None:
                geo_accessions.append(geo)
                studies.append(study)
                conditions.append(cond)
            else:
                log.append(f"{fname}: {geo} - No recognized condition in any characteristics line")

# Remove duplicates while preserving order
seen = set()
unique_geo = []
unique_study = []
unique_cond = []
for geo, study, cond in zip(geo_accessions, studies, conditions):
    if geo and geo not in seen:
        unique_geo.append(geo)
        unique_study.append(study)
        unique_cond.append(cond)
        seen.add(geo)

# Write output
with open("combined_geo_accession_study_condition.tsv", "w") as out:
    out.write("sample_geo_accession\tstudy\tcondition\n")
    for geo, study, cond in zip(unique_geo, unique_study, unique_cond):
        out.write(f"{geo}\t{study}\t{cond}\n")

# Write log
if log:
    with open("combine_meta_condition.log", "w") as lfile:
        for entry in log:
            lfile.write(entry + "\n")

print("Done. See combined_geo_accession_study_condition.tsv and combine_meta_condition.log")
