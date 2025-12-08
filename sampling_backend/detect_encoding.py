import chardet

filename = "taxonomy_all_new.jsonl"

with open(filename, "rb") as f:
    raw = f.read(200000)   # read first 200 KB
    print(chardet.detect(raw))
