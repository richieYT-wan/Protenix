# shard.py
import json, math, sys, pathlib

input_json = pathlib.Path(sys.argv[1])
n = int(sys.argv[2])
data = json.loads(input_json.read_text())
size = math.ceil(len(data) / n)
for i in range(n):
    out = input_json.with_name(f"{input_json.stem}_shard{i}.json")
    out.write_text(json.dumps(data[i*size:(i+1)*size]))
    print(out)