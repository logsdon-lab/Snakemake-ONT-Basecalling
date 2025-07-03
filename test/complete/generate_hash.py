

import hashlib 


with open("runs.txt") as fh:
    hasher = hashlib.sha256()
    for line in fh:
        line = '_'.join(line.strip().split("/"))
        hsh = hashlib.sha256(line.encode()).hexdigest()
        with open(f"{hsh}.done", "wt") as ofh:
            ofh.write("")