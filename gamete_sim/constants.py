"""Fixed reference values used throughout the pipeline."""

AUTOSOMES = tuple(range(1, 23))

# GRCh38 primary-assembly chromosome lengths (bp), GenBank GCA_000001405.15.
# These never change between GRCh38 patch releases. genetic_map.build() also
# checks them against the hg38 sizes recorded in the UCSC chain file.
GRCH38_LENGTHS = {
    1: 248956422,
    2: 242193529,
    3: 198295559,
    4: 190214555,
    5: 181538259,
    6: 170805979,
    7: 159345973,
    8: 145138636,
    9: 138394717,
    10: 133797422,
    11: 135086622,
    12: 133275309,
    13: 114364328,
    14: 107043718,
    15: 101991189,
    16: 90338345,
    17: 83257441,
    18: 80373285,
    19: 58617616,
    20: 64444167,
    21: 46709983,
    22: 50818468,
}

SEXES = ("male", "female")

# Namespaces for per-stream random seeds. Changing these changes every output.
SEX_SEED_KEY = {"male": 0, "female": 1}
