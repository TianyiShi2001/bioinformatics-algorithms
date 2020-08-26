import os

import re

chars = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
    "B",
    "Z",
    "X",
    "*",
]

with open("a.txt", "r") as f:
    vals = re.split(r"\s+", f.read())[1:-1]

vals.reverse()

matrix = {}
for i in chars:
    matrix[i] = {}
    for j in chars:
        matrix[i][j] = int(vals.pop())

import json

with open("blosum45.json", "w") as f:
    json.dump(matrix, f)

