# pydebruijn-graph-join
*by Adamas Aqsa Fahreza and Martianus Frederic Ezerman*

A proof-of-concept implementation of the Graph Joining Prefer-Opposite algorithm found in:
> Z. Chang, M.F. Ezerman, A.A. Fahreza, and Q. Wang, "A Graph Joining Greedy Approach to Binary de Bruijn Sequences." Available at (https://arxiv.org/pdf/2004.09810)

Requires the package [pydebruijn](https://github.com/adamasstokhorst/pydebruijn) to be installed, which in turn requires Python 2.

## Usage

Simply run:
```
python treejoin.py
```

A section of the file may be directly modified for different inputs.

Alternatively, `treejoin.py` may be imported. Its main functionality can be accessed through the `join_trees_from_anf()` function. The following lines can be included in a separate `.py` file in the same directory:

```python
import treejoin

for anf in treejoin.join_trees_from_anf(...):
    ...
```

or

```python
from treejoin import join_trees_from_anf

for anf in join_trees_from_anf(...):
	...
```

For details on the `anf` object emitted by `join_trees_from_anf()`, refer to the `pydebruijn` package.
