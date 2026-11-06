# CodingProblems

Practice solutions for sequence alignment and classic algorithm problems.

## Problems

### Alignment (`Alignment.py`)

BMI 214 coursework: Needleman–Wunsch-style global and Smith–Waterman-style local pairwise alignment with affine gap penalties, match matrix input, and recursive traceback (including forked equal-score paths).

```bash
python Alignment.py <input_file> <output_file>
```

Input file format (line-based): sequences A/B, alignment mode (`0` global / `1` local), gap parameters `dx ex dy ey`, alphabets, then match-matrix rows. Output is the score plus alignment string pairs.

Requires: `numpy`. Note: the script uses Python 2-era `print` statements in places; prefer a Python 2 environment or a light port if running on modern Python.

### Median of two sorted arrays (`MedianSortedArrays.py`)

LeetCode-style `Solution.findMedianSortedArrays` that concatenates both arrays, sorts, and returns the median.

### Perfect square (`PerfectSquare.py`)

LeetCode-style `Solution.isPerfectSquare` that uses Newton’s method–style iteration to test whether an integer is a perfect square.

## Layout

```
CodingProblems/
├── Alignment.py
├── MedianSortedArrays.py
├── PerfectSquare.py
└── README.md
```
