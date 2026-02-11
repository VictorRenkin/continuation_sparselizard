# Test Cases

This directory contains all benchmark configurations used to validate the nonlinear continuation framework.

Each test case is fully reproducible and can be executed using the corresponding NLFR or NNM scripts from the main project.

The benchmarks cover increasing model complexity, including stationary and rotating structures.

---

## Summary of Test Cases

| Structure                  | DOFs    | Method | Node Partition | Points | Newton Iter. | Time       |
|---------------------------|---------|--------|----------------|--------|--------------|------------|
| Clamped–clamped beam       | 12,375  | NLFR   | 1              | 137    | 1,232        | 0:41:07    |
|                            |         | NNM    | 1              | 56     | 218          | 0:27:32    |
| Cantilever beam            | 227,205 | NLFR   | 4              | 149    | 1,467        | 46:21:25   |
|                            |         | NNM    | 1              | 24     | 309          | 7:19:33    |
| Cantilever beam (rotation) | 50,325  | NLFR   | 1              | ≈272   | ≈40          | ≈1:13:42   |
|                            |         | NNM    | 1              | ≈46    | ≈13          | ≈0:28:53   |
| Fan blade                  | 81,795  | NLFR   | 1              | 29     | 260          | 2:20:06    |
|                            |         | NNM    | 1              | 9      | 41           | 0:24:45    |

---

