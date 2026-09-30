# abaqus-inp-to-freefem

A lightweight, dependency-free Python CLI tool to convert Abaqus 3D solid tetrahedral meshes (`.inp`) to the INRIA Medit (`.mesh`) format, natively supported by **FreeFem++** via `readmesh3`.

---

## Features

- **No external dependencies**: Pure Python 3 (requires standard library only).
- **Geometric Boundary Reconstruction**: Reconstructs external surface triangles automatically from 3D tetrahedral connectivity.
- **Parametric Label Mapping**: Maps Abaqus `*Nset` (Node Sets) and `*Elset` (Element Sets) to integer surface and volume region labels.
- **Robust Abaqus Parsing**: Parses both explicit comma-separated node/element lists and blocks utilizing the Abaqus `generate` syntax.
- **FreeFem++ Ready**: Guarantees full compatibility with `readmesh3()` and `medit()` in FreeFem++.

---

## Installation

Clone the repository and ensure you have Python 3 installed:

```bash
git clone [https://github.com/your-username/abaqus-inp-to-freefem.git](https://github.com/your-username/abaqus-inp-to-freefem.git)
cd abaqus-inp-to-freefem