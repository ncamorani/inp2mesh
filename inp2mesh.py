#!/usr/bin/env python3
"""
Abaqus .inp to INRIA Medit .mesh Converter for FreeFem++
Author: Nicola Camorani
License: MIT

Description:
    This script converts 3D tetrahedral solid meshes from Abaqus input files (.inp)
    into the INRIA Medit (.mesh) format, which is natively readable by FreeFem++ via readmesh3().
    
Key Technical Features:
    1. Scope isolation: Restricts node and element extraction strictly to the *Part level,
       preventing coordinate corruption from global Assembly reference nodes.
    2. Topological boundary reconstruction: Reconstructs 2D external boundary triangles
       from 3D tetrahedra using hash-based face occurrence counters.
    3. Parametric set mapping: Maps Abaqus *Nset (Node Sets) and *Elset (Element Sets)
       to user-defined integer boundary and volume region labels.
    4. Orientation enforcement: Uses an automated chiral majority voting algorithm to
       detect and fix negatively oriented tetrahedra (det < 0), preventing FreeFem++ runtime crashes.
"""

import argparse
import sys
from typing import Dict, List, Set, Tuple


def parse_args():
    """
    Configure and parse command-line interface (CLI) arguments.
    """
    parser = argparse.ArgumentParser(
        description="Convert Abaqus 3D tetrahedral (.inp) meshes to INRIA Medit (.mesh) format for FreeFem++."
    )
    
    # Required positional / path arguments
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the input Abaqus .inp file."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path to the destination INRIA Medit .mesh file."
    )
    
    # Boundary condition mapping arguments (can be specified multiple times)
    parser.add_argument(
        "--bc-nset",
        action="append",
        metavar="NSET_NAME=LABEL",
        help="Map an Abaqus Node Set (Nset) to an external surface boundary triangle label. "
             "Example: --bc-nset _PickedSet7=2 (Can be repeated for multiple sets).",
    )
    
    # Volume element mapping arguments (can be specified multiple times)
    parser.add_argument(
        "--vol-elset",
        action="append",
        metavar="ELSET_NAME=LABEL",
        help="Map an Abaqus Element Set (Elset) to a 3D tetrahedral volume region label. "
             "Example: --vol-elset _PickedSet2=1 (Can be repeated for multiple sets).",
    )
    
    # Default fallback labels
    parser.add_argument(
        "--default-bc-label",
        type=int,
        default=1,
        help="Default boundary label for external faces not matching any mapped BC Nset (default: 1).",
    )
    parser.add_argument(
        "--default-vol-label",
        type=int,
        default=1,
        help="Default volume region label for tetrahedra not matching any mapped Elset (default: 1).",
    )
    
    return parser.parse_args()


def parse_mapping_args(arg_list: List[str]) -> Dict[str, int]:
    """
    Convert a list of key-value mapping strings (e.g. ['_PickedSet7=2', 'TipSet=3'])
    into a lowercase dictionary mapping set names to integer labels.
    """
    mapping = {}
    if not arg_list:
        return mapping
    for item in arg_list:
        if "=" not in item:
            sys.exit(f"Error: Invalid mapping format '{item}'. Expected format is NAME=LABEL (e.g., _PickedSet7=2)")
        name, label_str = item.split("=", 1)
        name = name.strip()
        try:
            label = int(label_str.strip())
        except ValueError:
            sys.exit(f"Error: Label value in mapping '{item}' must be an integer, got '{label_str}'.")
        mapping[name.lower()] = label
    return mapping


def read_abaqus_inp(filepath: str):
    """
    Parse an Abaqus .inp file line by line with robust scope handling.
    
    Crucial Design Note:
        Only nodes and elements defined inside the '*Part' block are parsed.
        Any nodes defined later at the '*Assembly' level (such as reference points,
        rotation axis markers, or kinematic constraints) are deliberately skipped
        to avoid overwriting real mesh node coordinates.
    """
    nodes: Dict[int, Tuple[float, float, float]] = {}
    elements: Dict[int, Tuple[int, int, int, int]] = {}
    nsets: Dict[str, Set[int]] = {}
    elsets: Dict[str, Set[int]] = {}

    current_mode = None
    current_set_name = None
    is_generate = False
    inside_part = False

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            clean = line.strip()
            
            # Skip blank lines and pure Abaqus comment lines starting with '**'
            if not clean or clean.startswith("**"):
                continue

            # Check for keyword header lines starting with '*'
            if clean.startswith("*"):
                header_parts = [p.strip().lower() for p in clean.split(",")]
                keyword = header_parts[0]

                # Identify if the current set utilizes range expansion (generate keyword)
                is_generate = "generate" in header_parts

                # State tracking: Part boundary scope
                if keyword == "*part":
                    inside_part = True
                    current_mode = None
                    continue
                elif keyword == "*end part":
                    inside_part = False
                    current_mode = None
                    continue

                # Read nodes strictly within the Part
                if keyword == "*node" and inside_part:
                    current_mode = "node"
                    current_set_name = None
                    
                # Read solid elements strictly within the Part
                elif keyword == "*element" and inside_part:
                    el_type = None
                    for part in header_parts:
                        if part.startswith("type="):
                            el_type = part.split("=")[1].strip()
                    if el_type and el_type not in ["c3d4", "c3d4h"]:
                        print(f"Warning: Element type '{el_type}' detected. Optimized for linear tetrahedra (C3D4).")
                    current_mode = "element"
                    current_set_name = None
                    
                # Read node sets (can exist inside Part or at Assembly/Instance level)
                elif keyword == "*nset":
                    current_mode = "nset"
                    for part in clean.split(",")[1:]:
                        if "nset=" in part.lower():
                            current_set_name = part.split("=")[1].strip().lower()
                    if current_set_name and current_set_name not in nsets:
                        nsets[current_set_name] = set()
                        
                # Read element sets (can exist inside Part or at Assembly/Instance level)
                elif keyword == "*elset":
                    current_mode = "elset"
                    for part in clean.split(",")[1:]:
                        if "elset=" in part.lower():
                            current_set_name = part.split("=")[1].strip().lower()
                    if current_set_name and current_set_name not in elsets:
                        elsets[current_set_name] = set()
                else:
                    # Reset mode for unhandled keywords (e.g. *Material, *Boundary, *Step)
                    current_mode = None
                    current_set_name = None
                continue

            # Split comma-delimited data payload
            tokens = [p.strip() for p in clean.split(",") if p.strip()]
            if not tokens:
                continue

            # Store node coordinates: node_id -> (X, Y, Z)
            if current_mode == "node" and inside_part:
                node_id = int(tokens[0])
                coords = (float(tokens[1]), float(tokens[2]), float(tokens[3]))
                nodes[node_id] = coords

            # Store tetrahedral element connectivity: elem_id -> (n1, n2, n3, n4)
            elif current_mode == "element" and inside_part:
                elem_id = int(tokens[0])
                connectivity = (int(tokens[1]), int(tokens[2]), int(tokens[3]), int(tokens[4]))
                elements[elem_id] = connectivity

            # Store node indices belonging to an active *Nset
            elif current_mode == "nset" and current_set_name:
                if is_generate:
                    # Syntax: start_id, end_id, increment
                    start_id, end_id, step = int(tokens[0]), int(tokens[1]), int(tokens[2])
                    nsets[current_set_name].update(range(start_id, end_id + 1, step))
                else:
                    nsets[current_set_name].update(int(t) for t in tokens)

            # Store element indices belonging to an active *Elset
            elif current_mode == "elset" and current_set_name:
                if is_generate:
                    # Syntax: start_id, end_id, increment
                    start_id, end_id, step = int(tokens[0]), int(tokens[1]), int(tokens[2])
                    elsets[current_set_name].update(range(start_id, end_id + 1, step))
                else:
                    elsets[current_set_name].update(int(t) for t in tokens)

    return nodes, elements, nsets, elsets


def fix_tetrahedra_orientation(
    nodes: Dict[int, Tuple[float, float, float]],
    elements: Dict[int, Tuple[int, int, int, int]]
) -> Tuple[Dict[int, Tuple[int, int, int, int]], int]:
    """
    Enforce chiral orientation consistency across all tetrahedra.
    
    Theory:
        FreeFem++ strictly requires every tetrahedral element to exhibit a positive
        volume Jacobian:
            det( [p1 - p4,  p2 - p4,  p3 - p4] ) > 0
            
        Occasionally, meshing algorithms produce locally inverted elements (negative determinant),
        which triggers:
            'Fatal Error: number badly oriented element X'
            
        This function computes the signed scalar triple product for every element.
        Using a majority-voting scheme, it detects the prevailing chirality convention
        and swaps vertex 2 and vertex 3 exclusively for minority outlier elements,
        restoring strict positivity without distorting element geometry.
    """
    dets = {}
    pos_count = 0
    neg_count = 0

    for elem_id, (n1, n2, n3, n4) in elements.items():
        p1 = nodes[n1]
        p2 = nodes[n2]
        p3 = nodes[n3]
        p4 = nodes[n4]

        # Construct direction vectors relative to apex node p4
        ax, ay, az = p1[0] - p4[0], p1[1] - p4[1], p1[2] - p4[2]
        bx, by, bz = p2[0] - p4[0], p2[1] - p4[1], p2[2] - p4[2]
        cx, cy, cz = p3[0] - p4[0], p3[1] - p4[1], p3[2] - p4[2]

        # Calculate 3x3 determinant (signed volumetric determinant)
        det = (
            ax * (by * cz - bz * cy)
            - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx)
        )
        dets[elem_id] = det
        if det > 0:
            pos_count += 1
        else:
            neg_count += 1

    # Identify target chirality favored by the vast majority of elements
    target_positive = (pos_count >= neg_count)

    fixed_elements = {}
    inverted_count = 0

    for elem_id, conn in elements.items():
        det = dets[elem_id]
        is_pos = (det > 0)

        # Invert only elements that oppose the dominant majority
        if is_pos != target_positive:
            n1, n2, n3, n4 = conn
            fixed_elements[elem_id] = (n1, n3, n2, n4)  # Swap nodes 2 and 3 to reverse sign
            inverted_count += 1
        else:
            fixed_elements[elem_id] = conn

    return fixed_elements, inverted_count


def extract_boundary_faces(elements: Dict[int, Tuple[int, int, int, int]]) -> List[Tuple[int, int, int]]:
    """
    Extract exterior triangular boundary faces from 3D tetrahedral connectivity.
    
    Topological Algorithm:
        1. Deconstruct every 4-node tetrahedron into its 4 triangular constituent faces.
        2. Sort vertex indices of each face to create an order-invariant canonical key:
               canonical_key = tuple(sorted(face))
        3. Internal faces are shared between exactly two adjacent tetrahedra (frequency = 2).
        4. Exterior boundary faces belong to only one tetrahedron (frequency = 1).
        5. Retain only faces with an occurrence count of 1.
    """
    face_count: Dict[Tuple[int, int, int], Tuple[int, int, int]] = {}

    for conn in elements.values():
        a, b, c, d = conn
        tet_faces = [
            (a, b, c),
            (a, b, d),
            (a, c, d),
            (b, c, d),
        ]
        for face in tet_faces:
            canonical_key = tuple(sorted(face))
            if canonical_key not in face_count:
                face_count[canonical_key] = face  # First encounter: potentially boundary
            else:
                face_count[canonical_key] = None  # Second encounter: internal shared face

    # Filter out None values to keep only unique surface triangles
    return [face for face in face_count.values() if face is not None]


def write_medit_mesh(
    output_path: str,
    nodes: Dict[int, Tuple[float, float, float]],
    elements: Dict[int, Tuple[int, int, int, int]],
    boundary_faces: List[Tuple[int, int, int]],
    nset_to_bc: Dict[str, int],
    nsets: Dict[str, Set[int]],
    elset_to_vol: Dict[str, int],
    elsets: Dict[str, Set[int]],
    default_bc_label: int,
    default_vol_label: int,
):
    """
    Write the mesh in INRIA Medit (.mesh) format.
    
    Structure:
        - MeshVersionFormatted 1 / Dimension 3
        - Vertices:  X  Y  Z  vertex_label
        - Triangles: n1 n2 n3 boundary_label
        - Tetrahedra: n1 n2 n3 n4 volume_region_label
    """
    # Create a 1-based continuous index mapping for FreeFem++ compatibility
    sorted_node_ids = sorted(nodes.keys())
    id_map = {orig_id: new_idx for new_idx, orig_id in enumerate(sorted_node_ids, start=1)}

    # Filter active boundary node sets and volume element sets
    active_bc_sets = [(label, nsets[name]) for name, label in nset_to_bc.items() if name in nsets]
    active_vol_sets = [(label, elsets[name]) for name, label in elset_to_vol.items() if name in elsets]

    bc_stats: Dict[int, int] = {}
    vol_stats: Dict[int, int] = {}

    with open(output_path, "w", encoding="utf-8") as f:
        # File headers
        f.write("MeshVersionFormatted 1\n")
        f.write("Dimension 3\n\n")

        # Section 1: Vertices (Nodes)
        f.write("Vertices\n")
        f.write(f"{len(nodes)}\n")
        for orig_id in sorted_node_ids:
            x, y, z = nodes[orig_id]
            f.write(f"{x:.10e} {y:.10e} {z:.10e} 1\n")

        # Section 2: Triangles (2D Exterior Boundary Faces)
        f.write("\nTriangles\n")
        f.write(f"{len(boundary_faces)}\n")
        for face in boundary_faces:
            n1, n2, n3 = face
            face_label = default_bc_label

            # If all 3 vertices of the face belong to an active BC Nset, assign that label
            for label, node_set in active_bc_sets:
                if n1 in node_set and n2 in node_set and n3 in node_set:
                    face_label = label
                    break

            bc_stats[face_label] = bc_stats.get(face_label, 0) + 1
            f.write(f"{id_map[n1]} {id_map[n2]} {id_map[n3]} {face_label}\n")

        # Section 3: Tetrahedra (3D Solid Elements)
        f.write("\nTetrahedra\n")
        f.write(f"{len(elements)}\n")
        for elem_id, conn in elements.items():
            vol_label = default_vol_label
            
            # Check if this element belongs to a mapped volume Elset
            for label, el_set in active_vol_sets:
                if elem_id in el_set:
                    vol_label = label
                    break

            vol_stats[vol_label] = vol_stats.get(vol_label, 0) + 1
            e1, e2, e3, e4 = conn
            f.write(f"{id_map[e1]} {id_map[e2]} {id_map[e3]} {id_map[e4]} {vol_label}\n")

        # File termination keyword
        f.write("\nEnd\n")

    return bc_stats, vol_stats


def main():
    """
    Main driver pipeline: Argument parsing -> Extraction -> Orientation Fix -> Boundary Detect -> Medit Export.
    """
    args = parse_args()
    bc_mapping = parse_mapping_args(args.bc_nset)
    vol_mapping = parse_mapping_args(args.vol_elset)

    # 1. Parse Abaqus .inp file
    print(f"Reading Abaqus input: {args.input}")
    nodes, elements, nsets, elsets = read_abaqus_inp(args.input)

    print(f"Nodes loaded: {len(nodes)}")
    print(f"Elements loaded (C3D4): {len(elements)}")
    print(f"Node sets found: {len(nsets)} | Element sets found: {len(elsets)}")

    # 2. Check and correct tetrahedral orientations
    print("Checking and correcting tetrahedra orientations (FreeFem++ positive volume rule)...")
    elements, inverted_count = fix_tetrahedra_orientation(nodes, elements)
    if inverted_count > 0:
        print(f"-> Corrected {inverted_count} badly oriented tetrahedra.")

    # 3. Extract external boundary surface triangles
    print("Reconstructing external boundary surfaces...")
    boundary_faces = extract_boundary_faces(elements)
    print(f"Total boundary triangles identified: {len(boundary_faces)}")

    # 4. Write output Medit .mesh file
    print(f"Writing INRIA Medit output: {args.output}")
    bc_stats, vol_stats = write_medit_mesh(
        output_path=args.output,
        nodes=nodes,
        elements=elements,
        boundary_faces=boundary_faces,
        nset_to_bc=bc_mapping,
        nsets=nsets,
        elset_to_vol=vol_mapping,
        elsets=elsets,
        default_bc_label=args.default_bc_label,
        default_vol_label=args.default_vol_label,
    )

    # 5. Output summary diagnostics
    print("\nSummary of Exported Labels:")
    print("----------------------------")
    print("Boundary Triangles (Surfaces):")
    for lbl, count in sorted(bc_stats.items()):
        print(f"  - Label {lbl}: {count} faces")
    print("Tetrahedra (Volumes):")
    for lbl, count in sorted(vol_stats.items()):
        print(f"  - Region {lbl}: {count} elements")
    print(f"\nConversion successfully completed -> {args.output}")


if __name__ == "__main__":
    main()