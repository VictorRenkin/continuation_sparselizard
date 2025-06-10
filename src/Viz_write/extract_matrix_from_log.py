import numpy as np
import re
import os

def parse_matrix_blocks(file_path, output_dir):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    os.makedirs(output_dir, exist_ok=True)

    skip_ligne = 3
    found_mat = False
    number_of_blocks = 0

    for line in lines:
        match = re.match(r"[A-Z] block (\d+)x(\d+):", line.strip())
        if match:
            rows, cols = int(match[1]), int(match[2])
            matrix = np.zeros((rows, cols))
            found_mat = True
        if found_mat:
            if skip_ligne > 0:
                skip_ligne -= 1
                continue
            else:
                if not line.startswith("row"):
                    pass
            row_match = re.match(r"row\s+(\d+):", line)
            
            if row_match:
                row_idx = int(row_match.group(1))  # ex: 173

                # Enlève le "row xxx:" pour ne garder que les couples
                content = line.split(":", 1)[1]

                # Cherche tous les couples (col, val)
                matches = re.findall(r"\((\d+),\s*([^\)]+)\)", content)

                for col_str, val_str in matches:
                    col = int(col_str)
                    val = float(val_str)
                    matrix[row_idx, col] = val
                if row_idx == rows - 1:
                    found_mat = False
                    np.savetxt(f"{output_dir}/Jacobienne_matrix_{number_of_blocks}.csv", matrix, delimiter=",")
                    print(f"MATRIX {rows}x{cols} :\n", matrix)
                    number_of_blocks += 1
                    skip_ligne = 3

    
parse_matrix_blocks("log.txt", "matrices_csv")
