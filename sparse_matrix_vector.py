import numpy as np

def load_sparse_vectors_from_file(filepath):
    """
    Reads a text file containing sparse matrix elements.
    Expected line format: row col real_part imag_part

    Example: 
        % row  col  real   imag
        1    1   10.5  -40.2
        1    2   -2.1    5.0
        2    1   -2.1    5.0
        2    2    8.4  -30.1
        2    3   -3.3    8.2
        3    2   -3.3    8.2
        3    3    6.1  -20.0
    """
    rows = []
    cols = []
    values = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Ignore empty lines or comments
            if not line or line.startswith('%'):
                continue
            
            tokens = line.split()
            if len(tokens) < 3:
                continue
            
            r = int(tokens[0])
            c = int(tokens[1])
            
            # Parse value (supports separate real/imag columns or single complex string)
            if len(tokens) >= 4:
                val = complex(float(tokens[2]), float(tokens[3]))
            else:
                val = complex(tokens[2])
                
            rows.append(r)
            cols.append(c)
            values.append(val)
            
    rows = np.array(rows)
    cols = np.array(cols)
    values = np.array(values)
    
    # Determine matrix dimension N (assuming 1-based indexing as busses 
    # in this class have usually started at 1)
    N = max(rows.max(), cols.max())
    
    diag_vec = np.zeros(N, dtype=complex)
    offdiag = []
    col_list = []
    row_list = []
    
    # Separate diagonal and off-diagonal elements
    for r, c, val in zip(rows, cols, values):
        if r == c:
            diag_vec[r - 1] = val  # 1-based index adjustment
        else:
            row_list.append(r)
            col_list.append(c)
            offdiag.append(val)
            
    return diag_vec, np.array(offdiag), np.array(col_list), np.array(row_list), N

# Example usage:
if __name__ == "__main__":
    # DIAG, OFFDIAG, COL, ROW, N = load_sparse_vectors_from_file("ybus_elements.txt")
    pass