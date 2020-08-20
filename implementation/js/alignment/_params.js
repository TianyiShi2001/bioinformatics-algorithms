const MATCH = 2;
const MISMATCH = -1;
export const GAP = -1;
export const GAP_CHAR = "-";

export const SUBSTITUTION_MATRIX = {
  A: { A: MATCH, C: MISMATCH, G: MISMATCH, T: MISMATCH },
  C: { A: MISMATCH, C: MATCH, G: MISMATCH, T: MISMATCH },
  G: { A: MISMATCH, C: MISMATCH, G: MATCH, T: MISMATCH },
  T: { A: MISMATCH, C: MISMATCH, G: MISMATCH, T: MATCH },
};

// Set the constants that represent the three directions, which are used in the traceback matrix
export const UP = 1;
export const LEFT = 2;
export const DIAG = 3;
