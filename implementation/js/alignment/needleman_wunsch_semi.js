import { GAP, GAP_CHAR, SUBSTITUTION_MATRIX, UP, LEFT, DIAG } from "./_params.js";

export function needleman_wunsch_semi(s1, s2, substitution_matrix = SUBSTITUTION_MATRIX, gap_penalty = GAP, gap_char = GAP_CHAR) {
  let [score_matrix, traceback_matrix] = compute_score_and_traceback_matrices(s1, s2, substitution_matrix, gap_penalty);
  let [score, coords] = find_max_score_and_coords(score_matrix);
  let [aln1, aln2] = traceback(traceback_matrix, s1, s2, gap_char, coords);
  return [score, aln1, aln2];
}

function compute_score_and_traceback_matrices(s1, s2, substitution_matrix, gap_penalty) {
  let [nrow, ncol] = [s1.length + 1, s2.length + 1];
  let [score_matrix, traceback_matrix] = init_matrices([nrow, ncol], gap_penalty);
  for (let i = 1; i < nrow; i++) {
    for (let j = 1; j < ncol; j++) {
      let up = score_matrix[i - 1][j] + gap_penalty;
      let left = score_matrix[i][j - 1] + gap_penalty;
      let diag = score_matrix[i - 1][j - 1] + substitution_matrix[s1[i - 1]][s2[j - 1]];
      let [max_score, directions] = compute_max_score_and_direction(up, left, diag);
      score_matrix[i][j] = max_score;
      traceback_matrix[i][j] = directions;
    }
  }
  return [score_matrix, traceback_matrix];
}

function find_max_score_and_coords(score_matrix) {
  let max = 0;
  let coords = [0, 0];
  let [nrow, ncol] = [score_matrix.length, score_matrix[0].length];
  for (let i = 0; i < nrow; i++) {
    const n = score_matrix[i][ncol - 1];
    if (n > max) {
      max = n;
      coords = [i, ncol - 1];
    }
  }
  for (let j = 0; j < ncol; j++) {
    const n = score_matrix[nrow - 1][j];
    if (n > max) {
      max = n;
      coords = [nrow - 1, j];
    }
  }
  return [max, coords];
}

// from the traceback matrix, resolve one of the possible alignments
function traceback(traceback_matrix, s1, s2, gap_char, coords) {
  [s1, s2] = [Array.from(s1), Array.from(s2)];
  let [i, j] = coords;
  let [aln1, aln2] = [[], []];

  for (let m = 1; m < traceback_matrix.length - i; m++) {
    aln1.push(s1.pop());
    aln2.push(gap_char);
  }
  for (let n = 1; n < traceback_matrix[0].length - j; n++) {
    aln1.push(gap_char);
    aln2.push(s2.pop());
  }

  while (i !== 0 && j !== 0) {
    let direction = traceback_matrix[i][j];
    if (direction === UP) {
      aln1.push(s1.pop());
      aln2.push(gap_char);
      i--;
    } else if (direction === LEFT) {
      aln1.push(gap_char);
      aln2.push(s2.pop());
      j--;
    } else {
      aln1.push(s1.pop());
      aln2.push(s2.pop());
      i--;
      j--;
    }
  }

  for (let m = 0; m < i; m++) {
    aln1.push(s1.pop());
    aln2.push(gap_char);
  }
  for (let n = 0; n < j; n++) {
    aln1.push(gap_char);
    aln2.push(s2.pop());
  }

  return [aln1.reverse().join(""), aln2.reverse().join("")];
}

// initialise the score matrix and the traceback matrix
function init_matrices([nrow, ncol], gap_penalty) {
  let score_matrix = Array.from({ length: nrow }, (e) => Array(ncol).fill(0));
  let traceback_matrix = Array.from({ length: nrow }, (e) => Array(ncol).fill(0));
  return [score_matrix, traceback_matrix];
}

// compute the max score of one grid and one of the possible directions
function compute_max_score_and_direction(up, left, diag) {
  let max_score = Math.max(up, left, diag);
  let direction = up === max_score ? UP : left === max_score ? LEFT : DIAG;
  return [max_score, direction];
}
