import { show_alignment } from "./_show.js";
import { GAP, GAP_CHAR, SUBSTITUTION_MATRIX, UP, LEFT, DIAG } from "./_params.js";

// main
function main() {
  let [s1, s2] = process.argv.slice(2);
  let [s, a1, a2] = needleman_wunsch(s1, s2);
  console.log(show_alignment(a1, a2));
}

function needleman_wunsch(s1, s2, substitution_matrix = SUBSTITUTION_MATRIX, gap_penalty = GAP, gap_char = GAP_CHAR) {
  let [score_matrix, traceback_matrix] = compute_score_and_traceback_matrices(s1, s2, substitution_matrix, gap_penalty);
  let score = score_matrix[s1.length][s2.length];
  let [aln1, aln2] = traceback(traceback_matrix, s1, s2, gap_char);
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

// from the traceback matrix, resolve one of the possible alignments
function traceback(traceback_matrix, s1, s2, gap_char) {
  [s1, s2] = [Array.from(s1), Array.from(s2)];
  let [i, j] = [s1.length, s2.length];
  let [aln1, aln2] = [[], []];
  while (i > 0 && j > 0) {
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
  return [aln1.reverse().join(""), aln2.reverse().join("")];
}

// initialise the score matrix and the traceback matrix
function init_matrices([nrow, ncol], gap_penalty) {
  let score_matrix = Array.from({ length: nrow }, (e) => Array(ncol).fill(0));
  let traceback_matrix = Array.from({ length: nrow }, (e) => Array(ncol).fill(0));
  for (let j = 1; j < ncol; j++) {
    score_matrix[0][j] = score_matrix[0][j - 1] + gap_penalty;
    traceback_matrix[0][j] = LEFT;
  }
  for (let i = 1; i < nrow; i++) {
    score_matrix[i][0] = score_matrix[i - 1][0] + gap_penalty;
    traceback_matrix[i][0] = UP;
  }
  return [score_matrix, traceback_matrix];
}

// compute the max score of one grid and one of the possible directions
function compute_max_score_and_direction(up, left, diag) {
  let max_score = Math.max(up, left, diag);
  let direction = up === max_score ? UP : left === max_score ? LEFT : DIAG;
  return [max_score, direction];
}

main();
