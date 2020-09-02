use crate::alignment::common::*;
use ndarray::prelude::*;
use ndarray::Array;
use std::env;
pub type Score = i32;
pub type Direction = u8; // (0, 1, 2, 3) = (None, Up, Left, Diag)
pub type Matrix<T> = Array<T, Ix2>;
pub type ScoreMatrix = Matrix<Score>;
pub type TracebackMatrix = Matrix<Direction>;
pub type Seq = String;
pub struct Coords(pub usize, pub usize);

pub const MATCH: Score = 2;
pub const MISMATCH: Score = -1;
pub const GAP: Score = -1;

fn main() {
    let query: Vec<String> = env::args().skip(1).take(2).collect();
    let (score, aln1, aln2) = needleman_wunsch(&query[0], &query[1], '-');
    println!("{}\n{}\n{}", score, aln1, aln2);
}

fn init_score_and_traceback_matrices(nrow: usize, ncol: usize) -> (ScoreMatrix, TracebackMatrix) {
    let mut score_matrix: ScoreMatrix = Array::zeros((nrow, ncol));
    let mut traceback_matrix: TracebackMatrix = Array::zeros((nrow, ncol));
    for i in 1..nrow {
        score_matrix[[i, 0]] = score_matrix[[i - 1, 0]] + GAP;
        traceback_matrix[[i, 0]] = 1;
    }
    for j in 1..ncol {
        score_matrix[[0, j]] = score_matrix[[0, j - 1]] + GAP;
        traceback_matrix[[0, j]] = 2;
    }
    return (score_matrix, traceback_matrix);
}

fn compute_substitution_score(a: char, b: char) -> Score {
    if a == b {
        MATCH
    } else {
        MISMATCH
    }
}

fn compute_max_score_and_direction(up: Score, left: Score, diag: Score) -> (Score, Direction) {
    if left > up {
        if diag > left {
            return (diag, 3);
        }
        return (left, 2);
    } else {
        if diag > up {
            return (diag, 3);
        }
        return (up, 1);
    }
}

fn compute_score_and_traceback_matrices(s1: &str, s2: &str) -> (ScoreMatrix, TracebackMatrix) {
    let (nrow, ncol) = (s1.len() + 1, s2.len() + 1);
    let (mut score_matrix, mut traceback_matrix) = init_score_and_traceback_matrices(nrow, ncol);
    for i in 1..nrow {
        for j in 1..ncol {
            let up = score_matrix[[i - 1, j]] + GAP;
            let left = score_matrix[[i, j - 1]] + GAP;
            let diag = score_matrix[[i - 1, j - 1]]
                + compute_substitution_score(
                    s1.chars().nth(i - 1).unwrap(),
                    s2.chars().nth(j - 1).unwrap(),
                );
            let (max, dir) = compute_max_score_and_direction(up, left, diag);
            score_matrix[[i, j]] = max;
            traceback_matrix[[i, j]] = dir;
        }
    }
    return (score_matrix, traceback_matrix);
}

fn traceback(traceback_matrix: TracebackMatrix, s1: &str, s2: &str, gap_char: char) -> (Seq, Seq) {
    let mut s1: Vec<char> = s1.chars().collect();
    let mut s2: Vec<char> = s2.chars().collect();
    let (mut i, mut j) = (s1.len(), s2.len());
    let mut aln1: Vec<char> = Vec::with_capacity(i);
    let mut aln2: Vec<char> = Vec::with_capacity(j);
    while (i > 0) && (j > 0) {
        let dir = traceback_matrix[[i, j]];
        match dir {
            1 => {
                i = i - 1;
                aln1.push(s1.pop().unwrap());
                aln2.push(gap_char);
            }
            2 => {
                j = j - 1;
                aln1.push(gap_char);
                aln2.push(s2.pop().unwrap());
            }
            3 => {
                i = i - 1;
                j = j - 1;
                aln1.push(s1.pop().unwrap());
                aln2.push(s2.pop().unwrap());
            }
            _ => {}
        }
    }
    return (
        aln1.into_iter().rev().collect(),
        aln2.into_iter().rev().collect(),
    );
}

pub fn needleman_wunsch(s1: &str, s2: &str, gap_char: char) -> (Score, Seq, Seq) {
    let (score_matrix, traceback_matrix) = compute_score_and_traceback_matrices(s1, s2);
    let score = score_matrix[[s1.len(), s2.len()]];
    let (aln1, aln2) = traceback(traceback_matrix, s1, s2, gap_char);
    (score, aln1, aln2)
}
