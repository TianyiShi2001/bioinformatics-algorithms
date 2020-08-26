use std::env;

use crate::alignment::common::*;

fn main() {
    let query: Vec<String> = env::args().skip(1).take(2).collect();
    let (score, aln1, aln2) = needleman_wunsch(&query[0], &query[1], '-');
    println!("{}\n{}\n{}", score, aln1, aln2);
}

fn init_score_and_traceback_matrices(nrow: usize, ncol: usize) -> (ScoreMatrix, TracebackMatrix) {
    let mut score_matrix = vec![vec![0i32; ncol]; nrow];
    let mut traceback_matrix = vec![vec![Direction::None; ncol]; nrow];
    for i in 1..nrow {
        score_matrix[i][0] = score_matrix[i - 1][0] + GAP;
        traceback_matrix[i][0] = Direction::Up;
    }
    for j in 1..ncol {
        score_matrix[0][j] = score_matrix[0][j - 1] + GAP;
        traceback_matrix[0][j] = Direction::Left;
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
            return (diag, Direction::Diag);
        }
        return (left, Direction::Left);
    } else {
        if diag > up {
            return (diag, Direction::Diag);
        }
        return (up, Direction::Up);
    }
}

fn compute_score_and_traceback_matrices(s1: &str, s2: &str) -> (ScoreMatrix, TracebackMatrix) {
    let (nrow, ncol) = (s1.len() + 1, s2.len() + 1);
    let (mut score_matrix, mut traceback_matrix) = init_score_and_traceback_matrices(nrow, ncol);
    for i in 1..nrow {
        for j in 1..ncol {
            let up = score_matrix[i - 1][j] + GAP;
            let left = score_matrix[i][j - 1] + GAP;
            let diag = score_matrix[i - 1][j - 1]
                + compute_substitution_score(
                    s1.chars().nth(i - 1).unwrap(),
                    s2.chars().nth(j - 1).unwrap(),
                );
            let (max, dir) = compute_max_score_and_direction(up, left, diag);
            score_matrix[i][j] = max;
            traceback_matrix[i][j] = dir;
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
        let dir = traceback_matrix[i][j];
        match dir {
            Direction::Up => {
                i = i - 1;
                aln1.push(s1.pop().unwrap());
                aln2.push(gap_char);
            }
            Direction::Left => {
                j = j - 1;
                aln1.push(gap_char);
                aln2.push(s2.pop().unwrap());
            }
            Direction::Diag => {
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
    let score = score_matrix[s1.len()][s2.len()];
    let (aln1, aln2) = traceback(traceback_matrix, s1, s2, gap_char);
    (score, aln1, aln2)
}
