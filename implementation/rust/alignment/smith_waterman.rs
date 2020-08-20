use std::env;
mod commons;
use commons::*;

fn main() {
    let query: Vec<String> = env::args().skip(1).take(2).collect();
    let (score, aln1, aln2) = smith_waterman(&query[0], &query[1], '-');
    println!("{}\n{}\n{}", score, aln1, aln2);
}

fn init_score_and_traceback_matrices(nrow: usize, ncol: usize) -> (ScoreMatrix, TracebackMatrix) {
    let score_matrix = vec![vec![0isize; ncol]; nrow];
    let traceback_matrix = vec![vec![Direction::None; ncol]; nrow];
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
    let mut max: Score = 0;
    let mut dir: Direction = Direction::None;
    if up > max {
        max = up;
        dir = Direction::Up;
    }
    if left > max {
        max = left;
        dir = Direction::Left;
    }
    if diag > max {
        max = diag;
        dir = Direction::Diag;
    }
    (max, dir)
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

fn find_max_and_coords(score_matrix: &ScoreMatrix) -> (Score, Coords) {
    let mut max: Score = 0;
    let mut x = 1usize;
    let mut y = 1usize;
    for i in 1..(score_matrix.len()) {
        for j in 1..(score_matrix[0].len()) {
            let n = score_matrix[i][j];
            if n > max {
                max = n;
                x = i;
                y = j;
            }
        }
    }
    (max, Coords(x, y))
}

fn traceback(
    traceback_matrix: TracebackMatrix,
    s1: &str,
    s2: &str,
    gap_char: char,
    coords: Coords,
) -> (Seq, Seq) {
    let Coords(mut i, mut j) = coords;
    let mut aln1: Vec<char> = Vec::with_capacity(i);
    let mut aln2: Vec<char> = Vec::with_capacity(j);
    while (i > 0) && (j > 0) {
        let dir = traceback_matrix[i][j];
        match dir {
            Direction::Up => {
                i = i - 1;
                aln1.push(s1.chars().nth(i).unwrap());
                aln2.push(gap_char);
            }
            Direction::Left => {
                j = j - 1;
                aln1.push(gap_char);
                aln2.push(s2.chars().nth(j).unwrap());
            }
            Direction::Diag => {
                i = i - 1;
                j = j - 1;
                aln1.push(s1.chars().nth(i).unwrap());
                aln2.push(s2.chars().nth(j).unwrap());
            }
            _ => break,
        }
    }
    return (
        aln1.into_iter().rev().collect(),
        aln2.into_iter().rev().collect(),
    );
}

fn smith_waterman(s1: &str, s2: &str, gap_char: char) -> (Score, Seq, Seq) {
    let (score_matrix, traceback_matrix) = compute_score_and_traceback_matrices(s1, s2);
    let (score, coords) = find_max_and_coords(&score_matrix);
    let (aln1, aln2) = traceback(traceback_matrix, s1, s2, gap_char, coords);
    (score, aln1, aln2)
}
