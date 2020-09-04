use crate::alignment::common::*;

pub struct Aligner<'a, F: MatchFunc> {
    x: &'a Seq,
    y: &'a Seq,
    n: usize,
    m: usize,
    scoring: Scoring<F>,
    S: ScoreMatrix,
    T: TracebackMatrix,
    pub result: AlignmentResult<'a>,
}

impl<'a, F: MatchFunc> Aligner<'a, F> {
    pub fn new(x: &'a Seq, y: &'a Seq, scoring: Scoring<F>) -> Self {
        let n = x.len() + 1;
        let m = y.len() + 1;
        Aligner {
            x,
            y,
            n,
            m,
            S: vec![vec![0i32; m]; n],
            T: vec![vec![Direction::None; m]; n],
            scoring,
            result: AlignmentResult {
                alignment: Vec::new(),
                score: 0,
                x: x,
                y: y,
                i: 0usize,
                j: 0usize,
            },
        }
    }
    pub fn global(&mut self) {
        self.init_matrices();
        self.fill_matrices(false);
        self.traceback((self.x.len(), self.y.len()));
        self.result.score = self.S[self.n - 1][self.m - 1];
    }
    pub fn local(&mut self) {
        self.fill_matrices(true);
        let (score, coords) = self.find_max_and_coords();
        self.result.score = score;
        self.traceback(coords);
    }
    pub fn semiglobal(&mut self) {
        self.fill_matrices(false);
        let (score, coords) = self.find_max_and_coords_in_last_row_or_column();
        self.result.score = score;
        self.traceback(coords);
    }
    fn init_matrices(&mut self) {
        // init matrices for global alignment with terminal gap penalty
        // Semiglobal and local alignment does NOT need this initialisation
        for i in 1..self.n {
            self.S[i][0] = self.S[i - 1][0] + self.scoring.gap_extend;
            self.T[i][0] = Direction::Del;
        }
        for j in 1..self.m {
            self.S[0][j] = self.S[0][j - 1] + self.scoring.gap_extend;
            self.T[0][j] = Direction::Ins;
        }
    }
    fn fill_matrices(&mut self, local: bool) {
        for i in 1..self.n {
            for j in 1..self.m {
                let up = self.S[i - 1][j] + self.scoring.gap_extend;
                let left = self.S[i][j - 1] + self.scoring.gap_extend;
                let diag = self.S[i - 1][j - 1]
                    + self.scoring.match_fn.score(self.x[i - 1], self.y[j - 1]);
                let (mut max, mut dir) = compute_max_score_and_direction(up, left, diag);
                if local && (max < 0) {
                    max = 0;
                    dir = Direction::None;
                }
                self.S[i][j] = max;
                self.T[i][j] = dir;
            }
        }
    }
    fn traceback(&mut self, coords: Coords) {
        let (mut i, mut j) = coords;
        loop {
            let dir = self.T[i][j];
            match dir {
                Direction::Del => {
                    i -= 1;
                }
                Direction::Ins => {
                    j -= 1;
                }
                Direction::Sub => {
                    i -= 1;
                    j -= 1;
                }
                Direction::None => break,
            }
            self.result.alignment.push(dir);
        }
        self.result.alignment.reverse();
        self.result.i = i;
        self.result.j = j;
    }
    fn find_max_and_coords(&self) -> (Score, Coords) {
        // for local alignment
        let mut max = i32::MIN;
        let mut max_i = 0;
        let mut max_j = 0;
        for i in 0..self.n {
            for j in 0..self.m {
                let n = self.S[i][j];
                if n > max {
                    max = n;
                    max_i = i;
                    max_j = j;
                }
            }
        }
        return (max, (max_i, max_j));
    }
    fn find_max_and_coords_in_last_row_or_column(&self) -> (Score, Coords) {
        // for semiglobal alignment
        let mut max = i32::MIN;
        let mut i_max = 0usize;
        let mut j_max = 0usize;
        for j in 0..(self.m - 1) {
            let a = self.S[self.n - 1][j];
            if a > max {
                max = a;
                i_max = self.n - 1;
                j_max = j;
            }
        }
        for i in 0..(self.n - 1) {
            let a = self.S[i][self.m - 1];
            if a > max {
                max = a;
                i_max = i;
                j_max = self.m - 1;
            }
        }
        (max, (i_max, j_max))
    }
}
