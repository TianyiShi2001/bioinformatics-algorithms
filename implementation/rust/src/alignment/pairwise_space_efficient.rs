use crate::alignment::common::*;
use std::cmp::max;
use std::mem;

pub struct SpaceEfficientAligner<'a, F: MatchFunc> {
    x: &'a Seq,
    y: &'a Seq,
    scoring: Scoring<F>,
    pub result: AlignmentResult<'a>,
}

impl<'a, F: MatchFunc> SpaceEfficientAligner<'a, F> {
    pub fn new(x: &'a Seq, y: &'a Seq, scoring: Scoring<F>) -> Self {
        let m = x.len() + 1;
        let n = y.len() + 1;
        SpaceEfficientAligner {
            x,
            y,
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
        self.result.alignment = self.nw_recursive(0, self.x.len(), 0, self.y.len());
        self.result.score = self.cost_only_nw(self.x, self.y, false)[self.y.len()];
    }
    fn nw_recursive(&self, i1: usize, i2: usize, j1: usize, j2: usize) -> Vec<Direction> {
        if i1 + 1 == i2 {
            return self.nw_onerow(self.x[i1], &self.y[j1..j2]);
        } else {
            let imid = (i1 + i2 + 1) / 2;
            let jmid = self.find_mid(&self.x[i1..i2], &self.y[j1..j2]) + j1;
            return [
                self.nw_recursive(i1, imid, j1, jmid),
                self.nw_recursive(imid, i2, jmid, j2),
            ]
            .concat();
        }
    }
    fn find_mid(&self, x: &Seq, y: &Seq) -> usize {
        let m = x.len() + 1;
        let imid = m / 2;
        let upper = self.cost_only_nw(&x[..imid], y, false);
        let mut lower = self.cost_only_nw(&x[imid..], y, true);
        lower.reverse();
        let mut max = Score::MIN;
        let mut jmid = 0;
        for j in 0..=y.len() {
            let x = upper[j] + lower[j];
            if x > max {
                max = x;
                jmid = j;
            }
        }
        jmid
    }
    fn cost_only_nw(&self, x: &Seq, y: &Seq, rev: bool) -> Vec<Score> {
        let m = x.len() + 1;
        let n = y.len() + 1;
        let mut prev: Vec<Score> = vec![0; n];
        let mut curr: Vec<Score> = vec![0; n];
        for j in 1..n {
            curr[j] = curr[j - 1] + self.scoring.gap_extend; // 0th row
        }
        for i in 1..m {
            mem::swap(&mut prev, &mut curr);
            curr[0] = prev[0] + self.scoring.gap_extend;
            for j in 1..n {
                let up = prev[j] + self.scoring.gap_extend;
                let left = curr[j - 1] + self.scoring.gap_extend;
                let diag = if rev {
                    prev[j - 1] + self.scoring.match_fn.score(x[m - i - 1], y[n - j - 1])
                } else {
                    prev[j - 1] + self.scoring.match_fn.score(x[i - 1], y[j - 1])
                };
                curr[j] = max(up, max(left, diag));
            }
        }
        curr
    }
    fn nw_onerow(&self, x: u8, y: &Seq) -> Vec<Direction> {
        let mut S = Vec::<Score>::with_capacity(y.len() + 1);
        let mut T = Vec::<Direction>::with_capacity(y.len() + 1);
        S.push(self.scoring.gap_extend);
        for j in 0..y.len() {
            let up = (j as Score + 2) * self.scoring.gap_extend;
            let left = S[j] + self.scoring.gap_extend;
            let diag =
                (j as Score) * self.scoring.gap_extend + self.scoring.match_fn.score(x, y[j]);
            let (max, dir) = compute_max_score_and_direction(up, left, diag);
            S.push(max);
            T.push(dir);
        }
        match T.iter().rposition(|&x| x == Direction::Sub) {
            Some(n) => {
                for j in 0..y.len() {
                    T[j] = Direction::Ins;
                }
                T[n] = Direction::Sub
            }
            None => {
                for j in 0..y.len() {
                    T[j] = Direction::Ins
                }
                T.push(Direction::Del);
            }
        }
        T
    }
}
