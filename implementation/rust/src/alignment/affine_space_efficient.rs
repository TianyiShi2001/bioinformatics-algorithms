use crate::alignment::common::*;
use std::cmp::max;

pub struct SpaceEfficientAligner<'a, F: MatchFunc> {
    x: &'a Seq,
    y: &'a Seq,
    scoring: Scoring<F>,
    pub result: AlignmentResult<'a>,
}

impl<'a, F: MatchFunc> SpaceEfficientAligner<'a, F> {
    pub fn new(x: &'a Seq, y: &'a Seq, scoring: Scoring<F>) -> Self {
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
        self.result.alignment = self.nw_recursive(
            self.x,
            self.y,
            self.x.len(),
            self.y.len(),
            self.scoring.gap_open,
            self.scoring.gap_open,
        );
        self.result.score = self
            .cost_only_nw(self.x, self.y, false, self.scoring.gap_open)
            .0[self.y.len()];
    }
    /// Recursively compute alignments of sub-sequences and concatenating them
    fn nw_recursive(
        &self,
        x: &Seq,
        y: &Seq,
        m: usize,
        n: usize,
        tb: Score,
        te: Score,
    ) -> Vec<Direction> {
        // * m = x.len(); n = y.len()
        if n == 0 {
            return vec![Direction::Del; m];
        }
        if m == 0 {
            return vec![Direction::Ins; n];
        }
        if m == 1 {
            return self.nw_onerow(x[0], y, n, tb, te);
        }
        let (imid, jmid, join_by_deletion) = self.find_mid(x, y, m, n, tb, te);
        return if join_by_deletion {
            [
                self.nw_recursive(&x[..imid - 1], &y[..jmid], imid - 1, jmid, tb, 0),
                vec![Direction::Del; 2],
                self.nw_recursive(&x[imid + 1..], &y[jmid..], m - imid - 1, n - jmid, 0, te),
            ]
            .concat()
        } else {
            [
                self.nw_recursive(
                    &x[..imid],
                    &y[..jmid],
                    imid,
                    jmid,
                    tb,
                    self.scoring.gap_open,
                ),
                self.nw_recursive(
                    &x[imid..],
                    &y[jmid..],
                    m - imid,
                    n - jmid,
                    0,
                    self.scoring.gap_open,
                ),
            ]
            .concat()
        };
    }

    fn find_mid(
        &self,
        x: &Seq,
        y: &Seq,
        m: usize,
        n: usize,
        tb: Score,
        te: Score,
    ) -> (usize, usize, bool) {
        let imid = m / 2;
        let (cc_upper, dd_upper) = self.cost_only_nw(&x[..imid], y, false, tb);
        let (cc_lower, dd_lower) = self.cost_only_nw(&x[imid..], y, true, te);
        let mut max = Score::MIN;
        let mut jmid = 0;
        let mut join_by_deletion = false;
        for j in 0..=n {
            let c = cc_upper[j] + cc_lower[n - j];
            if c > max {
                max = c;
                jmid = j;
                join_by_deletion = false;
            }
            let d = dd_upper[j] + dd_lower[n - j] - self.scoring.gap_open; // subtract duplicating open!
            if d > max {
                max = d;
                jmid = j;
                join_by_deletion = true;
            }
        }
        (imid, jmid, join_by_deletion)
    }

    fn cost_only_nw(&self, x: &Seq, y: &Seq, rev: bool, tx: Score) -> (Vec<Score>, Vec<Score>) {
        let m = x.len() + 1;
        let n = y.len() + 1;
        let mut cc: Vec<Score> = vec![0; n]; // match/mismatch
        let mut dd: Vec<Score> = vec![0; n]; // deletion
        let mut e: Score; // I(i, j-1)
        let mut c: Score; // C(i, j-1)
        let mut s: Score; // C(i-1, j-1)
        let mut t: Score;
        t = self.scoring.gap_open; // originally self.scoring.gap_open;
        for j in 1..n {
            t += self.scoring.gap_extend;
            cc[j] = t;
            dd[j] = Score::MIN;
        }
        t = tx; // originally self.scoring.gap_open;
        for i in 1..m {
            s = cc[0];
            t += self.scoring.gap_extend;
            c = t;
            cc[0] = c;
            e = Score::MIN;
            for j in 1..n {
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = max(dd[j], cc[j] + self.scoring.gap_open) + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = if rev {
                    max(
                        max(dd[j], e),
                        s + self.scoring.match_fn.score(x[m - i - 1], y[n - j - 1]),
                    )
                } else {
                    max(
                        max(dd[j], e),
                        s + self.scoring.match_fn.score(x[i - 1], y[j - 1]),
                    )
                };
                s = cc[j];
                cc[j] = c;
            }
        }
        (cc, dd)
    }
    fn nw_onerow(&self, x: u8, y: &Seq, n: usize, tb: Score, te: Score) -> Vec<Direction> {
        let score_by_indels_only =
            max(tb, te) + self.scoring.gap_extend * (n as Score + 1) + self.scoring.gap_open;
        let mut max = score_by_indels_only;
        let score_with_one_substitution_BASE =
            (n as Score - 1) * self.scoring.gap_extend + self.scoring.gap_open; // plus substitution score and possibly one more gap_open
        let mut maxj_ = 0usize;
        for j_ in 0..n {
            // index of sequence instead of matrix; y[j] instead of j[j-1] is the jth character
            let score = score_with_one_substitution_BASE
                + self.scoring.match_fn.score(x, y[j_])
                + if j_ == 0 || j_ == n - 1 {
                    0
                } else {
                    self.scoring.gap_open
                };
            if score > max {
                max = score;
                maxj_ = j_;
            }
        }
        return if max == score_by_indels_only {
            let mut res = Vec::with_capacity(n + 1);
            res.push(Direction::Del);
            for i in 0..n {
                res.push(Direction::Ins)
            }
            res
        } else {
            let mut res = Vec::with_capacity(n);
            for i in 0..maxj_ {
                res.push(Direction::Ins)
            }
            res.push(Direction::Sub);
            for i in 0..(n - maxj_ - 1) {
                res.push(Direction::Ins)
            }
            res
        };
    }
}
