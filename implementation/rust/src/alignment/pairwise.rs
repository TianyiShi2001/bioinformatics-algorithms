use std::env;

fn main() {
    let mut args = env::args();
    let s1: Vec<u8> = (&mut args).skip(1).next().unwrap().into_bytes();
    let s2: Vec<u8> = args.next().unwrap().into_bytes();
    let aligner = Aligner::new(&s1, &s2, Scoring::from_scores(0, -1, 2, -1));
    let result = aligner.local();
    let (x, y) = result.as_strings('-');
    println!("x: {}\ny: {}\nscore: {}", x, y, result.score);
    //println!("{:?}\n{}", alignment, score);
}

pub struct Aligner<'a, F: MatchFunc> {
    x: &'a Seq,
    y: &'a Seq,
    n: usize,
    m: usize,
    scoring: Scoring<F>,
    S: ScoreMatrix,
    T: TracebackMatrix,
    result: AlignmentResult<'a>,
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
    pub fn global(mut self) -> AlignmentResult<'a> {
        self.init_matrices();
        self.fill_matrices(false);
        self.traceback((self.x.len(), self.y.len()));
        self.result.score = self.S[self.n - 1][self.m - 1];
        self.result
    }
    pub fn local(mut self) -> AlignmentResult<'a> {
        self.fill_matrices(true);
        let (score, coords) = self.find_max_and_coords();
        self.result.score = score;
        self.traceback(coords);
        self.result
    }
    pub fn semiglobal(mut self) -> AlignmentResult<'a> {
        self.fill_matrices(false);
        let (score, coords) = self.find_max_and_coords_in_last_row_or_column();
        self.result.score = score;
        self.traceback(coords);
        self.result
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
                let (mut max, mut dir) = Self::compute_max_score_and_direction(up, left, diag);
                if local && (max < 0) {
                    max = 0;
                    dir = Direction::None;
                }
                self.S[i][j] = max;
                self.T[i][j] = dir;
            }
        }
    }
    fn compute_max_score_and_direction(up: Score, left: Score, diag: Score) -> (Score, Direction) {
        let mut max = up;
        let mut dir = Direction::Del;
        if left > max {
            max = left;
            dir = Direction::Ins;
        }
        if diag > max {
            max = diag;
            dir = Direction::Sub;
        }
        (max, dir)
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

#[derive(Debug)]
pub struct AlignmentResult<'a> {
    pub alignment: Vec<Direction>,
    pub score: i32,
    pub x: &'a [u8],
    pub y: &'a [u8],
    i: usize,
    j: usize,
}

impl<'a> AlignmentResult<'a> {
    pub fn as_strings(&self, gap_char: char) -> (String, String) {
        let gap_char = gap_char as u8;
        let mut x: Vec<u8> = Vec::with_capacity(self.alignment.len());
        let mut y: Vec<u8> = Vec::with_capacity(self.alignment.len());
        let mut i = self.i;
        let mut j = self.j;
        for dir in &self.alignment {
            match dir {
                Direction::Del => {
                    x.push(self.x[i]);
                    y.push(gap_char);
                    i += 1;
                }
                Direction::Ins => {
                    x.push(gap_char);
                    y.push(self.y[j]);
                    j += 1;
                }
                Direction::Sub => {
                    x.push(self.x[i]);
                    y.push(self.y[j]);
                    i += 1;
                    j += 1;
                }
                Direction::None => {}
            }
        }
        unsafe {
            let x = String::from_utf8_unchecked(x);
            let y = String::from_utf8_unchecked(y);
            (x, y)
        }
    }
}

/// Trait required to instantiate a Scoring instance
pub trait MatchFunc {
    fn score(&self, a: u8, b: u8) -> i32;
}

/// A concrete data structure which implements trait MatchFunc with constant
/// match and mismatch scores
#[derive(Debug, Clone)]
pub struct MatchParams {
    pub match_score: i32,
    pub mismatch_score: i32,
}

impl MatchParams {
    /// Create new MatchParams instance with given match and mismatch scores
    ///
    /// # Arguments
    ///
    /// * `match_score` - the score for a match (should not be negative)
    /// * `mismatch_score` - the score for a mismatch (should not be positive)
    pub fn new(match_score: i32, mismatch_score: i32) -> Self {
        assert!(match_score >= 0, "match_score can't be negative");
        assert!(mismatch_score <= 0, "mismatch_score can't be positive");
        MatchParams {
            match_score,
            mismatch_score,
        }
    }
}

impl MatchFunc for MatchParams {
    #[inline]
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// The trait Matchfunc is also implemented for Fn(u8, u8) -> i32 so that Scoring
/// can be instantiated using closures and custom user defined functions
impl<F> MatchFunc for F
where
    F: Fn(u8, u8) -> i32,
{
    fn score(&self, a: u8, b: u8) -> i32 {
        (self)(a, b)
    }
}

/// Details of scoring are encapsulated in this structure.
/// An affine gap score model is used so that the gap score for a length 'k' is:
/// GapScore(k) = gap_open + gap_extend * k
#[derive(Debug, Clone)]
pub struct Scoring<F: MatchFunc> {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub match_fn: F,
    pub match_scores: Option<(i32, i32)>,
}

impl Scoring<MatchParams> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// match and mismatch scores.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_score` - the score for a match
    /// * `mismatch_score` - the score for a mismatch
    pub fn from_scores(
        gap_open: i32,
        gap_extend: i32,
        match_score: i32,
        mismatch_score: i32,
    ) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn: MatchParams::new(match_score, mismatch_score),
            match_scores: Some((match_score, mismatch_score)),
        }
    }
}

impl<F: MatchFunc> Scoring<F> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// and the score function. The clip penalties are set to MIN_SCORE by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions (also see bio::scores)
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn,
            match_scores: None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Direction {
    Del, // up
    Ins, // left
    Sub, // diagonal
    None,
}

pub type Score = i32;
pub type Matrix<T> = Vec<Vec<T>>;
pub type ScoreMatrix = Matrix<Score>;
pub type TracebackMatrix = Matrix<Direction>;
pub type Seq = Vec<u8>;
pub type Coords = (usize, usize);
