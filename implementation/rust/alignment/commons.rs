// The enum that represent the three directions, which are used in the traceback matrix
#[derive(Debug, Clone, Copy)]
pub enum Direction {
    Up,
    Left,
    Diag,
    None,
}
pub type Score = isize;
pub type Matrix<T> = Vec<Vec<T>>;
pub type ScoreMatrix = Matrix<Score>;
pub type TracebackMatrix = Matrix<Direction>;
pub type Seq = String;
pub struct Coords(pub usize, pub usize);

pub const MATCH: Score = 2;
pub const MISMATCH: Score = -1;
pub const GAP: Score = -1;
