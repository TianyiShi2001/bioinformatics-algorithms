// The enum that represent the three directions, which are used in the traceback matrix
use ndarray::prelude::*;
//use num_traits::Zero;

#[derive(Debug, Clone, Copy)]
pub enum Direction {
    Up,
    Left,
    Diag,
    None,
}
// impl Zero for Direction {
//     zero() -> Self {
//         Direction::None;
//     }
// }
pub type Score = i32;
pub type Matrix<T> = Vec<Vec<T>>;
pub type ScoreMatrix = Matrix<Score>;
pub type TracebackMatrix = Matrix<Direction>;
pub type Seq = String;
pub struct Coords(pub usize, pub usize);

pub const MATCH: Score = 2i32;
pub const MISMATCH: Score = -1i32;
pub const GAP: Score = -1i32;
