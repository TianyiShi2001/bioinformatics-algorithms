use bioinformatics_algorithms::alignment::common::*;
use bioinformatics_algorithms::alignment::pairwise::*;
use bioinformatics_algorithms::alignment::pairwise_space_efficient::*;
// use bio::alignment::pairwise::*;
use std::env;

fn main() {
    let mut args = env::args();
    let s1: Vec<u8> = (&mut args).skip(1).next().unwrap().into_bytes();
    let s2: Vec<u8> = args.next().unwrap().into_bytes();
    let mut aligner = Aligner::new(&s1, &s2, Scoring::from_scores(0, -1, 2, -1));
    let mut aligner_se = SpaceEfficientAligner::new(&s1, &s2, Scoring::from_scores(0, -1, 2, -1));
    aligner.global();
    let (x, y) = aligner.result.as_strings('-');
    println!("x: {}\ny: {}\nscore: {}", x, y, aligner.result.score);

    aligner_se.global();
    let (x, y) = aligner_se.result.as_strings('-');
    println!("x: {}\ny: {}\nscore: {}", x, y, aligner_se.result.score);
    //println!("{:?}\n{}", alignment, score);
}
