use bioinformatics_algorithms::alignment::common::Scoring;
use bioinformatics_algorithms::alignment::common::*;
use bioinformatics_algorithms::alignment::pairwise::*;
//use bioinformatics_algorithms::alignment::pairwise_space_efficient::*;
use bioinformatics_algorithms::alignment::affine_space_efficient::SpaceEfficientAligner;
use kodama::{linkage, Method};
// use bio::alignment::pairwise::*;
use std::env;

fn main() {
    run_aln();
}

fn run_aln() {
    let mut args = env::args();
    let s1: Vec<u8> = (&mut args).skip(1).next().unwrap().into_bytes();
    let s2: Vec<u8> = args.next().unwrap().into_bytes();
    // let mut aligner = Aligner::new(&s1, &s2, Scoring::from_scores(-5, -1, 2, -1));
    let mut aligner_se = SpaceEfficientAligner::new(&s1, &s2, Scoring::from_scores(-5, -1, 2, -1));
    // aligner.global();
    //println!("x: {}\ny: {}\nscore: {}", x, y, aligner.result.score);

    aligner_se.global();
    let (x, y) = aligner_se.result.as_strings('-');
    println!("x: {}\ny: {}\nscore: {}", x, y, aligner_se.result.score);
    //println!("{:?}\n{}", alignment, score);
}

fn run_cluster() {
    let res = linkage(
        &mut [
            0.111, 0.25, 0.555, 0.444, 0.375, 0.222, 0.111, 0.5, 0.5, 0.111,
        ],
        5usize,
        Method::Average,
    );
    println!("{:?}", res);
    for step in res.steps() {
        println!("{:?}", step);
    }
}
