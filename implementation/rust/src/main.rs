pub mod alignment;

fn main() {
    let res = alignment::needleman_wunsch_1::needleman_wunsch("ACAATCCGGTACTAGCTAGCTAGCGCGTTGCGTAAAGTCCGTTAGCCAAATGGAAATTTCCCGGGCGCGTATAATGCATGCATGCATGCATGCATTTTTTGGGGGGCCCCCCCTAGCATCGACTACGATCAGCTACGACTCGATCAGCTACGACTAGCATCGACTACGCTACGATCCGAGTCAGCTAGCTAGCATCGATCGACTACGACTCAGATCGATCGATCACGTCAGCTACTACAGTCAGCTACGCTAGCTACGCTAGCTACGTCAGCTACGTATAGATCAGCTCCGATCAGCGATCAGCTCAGCTACGACTGCTACGATCGACTACGTCGACTAGCTACGCTAGCTACGCTACGCTACGTCAGCTACGACTCTTTTAGCGGCTATGCATTCACGGCGGATGGCATGCGATGCGATTCGGGCATTGCATTCTAGGCTAGCGATCGCATCTGGCGTACTACTCACGAGCTACGGCGATCACGACATGCATCTGCATCGACTCCACTCGACTACGACAGTATGCATCAGCTCTCTCACGACTGCAT", "AGCATGCTGCCTGTAATGCATCAAACGCTAGCTAGCATCGATCAAATTTCCCGGGCGCGTATATAGCTAGCATCGATAGCATCAGCTTATTTGGCGGGCCACCCGCATCGACTACGATCGCTACGACTCGATCAGCTCAGCATCGACTACGCATCGTACGACTGCTCATCGACATCTCTCAGCTTAGCATCGAGACTACGACTGATCGATCGAGTCAGCTACTACGCTAGCTACGTCACAGTCAGCTACGCTAGCTACGCTAGCTACGTCACTACGGATCAGCCGATCGCGATCGCATGATCAGCTCAGCTACGACTGCTACCGACTACGTCGACTAGCTACGCTAGCTACGCGCTACGTCAGCTACGTAGCGGCTATGCATTCAGGATGGCATGCGATGCGATTCGGGCATTTCTAGGCTAGCGACGATGCTCGAGCTATTCGACTCGCATCTGACGCGTACTACTCACGAGCTACGGCGATCAGCTACGACTGCATCGACTCGCATCACTCGACTACGACTCGCTCTCTCAGCGATCGATGTCAGT", '-');
    println!("{:?}", res);
}
