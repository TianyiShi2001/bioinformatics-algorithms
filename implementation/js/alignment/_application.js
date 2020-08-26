import { show_alignment } from "./_show.js";
import { needleman_wunsch } from "./needleman_wunsch.js";
import { smith_waterman } from "./smith_waterman.js";
import { needleman_wunsch_semi } from "./needleman_wunsch_semi.js";
import { BLOSUM45 } from "./matrices/matrices.js";

// main
function main() {
  let [alg, s1, s2] = process.argv.slice(2);
  let s, a1, a2;
  if (alg === "nw") {
    [s, a1, a2] = needleman_wunsch(s1, s2, BLOSUM45, -5);
  } else if (alg === "sw") {
    [s, a1, a2] = smith_waterman(s1, s2, BLOSUM45, -5);
  } else if (alg === "nws") {
    [s, a1, a2] = needleman_wunsch_semi(s1, s2, BLOSUM45, -5);
  }
  console.log(show_alignment(a1, a2));
  console.log(s, a1, a2);
}

main();
