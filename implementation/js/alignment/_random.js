const chars = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X"];

function choose(choices) {
  let index = ~~(Math.random() * choices.length);
  return choices[index];
}

function randomSeq(length) {
  let res = [];
  for (let i = 0; i < length; i++) {
    res.push(choose(chars));
  }
  return res.join("");
}

function mutate(string, chars = chars, del = 0.1, ins = 0.1, replace = 0.1) {
  string = Array.from(string);
  for (let i = 0; i < ~~(del * string.length); i++) {
    let j = ~~(Math.random() * string.length);
    string = string.slice(0, j).concat(string.slice(j + 1));
  }
  for (let i = 0; i < ~~(ins * string.length); i++) {
    let j = ~~(Math.random() * string.length);
    string = string
      .slice(0, j)
      .concat([choose(chars)])
      .concat(string.slice(j));
  }
  for (let i = 0; i < ~~(replace * string.length); i++) {
    let j = ~~(Math.random() * string.length);
    string = string
      .slice(0, j)
      .concat([choose(chars)])
      .concat(string.slice(j + 1));
  }
  return string.join("");
}

function main() {
  let [m, a, b, c] = process.argv.slice(2);
  let s1 = randomSeq(m);
  let s2 = mutate(s1, chars, a, b, c);
  console.log(s1 + " " + s2);
}
main();
