Array.prototype.clone = function () {
  return this.slice(0);
};

const MATCH = 2;
const MISMATCH = -1;
const GAP = -1;
const GAP_CHAR = "-";
const UP = 1;
const LEFT = 2;
const DIAG = 3;

const SUBSTITUTION_MATRIX = {
  A: { A: MATCH, C: MISMATCH, G: MISMATCH, T: MISMATCH },
  C: { A: MISMATCH, C: MATCH, G: MISMATCH, T: MISMATCH },
  G: { A: MISMATCH, C: MISMATCH, G: MATCH, T: MISMATCH },
  T: { A: MISMATCH, C: MISMATCH, G: MISMATCH, T: MATCH },
}; 

function substitution(a, b) {
  return SUBSTITUTION_MATRIX[a][b];
}

/**
 * Find the midpoint `(n/2, j); 1 <= j <= m`, which maximizes
 * `AS(S[1..n/2], T[1..j]) + AS(S[n/2+1..n], T[j+i..m])`
 * # Steps
 * 1. Calculate `AS(S[1..n/2], T[1..j]) for 1 <= j <= m`
 * 2. Calculate `AS(S[n/2+1..n], T[j+1..m]) for 1 <= j <= m`
 * 3. Determine the mid point
 * @param {string} s sequence S[1..n]
 * @param {string} t sequence T[1..m]
 */
function findMid(s, t) {
  let n = s.length + 1;
  let mid = ~~(n / 2);
  let upper = costOnlyNw(s.slice(0, mid), t);
  let lower = costOnlyNw(Array.from(s.slice(mid)).reverse(), Array.from(t).reverse()).reverse();
  let max = upper[0] + lower[0];
  let jmid = 0;
  for (let j = 1; j < upper.length; j++) {
    const sum = upper[j] + lower[j];
    if (sum > max) {
      max = sum;
      jmid = j;
    }
  }
  return jmid;
}

/**
 * Compute the optimal alignment score between S[1..n] and T[1..m] using O(n) space
 * @param {string} s sequence S
 * @param {string} t sequence T
 */
function costOnlyNw(s, t) {
  let n = s.length + 1;
  let m = t.length + 1;
  let prev = [0];
  for (let j = 1; j < m; j++) {
    prev[j] = prev[j - 1] + GAP;
  }
  let curr = Array.from({ length: m });
  for (let i = 1; i < n; i++) {
    curr[0] = prev[0] + GAP;
    for (let j = 1; j < m; j++) {
      let up = prev[j] + GAP;
      let left = curr[j - 1] + GAP;
      let diag = prev[j - 1] + substitution(s[i - 1], t[j - 1]);
      curr[j] = Math.max(up, left, diag);
    }
    prev = [...curr];
  }
  return curr;
}

/**
 *
 * @param {string} s sequence S with a single item
 * @param {string} t sequence T
 */
function nw_one_row(s, t) {
  let a = Array.from(t).reduce(
    (m, curr, idx) => {
      let up = GAP * (idx + 2);
      let left = m[0][idx] + GAP;
      let diag = GAP * idx + substitution(s, curr);
      let max = Math.max(up, left, diag);
      m[0].push(max);
      m[1].push([up, left, diag].indexOf(max) + 1);
      return m;
    },
    [[GAP], []]
  )[1];
  let diagIdx = a.reverse().indexOf(3);
  if (diagIdx === -1) {
    return [UP].concat(Array.from({ length: t.length }).fill(LEFT));
  } else {
    return Array.from({ length: t.length - diagIdx - 1 })
      .fill(LEFT)
      .concat([DIAG])
      .concat(Array.from({ length: diagIdx }).fill(LEFT));
  }
}

function space_efficient_nw_recursive(s, t, i1, i2, j1, j2) {
  if (i1 === i2 - 1) {
    return nw_one_row(s[i1], t.slice(j1, j2)); // compute the alignment using nw
  }
  let imid = ~~((i1 + i2 + 1) / 2);
  let jmid = findMid(s.slice(i1, i2), t.slice(j1, j2)) + j1;
  return space_efficient_nw_recursive(s, t, i1, imid, j1, jmid).concat(space_efficient_nw_recursive(s, t, imid, i2, jmid, j2));
}

function space_efficient_nw(s, t) {
  return space_efficient_nw_recursive(s, t, 0, s.length, 0, t.length);
}

//console.log(findMid("ACAATCC", ["A", "G", "C", "A", "T", "G", "C"]));
//console.log(space_efficient_nw("ACAATCC", ["A", "G", "C", "A", "T", "G", "C"], 0, 7, 0, 7));
//console.log(nw_one_row("C", "GC"));
let [s1, s2] = process.argv.slice(2);
//console.log(space_efficient_nw("ACAATCC", ["A", "G", "C", "A", "T", "G", "C"], 0, 7, 0, 7));
console.log(space_efficient_nw(s1, s2));
