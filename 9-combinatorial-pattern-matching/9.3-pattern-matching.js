/**
 * Brutal force exact pattern matching taking O(m) to O(nm) time
 * @param {string} p pattern
 * @param {string} t text
 */
function patternMatching(p, t) {
    let n = p.length;
    let m = t.length;
    let res = []
    for (let i = 0; i <= (m-n); i++) {
        let s = t.slice(i, i+n)
        if (s === p) {
            res.push(i)
        }
    }
    return res
}

console.log(patternMatching("GGT", "ATGGTCGGT"))