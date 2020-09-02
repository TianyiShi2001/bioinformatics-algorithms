function duplicateRemoval(arr) {
    let max = Math.max(arr);
    let b = Array.from({length: max}).fill(0);
    let res = []
    for (const num of arr) {
        b[num] = 1
    }
    b.map((exist, n) => exist && res.push(n))
    return res;
}

console.log(duplicateRemoval([0,10,8,10,3,4,5,4]))