function recursiveChange(M, c) {
  if (M === 0) {
    return 0;
  }
  let bestNumCoins = Infinity;
  for (const i of c) {
    if (M >= i) {
      numCoins = recursiveChange(M - i, c);
      if (numCoins + 1 < bestNumCoins) {
        bestNumCoins = numCoins + 1;
      }
    }
  }
  console.log(bestNumCoins);
  return bestNumCoins;
}

// console.log(recursiveChange(29, [1, 3, 7]));

function change(M, c) {
  let res = [0];
  for (let m = 1; m <= M; m++) {
    res[m] = Infinity;
    for (const i of c) {
      if (m >= i) {
        if (res[m - i] + 1 < res[m]) {
          res[m] = res[m - i] + 1;
        }
      }
    }
  }
  return res[M];
}

console.log(change(77, [1, 3, 7]));
