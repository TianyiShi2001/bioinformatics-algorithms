/**
 * merge two sorted lists in O(a + b) time
 * @param  {Array} arr1 Array 1
 * @param  {Array} arr2 Array 2
 * @return {Array} sorted array
 */
function merge(arr1, arr2) {
  let res = [];
  let i = arr1.length - 1;
  let j = arr2.length - 1;
  while (i >= 0 && j >= 0) {
    if (arr1[i] > arr2[j]) {
      res.push(arr1.pop());
      i--;
    } else {
      res.push(arr2.pop());
      j--;
    }
  }
  while ((a = arr1.pop())) {
    res.push(a);
  }
  while ((b = arr2.pop())) {
    res.push(b);
  }
  return res.reverse();
}

/**
 * Sort an array
 * @param {Array} arr the array to be sorted
 * @return {Array} sorted array
 */
function mergeSort(arr) {
  let n = arr.length;
  if (n === 1) {
    return arr;
  }
  let mid = ~~(n / 2);
  let left = arr.slice(0, mid);
  let right = arr.slice(mid);
  return merge(mergeSort(left), mergeSort(right));
}

console.log(mergeSort([5, 1, 3, 5, 7, 6, 5, 43, 3, 5, 7, 78, 8, 12]));
