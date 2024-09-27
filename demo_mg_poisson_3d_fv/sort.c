#include "cstd.h"

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(int arr[], int index[], int l, int m, int r)
{
  int i, j, k;
  int n1 = m - l + 1;
  int n2 = r - m;

  /* create temp arrays */
  int *L = alloc1int(n1);
  int *R = alloc1int(n2);
  int *L2 = alloc1int(n1);
  int *R2 = alloc1int(n2);

  /* Copy data to temp arrays L[] and R[] */
  for (i = 0; i < n1; i++){
    L[i] = arr[l + i];
    L2[i] = index[l + i];
  }
  for (j = 0; j < n2; j++){
    R[j] = arr[m + 1 + j];
    R2[j] = index[m + 1 + j];
  }
  /* Merge the temp arrays back into arr[l..r]*/
  i = 0; // Initial index of first subarray
  j = 0; // Initial index of second subarray
  k = l; // Initial index of merged subarray
  while (i < n1 && j < n2) {
    if (L[i] > R[j]) {
      arr[k] = L[i];
      index[k] = L2[i];
      i++;
    } else {
      arr[k] = R[j];
      index[k] = R2[j];
      j++;
    }
    k++;
  }

  /* Copy the remaining elements of L[], if there
     are any */
  while (i < n1) {
    arr[k] = L[i];
    index[k] = L2[i];
    i++;
    k++;
  }

  /* Copy the remaining elements of R[], if there
     are any */
  while (j < n2) {
    arr[k] = R[j];
    index[k] = R2[j];
    j++;
    k++;
  }

  free1int(L);
  free1int(R);
  free1int(L2);
  free1int(R2);
}

/* l is for left index and r is right index of the
   sub-array of arr to be sorted */
void mergeSort(int arr[], int index[], int l, int r)
{
  if (l < r) {
    // Same as (l+r)/2, but avoids overflow for
    // large l and h
    int m = l + (r - l) / 2;

    // Sort first and second halves
    mergeSort(arr, index, l, m);
    mergeSort(arr, index, m + 1, r);

    merge(arr, index, l, m, r);
  }
}
