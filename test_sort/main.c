#include "cstd.h"

int main()
{
  int i, j, lambda_i;
  int minval, maxval;
  int n = 10;


  int *lambda = alloc1int(n);
  srand(time(0));
  for(i=0; i<n; i++) {
    lambda[i] = rand()%n + 100;
    if(i==0){
      minval = lambda[i];
      maxval = lambda[i];
    }else{
      minval = MIN(minval, lambda[i]);
      maxval = MAX(maxval, lambda[i]);
    }
    printf("%d %d\n", i, lambda[i]);
  }
  printf("[min,max]=[%d, %d]\n", minval, maxval);

  int d = maxval-minval+1;
  int *interval_count = alloc1int(d);
  for(i=0; i<n; i++) {
    lambda_i = lambda[i]-minval;//lambda_i shifted by minval, such that index starts to 0
    interval_count[lambda_i]++;
  }
  
  /* int s = 0; */
  /* for(i=0; i<d; i++){ */
  /*   s += interval_count[i]; */
  /*   interval_count[i] = s; */
  /* } */

  /* int *sortedLambda = alloc1int(n); */
  /* for(i=n-1; i>=0; i--){ */
  /*   lambda_i = lambda[i]-minval; */
  /*   sortedLambda[interval_count[lambda_i]-1] = lambda[i]; */
  /*   interval_count[lambda_i]--; */
  /* } */

  /* for(i=0; i<n; i++){ */
  /*   printf("%d\n", sortedLambda[i]); */
  /* } */

  int *interval_ptr = alloc1int(d);
  int s = 0;
  for(i=0; i<d; i++){
    interval_ptr[i] = s;
    s += interval_count[i];
    interval_count[i] = 0;
  }

  int *index_to_node = alloc1int(n);
  int *node_to_index = alloc1int(n);
  for(i=0; i<n; i++){
    int lambda_i = lambda[i]-minval;
    j = interval_ptr[lambda_i] + interval_count[lambda_i];//index of lambda_i in sorted array

    index_to_node[j] = i;
    node_to_index[i] = j;
    interval_count[lambda_i]++;
  }
  for(j=0; j<n; j++){
    i = index_to_node[j];
    printf("%d %d\n", i, lambda[i]);
  }


  
  free1int(lambda);
}
