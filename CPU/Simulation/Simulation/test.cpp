#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double get_random() { return ((double)rand() / (double)RAND_MAX); }
int main() {
  double n = 0;
  srand(1);
  for (int i = 0; i < 100; i++) {

    n = (double)rand() / (double)RAND_MAX;
    
    printf("%f\n", n);
  }
  return 0;
}
