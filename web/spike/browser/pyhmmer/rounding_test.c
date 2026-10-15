#include <stdio.h>
#include <math.h>
int main(void) {
  volatile double xs[] = {-2.5, -1.5, -0.5, -0.49999, 0.5, 0.50001, 1.5, 2.5, 3.5, 7.2, -7.2, 1e9 + 0.5};
  for (int i = 0; i < 12; ++i) {
    double x = xs[i]; float xf = (float)x;
    printf("%d lrint=%ld lround=%ld round=%.1f roundf=%.1f floor=%.1f ceil=%.1f truncf=%.1f lrintf=%ld\n",
           i, lrint(x), lround(x), round(x), roundf(xf), floor(x), ceil(x), truncf(xf), lrintf(xf));
  }
  return 0;
}
