#include <stdio.h>
#include <math.h>
int main(void) {
  volatile float xs[] = {-30.5f, -4.2f, -0.7f, 0.0f, 0.3f, 2.5f, 9.75f, 30.0f};
  for (int i = 0; i < 8; ++i) {
    float x = xs[i];
    printf("%d expf=%.9g exp2f=%.9g logf=%.9g log2f=%.9g exp=%.17g log=%.17g exp2=%.17g\n", i,
           expf(x), exp2f(x), logf(fabsf(x) + 1e-3f), log2f(fabsf(x) + 1e-3f),
           exp((double)x), log(fabs((double)x) + 1e-3), exp2((double)x));
  }
  return 0;
}
