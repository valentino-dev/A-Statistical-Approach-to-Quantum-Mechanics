#include <stdio.h>
#include <time.h>

int main ()
{
	clock_t begin = clock();
	for (int i = 0; i < (1<<20); i++) printf("%i\n", i);

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("%6.3f", time_spent);
	return 0;
}
