#include <stdio.h>

int main(){
	FILE *data_file;
	data_file = fopen("test.csv", "w");
	fprintf(data_file, "test\n");
	fclose(data_file);
	return 0;
}
