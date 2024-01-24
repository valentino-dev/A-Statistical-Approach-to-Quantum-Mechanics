#include <curand.h>
#include <curand_kernel.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <cjson/cJson.h>

class Settings{
	int MCS = 0;
	int MCI = 0;
	int runs = 0;
	int xj_perThread = 0;
	int threadsPerBlock = 0;
	int n = 0;
	int StatisticalIndependentIterations = 0;

	void readJson(char* name[]){
		FILE *jasonPtr = fopen(name, "r");
			if (jasonPtr == NULL){
			printf("ERROR: Unable to open the file.\n")
		}

		char buffer[1024];
		int len = fread(buffer, 1, sizeof(buffer), jasonPtr);

		cJSON *json  = cJSON_Parse(buffer);
		if (json == NULL) { 
        	const char *error_ptr = cJSON_GetErrorPtr(); 
        	if (error_ptr != NULL) { 
            	printf("Error: %s\n", error_ptr); 
        	} 
        	cJSON_Delete(json); 
    	} 

		this.MCS = cJSON_GetObjectItemCaseSensitive(json, "MCS")->valueint
		this.MCI = cJSON_GetObjectItemCaseSensitive(json, "MCI")->valueint
		this.xj_perThread = cJSON_GetObjectItemCaseSensitive(json, "xj_perThread")->valueint
		this.threadsPerBlock = cJSON_GetObjectItemCaseSensitive(json, "threadsPerBlock")->valueint
		this.n = cJSON_GetObjectItemCaseSensitive(json, "n")->valueint
		this.runs = cJSON_GetObjectItemCaseSensitive(json, "runs")->valueint
		this.StatisticalIndependentIterations = cJSON_GetObjectItemCaseSensitive(json, "StatisticalIndependentIterations")->valueint

        cJSON_Delete(json); 
	}



}

int main(){
	FILE *file;
	file = fopen("data/data_fig6.csv", "r");

	char data[100]

	

	return 0;
}
