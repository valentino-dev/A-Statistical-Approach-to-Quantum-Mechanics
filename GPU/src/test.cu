#include <stdio.h>
#include <time.h>

int main() {
  int log_length = 20;
  int log_size = log_length * sizeof(int);
  int *log;
  log = (int *)malloc(log_size);
  for (int i = 0; i < log_length; i++) {
    log[i] = (int)'_';
  }
  // int message_length = 10;
  char message[10];
  // for (int i = 0; i < 14; i++) {
  //  message[i] = '_';
  //}
  sprintf(message, "%2.6lf|", 420.32);


  printf("%s\n", message);
  for (int i = 0; i < sizeof(message); i++) {
    printf("|%c|\n", message[i]);
    log[i] = (int)message[i];
  }

  for (int i = 0; i < log_length; i++) {
    printf("%c", log[i]);
  }
  printf("|");
  printf("\n%lo", sizeof(message));
  free(log);
  return 0;
}
