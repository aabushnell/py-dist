#include <errno.h>
#include <netinet/in.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <unistd.h>

#define PORT 8080
#define BUFFER_SIZE 1000

void *handle_client(void *arg) {
  int client_fd = *((int *)arg);
  char *buffer = (char *)malloc(BUFFER_SIZE * sizeof(char));

  ssize_t bytes_recieved = recv(client_fd, buffer, BUFFER_SIZE, 0);
  printf("BUFFER: %s\n", buffer);

  close(client_fd);
  free(arg);
  free(buffer);
  return NULL;
}

int main() {
  int server_fd;
  struct sockaddr_in server_addr;

  if ((server_fd = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
    fprintf(stderr, "ERROR: socket failed %d -> %s\n", errno, strerror(errno));
    exit(EXIT_FAILURE);
  }

  server_addr.sin_family = AF_INET;
  server_addr.sin_addr.s_addr = INADDR_ANY;
  server_addr.sin_port = htons(PORT);

  if (bind(server_fd, (struct sockaddr *)&server_addr, sizeof(server_addr)) <
      0) {
    fprintf(stderr, "ERROR: bind failed %d -> %s\n", errno, strerror(errno));
    exit(EXIT_FAILURE);
  }

  if (listen(server_fd, 10) < 0) {
    fprintf(stderr, "ERROR: listen failed %d -> %s\n", errno, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while (1) {
    struct sockaddr_in client_addr;
    socklen_t client_addr_len = sizeof(client_addr);
    int *client_fd = malloc(sizeof(int));

    if ((*client_fd = accept(server_fd, (struct sockaddr *)&server_addr,
                             &client_addr_len)) < 0) {
      fprintf(stderr, "ERROR: accept failed %d -> %s\n", errno,
              strerror(errno));
      continue;
    }

    pthread_t thread_id;
    pthread_create(&thread_id, NULL, handle_client, (void *)client_fd);
    pthread_detach(thread_id);
  }

  close(server_fd);
  return EXIT_SUCCESS;
}
