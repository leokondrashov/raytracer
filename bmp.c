#include "bmp.h"
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <assert.h>

void dump_bmp(const char *file, const char *data, unsigned int h, unsigned int w) {
	assert(file);
	assert(data);

	unsigned int line_len = (3 * w + 3) & (~3); // length of line in bytes with padding
	char *img = (char *) calloc(h * line_len, sizeof(char));

	char header[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	char info[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};

	unsigned int filesize = 54 + h * line_len;
	header[2] = (char)(filesize      );
	header[3] = (char)(filesize >>  8);
	header[4] = (char)(filesize >> 16);
	header[5] = (char)(filesize >> 24);

	info[ 4] = (char)(w      );
	info[ 5] = (char)(w >>  8);
	info[ 6] = (char)(w >> 16);
	info[ 7] = (char)(w >> 24);
	info[ 8] = (char)(h      );
	info[ 9] = (char)(h >>  8);
	info[10] = (char)(h >> 16);
	info[11] = (char)(h >> 24);

	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			char val = data[j * w + i];
			img[(h - 1 - j) * line_len + 3 * i + 0] = val;
			img[(h - 1 - j) * line_len + 3 * i + 1] = val;
			img[(h - 1 - j) * line_len + 3 * i + 2] = val;
		}
	}

	int fd = open(file, O_WRONLY | O_CREAT, 0666);
	write(fd, header, 14);
	write(fd, info, 40);
	write(fd, img, h * line_len);
	close(fd);
	free(img);
}
