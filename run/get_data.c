#include <stdio.h>
#include <stdlib.h>

// function to load data from psim5 bin file
void get_data(char *filename, int rows, int cols, float *data) {
    FILE *fh = NULL;
    float *buf = NULL;

    if ((rows == 0) || (cols == 0)) {
        fprintf(stderr, "Error: rows and cols must be nonzero\n");
        exit(1);
    }

    // open file
    fh = fopen(filename, "rb");
    if (fh == NULL) {
        fprintf(stderr, "Error opening input file: %s\n", filename);
        exit(1);
    }

    // allocate column buffer
    buf = malloc(rows * sizeof(float));
    if (buf == NULL) {
        fprintf(stderr, "Error allocating buf\n");
        exit(1);
    }

    // load data
    for (int c = 0; c < cols; c++) {
        float colmin, colmax;
        int nsuc = fread(buf, sizeof(float), rows, fh);
        if (nsuc != rows) {
            fprintf(stderr, "Error reading column %d\n", c);
            exit(1);
        }

        // find min/max for this column
        colmin = buf[0];
        colmax = buf[0];
        for (int r = 1; r < rows; r++) {
            colmin = (buf[r] < colmin) ? buf[r] : colmin;
            colmax = (buf[r] > colmax) ? buf[r] : colmax;
        }

        // store result
        data[c] = colmin;
        data[cols + c] = colmax;
    }

    free(buf);
    fclose(fh);
}
