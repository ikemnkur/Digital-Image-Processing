/*
 * hw01.c:
 *
 * - Read in two size x size BYTE images I1 and I2.
 * - Make image J so left half is equal to left half of I1 and right half
 * is equal to right half of I2. Write J to disk.
 * - Make image K so it is equal to J with right and left halves swapped.
 * Write K to disk.
 *
 * - The input images must be square.
 *
 * 02/04/2002 jph
 *
 * 02/05/2003: histograms of J and K are no longer required. jph.
 * 01/15/2013: this code is for problem 1 only. jph.
 *
 */
#include <math.h>

#include <stdio.h>

#include <stdlib.h>

#include <fcntl.h>

#include <string.h>

#define BYTE unsigned char
/*
 * Function Prototypes (forward declarations)
 */
void disk2byte();
void byte2disk();
/*----------------------------------------------------------------------*/
/* MAIN */
/*----------------------------------------------------------------------*/
main(argc, argv)
int argc;
char * argv[];

int size; /* num rows/cols in images */
int so2; /* size / 2 */
int i; /* counter */
int row; /* image row counter */
int col; /* image col counter */
BYTE * I1; /* input image I1 */
BYTE * I2; /* input image I2 */
BYTE * J; /* output image J */
BYTE * K; /* output image K */
char * InFn1; /* input filename for image I1 */
char * InFn2; /* input filename for image I2 */
char * OutFnJ; /* output filename for image J */
char * OutFnK; /* output filename for image K */
/*
 * Check for proper invocation, parse args
 */
if (argc != 6) {
  printf("\n%s: Swap image halves for hw01.", argv[0]);
  printf("\nUsage: %s size InFn1 InFn2 OutFnJ OutFnK\n",
    argv[0]);
  exit(0);
}

size = atoi(argv[1]);
if (size % 2) {
  printf("\n%s: size must be divisible by 2.\n", argv[0]);
  exit(0);
}
InFn1 = argv[2];
InFn2 = argv[3];
OutFnJ = argv[4];
OutFnK = argv[5];
so2 = size >> 1;
/*
 * Allocate image arrays
 */
if ((I1 = (BYTE * ) malloc(size * size * sizeof(BYTE))) == NULL) {
  printf("\n%s: free store exhausted.\n", argv[0]);
  exit(-1);
}
if ((I2 = (BYTE * ) malloc(size * size * sizeof(BYTE))) == NULL) {
  printf("\n%s: free store exhausted.\n", argv[0]);
  exit(-1);
}
if ((J = (BYTE * ) malloc(size * size * sizeof(BYTE))) == NULL) {
  printf("\n%s: free store exhausted.\n", argv[0]);
  exit(-1);
}
if ((K = (BYTE * ) malloc(size * size * sizeof(BYTE))) == NULL) {
  printf("\n%s: free store exhausted.\n", argv[0]);
  exit(-1);
}
/*
 * Read input images
 */
disk2byte(I1, size, size, InFn1);
disk2byte(I2, size, size, InFn2);
/*
 * Make output image J: left half is I1, right half is I2
 */
for (i = row = 0; row < size; row++) {
  for (col = 0; col < so2; col++, i++) {
    J[i] = I1[i];
  }
  for (; col < size; col++, i++) {
    J[i] = I2[i];
  }
}
/*
 * Make output image K: swap left and right halves of J
 */
for (row = 0; row < size; row++) {
  for (col = 0; col < so2; col++) {
    K[row * size + col] = J[row * size + col + so2];
  }
  for (; col < size; col++) {
    K[row * size + col] = J[row * size + col - so2];
  }
}
/*
 * Write the output images
 */
byte2disk(J, size, size, OutFnJ);
byte2disk(K, size, size, OutFnK);