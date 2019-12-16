#include <zlib.h>  
#include <stdio.h>
#include <string.h>  

#include "kseq.h"  
// STEP 1: declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread)  



  
int main(int argc, char *argv[])  
{  
    gzFile fp;  
    kseq_t *seq;
    long bases   = 0;
    long q20_cnt = 0;
    long q30_cnt = 0;
    int l;  
    if (argc != 2) {  
        fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);  
        return 1;  
    }  
    fp = gzopen(argv[1], "r"); // STEP 2: open the file handler  
    seq = kseq_init(fp); // STEP 3: initialize seq  
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence 
        char *q = seq->qual.s;
        int   c = 0;
        long  f = strlen(seq->qual.s);
        while (c < f) {
            if (*q >= 53) { q20_cnt++;}
            if (*q >= 63) { q30_cnt++;}
            q++;
            c++;
        }
        bases += f;
    }
    printf("%ld\t%ld\t%ld\n", bases, q20_cnt, q30_cnt);
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp); // STEP 6: close the file handler
    return 0;
}
