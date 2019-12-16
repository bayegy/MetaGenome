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
    long seqs    = 0;
    long bases   = 0;
    long q20_cnt = 0;
    long q30_cnt = 0;
    long gc_cnt  = 0;
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
        while (c < strlen(seq->qual.s)) {
            if (*q - 33 >= 20) { q20_cnt++;}
            if (*q - 33 >= 30) { q30_cnt++;}
            q++;
            c++;
        }

        char *s = seq->seq.s;
        int   d = 0;
        while (d < strlen(seq->seq.s)) {
            if (*s == 'C' || *s == 'G') { gc_cnt++; }
            s++;
            d++;
        }

       bases += strlen(seq->seq.s);
       seqs += 1;  
    }  
    printf("%ld\t%ld\t%ld\t%ld\t%ld\n", seqs, bases, q20_cnt, q30_cnt, gc_cnt);      
    kseq_destroy(seq); // STEP 5: destroy seq  
    gzclose(fp); // STEP 6: close the file handler  
    return 0;  
}
