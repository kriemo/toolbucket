extern "C" {
#include "kseq.h"
}

#include <stdio.h>
#include <zlib.h>
#include <stdio.h>
#include <ostream>
#include <iostream>

void usage(){
  std::cerr << "checkfq read1.fq.gz read2.fq.gz"  
              << std::endl;
}

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
  gzFile fp, fp2;
  kseq_t *seq_1;
  kseq_t *seq_2;
  int l ;
  
  if (argc != 3) {
    usage() ;
    return 1 ;
  }

  fp = gzopen(argv[1], "r"); 
  fp2 = gzopen(argv[2], "r");
  
  seq_1 = kseq_init(fp); 
  seq_2 = kseq_init(fp2);
  
  int read_count = 0;
  while ((l = kseq_read(seq_1)) >= 0) {
    
      std::string fq1_name = seq_1->name.s ;

    int k = kseq_read(seq_2) ;
    if (k < 0) {
      std::cerr << "problem parsing fastq 2 at read " 
                << read_count 
                << " : " 
                << fq1_name 
                << std::endl;
      return 1;
    }

    std::string fq2_name = seq_2->name.s ;
    if (fq1_name != fq2_name) {
      std::cerr << "name mismatch between read_1 " << fq1_name
                << " and read 2 " << fq2_name  
                << std::endl 
                << "at read " 
                << read_count 
                << std::endl ;
      return 1 ;
    }
    
    read_count += 1 ;
  }

  // check if any records remain in seq_2
  int k = kseq_read(seq_2) ;
  if ( k >= 0 ) {
    std::cerr << "fastq 2 has more reads than fastq 1" << std::endl;
    return 1;
  }

  kseq_destroy(seq_1); 
  kseq_destroy(seq_2); 
  gzclose(fp); 
  gzclose(fp2); 

  return 0;
}

