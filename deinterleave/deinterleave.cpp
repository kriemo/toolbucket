#include <zlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "gzstream/gzstream.h"

extern "C" {
#include "kseq.h"
}


KSEQ_INIT(gzFile, gzread)

int write_error(const char* seq_id){
    fprintf(stderr, "Error writing read %s \n", seq_id) ;
    return 1;

}

void write_record(ogzstream& fp, kseq_t *seq){
  
    fp << "@" << seq->name.s << " " << seq->comment.s << "\n"
      << seq->seq.s << "\n"
      << "+\n"
      << seq->qual.s << std::endl ;
}

int main(int argc, char *argv[])
{
  gzFile fp;
  kseq_t *seq;
  int l ;
  
  if (argc == 1) {
    fprintf(stderr, "Usage: %s <in.seq1> <in.seq2> <cbc_start> <cbc_end> <umi_start> <umi_end> <delim> \n", argv[0]);
    return 1;
  }

  fp = gzopen(argv[1], "r"); 
  ogzstream fp1("1.fq.gz");
  ogzstream fp2("2.fq.gz");

  seq = kseq_init(fp); 

  while ((l = kseq_read(seq)) >= 0) {
    
    std::string seq_id = seq->name.s ;
    write_record(fp1, seq) ;

    if ((l = kseq_read(seq)) >= 0){
      if(seq_id == seq->name.s) {
        write_record(fp2, seq) ; 
      } else {
        fprintf(stderr, "Error names are mismatched %s and %s", seq->name.s, seq_id.c_str()); 
        return 1 ;
      }
    }
  }

  
  kseq_destroy(seq); 
  gzclose(fp); 
  fp1.close();
  fp2.close();
  return 0;
}

