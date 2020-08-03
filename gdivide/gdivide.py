import argparse
import pysam
import gzip as gz
import os
import sys
import binascii
import math

def get_filename_extension(file_path, ndots = 1):
    file_basename = os.path.basename(file_path)
    ext = file_basename.split('.')[-ndots:]
    return '.' + '.'.join(ext)

def get_filename_without_extension(file_path):
    file_basename = os.path.basename(file_path)
    filename_without_extension = file_basename.split('.')[0]
    return filename_without_extension

def is_gz_file(filepath):
    """
    "https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed"
    """
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

def gzopen(fn, *args):
    
    if(is_gz_file(fn)):
      f = gz.open(fn, *args)
    else:
      f = open(fn, *args)
    
    return f

def is_bam(fn):
    try:
        b = pysam.AlignmentFile(fn)
        is_valid_bam = True
        b.close()
    except:
        is_valid_bam = False

    return is_valid_bam

def guess_text_file_type(fn):
  ftype = ''
  with gzopen(fn, 'rt') as f:
    
    try:
      hdr = f.readline()
      line2 = f.readline()
      fields = hdr.split("\t")  
    except:
      sys.exit('Unknown file type, please specify file type')

    if hdr.startswith('>'):
      ftype = 'FASTA'
    
    elif hdr.startswith('@'):
      if line2.startswith('@'):
        sys.exit('Unknown file type, please specify file type')
      elif any([line2.startswith(x) for x in ['A','T','C', 'G', 'N']]):
        ftype = 'FASTQ'
      else:
        sys.exit('Unknown file type, please specify file type')
    
    elif len(fields) >= 3:
      try:
        start = int(fields[1])
        end = int(fields[2])
        ftype = 'BED'
      except ValueError:
        sys.exit('Unknown file type, please specify file type')
    else:
      sys.exit('Unknown file type, please specify file type')
    
  return ftype

def guess_file(fn):
    'returns file type'
    
    if is_bam(fn):
      return 'BAM'
    else:
      return guess_text_file_type(fn)

def fastx_split(fn, 
           is_pe = False,
           nrecord = None, 
           nfile = None, 
           prefix = None,
           directory = '.'):
  
  fns = fn.split(',')
  if is_pe:
    suffix = ["R1.fastq", "R2.fastq"]
    fq1 =  pysam.FastxFile(fns[0])
    fq2 =  pysam.FastxFile(fns[1])
  else:
    fq1 = pysam.FastxFile(fns[0])
    
    fq1_fname = fq1.filename
    rec1 = next(fq1)
    if rec1.quality is None and rec1.comment is None:
        suffix = ["fasta"]
    else:
        suffix = ["fastq"]
    fq1.close()
    fq1 = pysam.FastxFile(fq1_fname)
  
  if is_gz_file(fns[0]):
    is_compressed = True
    suffix = [x + ".gz" for x in suffix]
  else:
    is_compressed = False
  
  if nfile:
    print("determining # of reads in fastx")
    
    n = 0
    for rec in fq1:
        n += 1
    # no seek method so close and reopen to regenerate iterator
    fq1_fname = fq1.filename
    fq1.close()
    fq1 = pysam.FastxFile(fq1_fname)
    nrecord = math.ceil(n / float(nfile))
    print("there are %d reads in fastx" % n)
    print("splitting fastxs into files with %d reads" % nrecord)
  
  
  record_count = 0
  fout_fns = []
  n = 1
  fout_str = os.path.join(directory, prefix + '.{}.' + '{}')
  for suf in suffix:
    split = format(n, '03')
    fout_fn = fout_str.format(split, suf)
    fout_fns.append(fout_fn)

  
  fouts = []
  for fout_fn in fout_fns:
    if is_compressed:
      fout = gz.open(fout_fn, mode = 'wt', compresslevel = 6)
    else:
      fout = open(fout_fn, mode='w')
    fouts.append(fout)
  
  print("splitting into file #: %d" % n)

  for rec in fq1:
    record_count += 1
    if record_count > nrecord:
        record_count = 1
        n += 1
        print("splitting into file #: %d" % n)
        # close previous files and make new files
        for idx,fout in enumerate(fouts):
            fout.close()
            split = format(n, '03')
            fout_fn = fout_str.format(split, suffix[idx])
            if is_compressed:
              fout = gz.open(fout_fn, mode = 'wt', compresslevel = 6)
            else:
              fout = open(fout_fn, mode='w')
            fouts[idx] = fout

    fouts[0].write(str(rec) + '\n')
    if is_pe:
      rec2 = next(fq2)
      fouts[1].write(str(rec2) + '\n')

  for f in fouts:
    f.close()
  

def bam_split(fn, 
           by_chrom = False,
           nrecord = None, 
           nfile = None, 
           prefix = None,
           directory = '.'):
  
  bam = pysam.AlignmentFile(fn)

  if by_chrom or nfile:
    if not nfile:
      chrom_splits = list(bam.references)
    else:
      idxstats = bam.get_index_statistics()
      idxstats.sort(key=lambda x: x[3], reverse = True)
      chrom_splits = [x[0] for x in idxstats[:nfile - 1]]
      remaining_chroms = [x[0] for x in idxstats[nfile - 1:]]
      if remaining_chroms:
        if len(remaining_chroms) == 1:
            remaining_chroms = remaining_chroms[0]
        chrom_splits.append(remaining_chroms)
    print(chrom_splits) 
  else:
    chrom_splits = list(bam.references)
    print("counting # of records")

def fsplit(fn, ftype, **kwargs):

  if ftype == "FASTQ" or ftype == "PE-FASTQ":
    del kwargs["by_chrom"]
    fastx_split(fn, is_pe = ftype == "PE-FASTQ", **kwargs)
  elif ftype == "FASTA":
    del kwargs["by_chrom"]
    fastx_split(fn, is_pe = False, **kwargs)
  elif ftype == "BAM":
    pass
    #bam_split(fn, **kwargs)
  else:
    pass

def main():
  parser = argparse.ArgumentParser(description="""
  General utility for splitting genomics data into
  multiple files.""")
  parser.add_argument('-i',
                      '--input',
                      help = """input file, set to - to read from stdin,
                      if paired end fastq reads then supply in csv format, e.g.
                      r1.fastq.gz,r2.fastq.gz
                      """,
                      required = True)
  parser.add_argument('-n',
                      '--nfiles',
                      help = """desired # of files to generate, specify
                      either -n or -l""",
                      type = int,
                      required = False)
  parser.add_argument('-l',
                      '--nrecords',
                      help = """desired # of records in files to generate,
                      specify either -n or -l""",
                      type = int,
                      required = False)
  parser.add_argument('-p',
                      '--prefix',
                      help = """prefix for output, defaults to filename
                      without extension""",
                      required = False)
  parser.add_argument('-o',
                      '--outdir',
                      help = 'output directory, defaults to ./',
                      required = False,
                      default = "./")
  parser.add_argument('-c',
                      '--chrom',
                      help = """if set, split by chromosome, if -n is set,
                      then will split by chrom, starting with chrom with
                      most reads, the final file will contain remaining
                      chroms. e.g. -n 20 -c will produce chr1->chr19, then
                      group remaining reads into final file named
                      other.bam. 
                      If -l and -c are set then will split per chromosome
                      files into multiple files if more than -l records
                      are present.
                      """,
                      action = 'store_true')
  parser.add_argument('-t',
                      '--type',
                      help = """format is autodetected, however if not
                      correct, then specify file type, one of BED, BAM, FASTQ,
                      PE-FASTQ, FASTA""")

  args = parser.parse_args()
  fin = args.input.split(',')
  if len(fin) > 1:
    ftype = 'PE-FASTQ'

  else:
    if args.type:
      ftype = args.type
    else:
      ftype = guess_file(args.input)
  
  print("detected %s file" % ftype) 
  
  if not args.nfiles and not args.nrecords:
    if ftype == 'BAM' and args.chrom:
        pass
    else:
        sys.exit('specify either -n or -l')
  
  outdir = args.outdir
  if not os.path.exists(outdir):
    os.makedirs(outdir)
  
  if args.prefix:
      prefix = args.prefix
  else:
      prefix = get_filename_without_extension(fin[0])
  
  if args.chrom:
      by_chrom = True
  else:
      by_chrom = False
  
  other_args = {
          'directory' : outdir,
          'prefix' : prefix,
          'by_chrom' : by_chrom}
  
  if args.nfiles:
    if args.nrecords:
      sys.exit('specify either -n or -l, not both')
    nfiles = max(1, args.nfiles)
    fsplit(args.input, ftype, nfile = nfiles, **other_args)

  elif args.nrecords:
    if args.nfiles:
      sys.exit('specify either -n or -l, not both')
    nrecords = max(1, args.nrecords)
    fsplit(args.input, ftype, nrecord = nrecords, **other_args)

  elif ftype == 'BAM' and args.chrom:
    pass
    #fsplit(args.input, ftype, **other_args)





if __name__ == '__main__':  main()

