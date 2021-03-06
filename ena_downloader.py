#! /usr/bin/env python3

import os
import sys
import argparse
import urllib.request
import urllib.error
from subprocess import call

""" Utility to download fastq and metadata files from ENA database using
just the project identifier (i.e SRX.... or PRJN....)

Uses either wget (slow) or aspera (fast)
"""

def get_study_metadata(study_id, logfile):
  
  baseurl = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="
  fields_url = "&result=read_run&fields=run_accession,sample_title,fastq_ftp,fastq_aspera"
  
  full_url = baseurl + study_id + fields_url

  # get study info
  try: 
    fp = urllib.request.urlopen(full_url)
  except urllib.error.HTTPError as e:
      sys.exit("unable to retrieve study from ENA: error {}".format(e.code)) 
  except urllib.error.URLError as e:
      sys.exit("unable to retrieve study from ENA: error {}".format(e.reason)) 
  
  mdata = []
  for line in fp:
    line = line.decode("utf8")
    logfile.write(line)
    mdata.append(line)
  
  fp.close()
  
  if len(mdata) < 2:
      # ENA will return only headers with some input strings
      sys.exit("unable to retrieve study from ENA")
  
  return mdata

def download_files(metadata, fq_ids, dl_prog, dl_prog_path, ssh_key,
        output_dir):

  ascp_cmd = [
        dl_prog_path,
        "-QT",  
        "-l",
        "200m",
        "-P",
        "33001",
        "-i",
        ssh_key]

  wget_cmd = [
    dl_prog_path
  ]

  header = metadata.pop(0)
  
  if fq_ids:
    filter_fqs = True
  else:
    filter_fqs = False

  for line in metadata:
     accession, title, ftp_url, ascp_url = line.rstrip().split("\t")

     if filter_fqs:
       if accession not in fq_ids:
         continue

     ftp_urls = ftp_url.split(";")
     for i,furl in enumerate(ftp_urls):
         if furl.startswith('ftp://'):
             pass
         else:
             ftp_urls[i] =  "ftp://" + furl 
     
     ascp_urls = ascp_url.split(";")
     
     if len(ftp_urls) == 2 and len(ascp_urls) == 2:
         libtype = "paired_end"
     elif len(ftp_urls) == 1 and len(ascp_urls) == 1:
         libtype = "single_end"
     else:
         sys.exit("unknown ftp or ascp urls in mdata")
     
     dl_cmds = []
     if dl_prog == 'ascp':
       fqs = [os.path.join(output_dir, os.path.basename(x)) for x in ascp_urls]
          
       for idx,fq in enumerate(fqs):
         dl_cmd = ascp_cmd + ["era-fasp@" + ascp_urls[idx], fq]
         dl_cmds.append(dl_cmd)

     elif dl_prog == 'wget':
       fqs = [os.path.join(output_dir, os.path.basename(x)) for x in ftp_urls]

       for idx,fq in enumerate(fqs):
         dl_cmd = wget_cmd + ["-O", fq, ftp_urls[idx]]
         dl_cmds.append(dl_cmd)

     else:
       sys.exit("unknown download program {}".format(dl_prog))

     
     for cmd in dl_cmds:
       print("downloading sample {} from accession {}".format(accession, title))
       print(" ".join(cmd))
       call(cmd)



def main():
  parser = argparse.ArgumentParser(description = """
    Utility to download fastqs from European Nucleotide Archive
    """)

  parser.add_argument('-s',
                      '--study',
                      help = """
                      Study accessions numbers(ERP, SRP, DRP, PRJ prefixes)
                      e.g. SRX017289 or PRJNA342320
                      """,
                      required = True)

  parser.add_argument('-i',
                      '--ids',
                      help = """
                      Specific fastqs to download, defaults to all fastqs
                      in study
                      """,
                      required = False, nargs = "+")
  
  parser.add_argument('-d',
                      '--downloader',
                      help ="""either 'ascp' or 'wget' (default) """,
                   required = False,
                   default = "wget")
  
  parser.add_argument('-p',
                      '--prog_path',
                      help ="""path to downloader binary, defaults to name of program,
                      i.e. 'ascp', 'wget' '""",
                   required = False)
  
  parser.add_argument('-k',
                      '--sshkey',
                      help ="""path to aspera openssh key file (i.e something like
                      ".aspera/connect/etc/asperaweb_id_dsa.openssh") """,
                   required = False)
  
  parser.add_argument('-o',
                      '--outputdir',
                      help ="""output directory to place fastqs, defaults
                      to '.' """,
                      default = ".",
                   required = False)
  
  parser.add_argument('-m',
                      '--metadatafile',
                      help ="""
                      preparsed metadata file for downloading
                      """,
                   required = False)

  args=parser.parse_args()
  
  study_id = args.study
  fq_ids = args.ids
  dl_prog = args.downloader
  dl_prog_path = args.prog_path
  ssh_key = args.sshkey
  output_dir = args.outputdir 

  if dl_prog == "ascp" and ssh_key is None:
    sys.exit("-k required for using aspera")
  if dl_prog_path is None:
    dl_prog_path = dl_prog
  
  if output_dir:
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
  else:
    output_dir = "."
  
  logfile = os.path.join(output_dir, study_id + "_download_log.txt")

  log_fp = open(logfile, 'w')

  if args.metadatafile:
      mdata = []
      with open(args.metadatafile) as f:
          for line in f:
              mdata.append(line)
  else:
      mdata = get_study_metadata(study_id, log_fp)
  
  download_files(mdata,
          fq_ids,
          dl_prog, 
          dl_prog_path, 
          ssh_key,
          output_dir)

  log_fp.close()

if __name__ == '__main__': main()

