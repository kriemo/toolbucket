# toolbucket

A collection of utility scripts

## ena_downloader.py

Download metadata and fastqs from the European Nucleotide Archive using `aspera` or `wget`
 
```bash
python3 ena_downloader.py -h

usage: ena_downloader.py [-h] -s STUDY [-i IDS [IDS ...]] [-d DOWNLOADER]
                         [-p PROG_PATH] [-k SSHKEY] [-o OUTPUTDIR]
                         [-m METADATAFILE]

Utility to download fastqs from European Nucleotide Archive

optional arguments:
  -h, --help            show this help message and exit
  -s STUDY, --study STUDY
                        Study accessions numbers(ERP, SRP, DRP, PRJ prefixes)
                        e.g. SRX017289 or PRJNA342320
  -i IDS [IDS ...], --ids IDS [IDS ...]
                        Specific fastqs to download, defaults to all fastqs in
                        study
  -d DOWNLOADER, --downloader DOWNLOADER
                        either 'ascp' or 'wget' (default)
  -p PROG_PATH, --prog_path PROG_PATH
                        path to downloader binary, defaults to name of
                        program, i.e. 'ascp', 'wget' '
  -k SSHKEY, --sshkey SSHKEY
                        path to aspera openssh key file (i.e something like
                        ".aspera/connect/etc/asperaweb_id_dsa.openssh")
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        output directory to place fastqs, defaults to '.'
  -m METADATAFILE, --metadatafile METADATAFILE
                        filename for log file containing metadata downloaded
                        for study. defaults to study name + "_metadata.txt" in
                        outdir
```
