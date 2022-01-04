# ramdaq: troubleshooting

## Unable to untar index file

Errors like the following:

```bash
Caused by:
  Process `untar_hisat2_idx` terminated with an error exit status (2)
```

```bash
Command error:
  
  gzip: stdin: unexpected end of file
  tar: Unexpected EOF in archive
  tar: Unexpected EOF in archive
  tar: Error is not recoverable: exiting now
```

By default, ramdaq downloads the large-size index files used by hisat2 and RSEM. In rare cases, this download process terminates incompletely and ramdaq exits with decompression fail. Rerunning ramdaq cannot remove this error.

To fix this, you have to remove the run's ```work/``` directory (e.g., with ```rm -rf work/```) and re-run ramdaq from scratch.
