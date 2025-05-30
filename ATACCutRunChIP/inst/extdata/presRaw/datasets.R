cnr_paths <- c("ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/018/SRR20110418/SRR20110418_1.fastq.gz",
"ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/018/SRR20110418/SRR20110418_2.fastq.gz",
"ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/017/SRR20110417/SRR20110417_1.fastq.gz",
"ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/017/SRR20110417/SRR20110417_2.fastq.gz",
"ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/009/SRR20110409/SRR20110409_1.fastq.gz",
"ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/009/SRR20110409/SRR20110409_2.fastq.gz",
"ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/007/SRR20110407/SRR20110407_1.fastq.gz",
"ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/007/SRR20110407/SRR20110407_2.fastq.gz")

cnr_names <- c("SOX9CNR_D0_rep1_R1.fastq.gz",
  "SOX9CNR_D0_rep1_R2.fastq.gz",
  "SOX9CNR_D0_rep2_R1.fastq.gz",
  "SOX9CNR_D0_rep2_R2.fastq.gz",
  "SOX9CNR_W6_rep1_R1.fastq.gz",
  "SOX9CNR_W6_rep1_R2.fastq.gz",
  "SOX9CNR_W6_rep2_R1.fastq.gz",
  "SOX9CNR_W6_rep2_R2.fastq.gz")

atac_paths <- c("ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/022/SRR20111522/SRR20111522_1.fastq.gz",
                "ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/022/SRR20111522/SRR20111522_2.fastq.gz",
                "ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/021/SRR20111521/SRR20111521_1.fastq.gz",
                "ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/021/SRR20111521/SRR20111521_2.fastq.gz",
                "ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/016/SRR20111516/SRR20111516_1.fastq.gz",
                "ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/016/SRR20111516/SRR20111516_2.fastq.gz",
                "ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/015/SRR20111515/SRR20111515_1.fastq.gz",
                "ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/015/SRR20111515/SRR20111515_2.fastq.gz")

atac_names <- c("D0_ATAC_rep1_R1.fastq.gz",
                "D0_ATAC_rep1_R2.fastq.gz",
                "D0_ATAC_rep2_R1.fastq.gz",
                "D0_ATAC_rep2_R2.fastq.gz",
                "W6_ATAC_rep1_R1.fastq.gz",
                "W6_ATAC_rep1_R2.fastq.gz",
                "W6_ATAC_rep2_R1.fastq.gz",
                "W6_ATAC_rep2_R2.fastq.gz")

dir.create("/rugpfs/fs0/brc/scratch/brc_pipeline/external_fastqs/BRC/atac_sox9_fuchs_2023_mpaul")

sapply(1:length(atac_paths), function(x){
#sapply(5, function(x){
  print(x)
  if(!file.exists(file.path("/rugpfs/fs0/brc/scratch/brc_pipeline/external_fastqs/BRC/atac_sox9_fuchs_2023_mpaul", atac_names[x]))){
  download.file(atac_paths[x], file.path("/rugpfs/fs0/brc/scratch/brc_pipeline/external_fastqs/BRC/atac_sox9_fuchs_2023_mpaul", atac_names[x]))
}
  })

dir.create("/rugpfs/fs0/brc/scratch/brc_pipeline/external_fastqs/BRC/cnr_sox9_fuchs_2023_mpaul")

sapply(1:length(cnr_paths), function(x){
  download.file(cnr_paths[x], file.path("/rugpfs/fs0/brc/scratch/brc_pipeline/external_fastqs/BRC/cnr_sox9_fuchs_2023_mpaul", cnr_names[x]))
})

datasets <- rbind(
c("SOX9CNR_D0_rep1",	"YHBLDAMZIWFJXNRSVKOT"),
c("SOX9CNR_D0_rep2",	"ZYCLMWGQSVOHTNAIFPDX"),
  c("SOX9CNR_W6_rep1",	"LIQKNMDFWPVBJSTHOEYG"),
    c("SOX9CNR_W6_rep2",	"TWHECASNYQDIKUOXPVGB"),
      c("W6_ATAC_rep1",	"FOHZKJLBPGQUWNYMVSTC"),
        c("D0_ATAC_rep2",	"AJUTISKQHGLPWOMBVYDE"),
          c("D0_ATAC_rep1",	"WMVPNEXUHSFKALQZGOJD"),
            c("W6_ATAC_rep2",	"AFSWRTMKPEIDZBQLNHGY"))


