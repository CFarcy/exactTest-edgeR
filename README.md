# exactTest-edgeR
# launch command for the script: 
Rscript exactTest.R countTable.txt working/directory threshold p-value

/*- countTable.txt is the complete path of your file containing the read count
 - working/directory is the output directory where R will work and create the output
- threshold is the minimum of read you are allowing for a gene (all condition included)
 - p-value used to filter your genes predicted as differentially expressed
