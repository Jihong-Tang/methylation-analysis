#This script is used to transfer the bismark result .cov file to facilitate the downstream analysis using package DSS. The original .cov file should be transfered to .txt file to facilite the data importing.
#Argument: Command argument:[1] The path of the file to be transfered
#Author: Jihong Tang
#Date: April 25, 2019

library(tidyr)
library(dplyr)

file_names <- commandArgs(T)

func_read_file <- function(file_name){
	dir_vec <- strsplit(file_name, split = "/")[[1]]
	len <- length(dir_vec)
	file_prefix = substring(dir_vec[len], 0, nchar(dir_vec[len]) - 4)
	file_save_path = substring(file_name, 0, nchar(file_name) - nchar(dir_vec[len]))
	print(paste("File", file_name, "is being importing and this may take a while..."), sep = "")
	rawdata_df <- read.table(file_name, header = F, stringsAsFactors = F)
	print("Importing file is finished!")
	colnames(rawdata_df) <- c("chr", "start", "end", "methyl%", "methyled", "unmethyled")
	write.table(rawdata_df, paste(file_save_path, file_prefix, "_transfered.txt", sep = ""), row.names = F )
}

lapply(file_names, func_read_file)

