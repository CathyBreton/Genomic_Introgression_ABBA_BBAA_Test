# m_matschiner Mon Apr 16 12:47:43 CEST 2018

# Define functions to calculate the numbers of abba and baba patterns.
bbaa = function(p1, p2, p3, p4) p1 * p2 * (1 - p3) * (1 - p4)
abba = function(p1, p2, p3, p4) (1 - p1) * p2 * p3 * (1 - p4)
baba = function(p1, p2, p3, p4) p1 * (1 - p2) * p3 * (1 - p4)
D.stat = function(dataframe) (sum(dataframe$ABBA) - sum(dataframe$BABA)) / (sum(dataframe$ABBA) + sum(dataframe$BABA))
fd.stat = function(p1, p2, p3, p4) {
  pd = pmax(p2, p3)
  (sum(abba(p1, p2, p3, p4)) - sum(baba(p1, p2, p3, p4))) / (sum(abba(p1, pd, pd, p4)) - sum(baba(p1, pd, pd, p4)))
}

# Define functions for jackkniving.
get_genome_blocks <- function(block_size, lg_lengths) {
  block_starts <- sapply(lg_lengths, function(l) seq(1, l, block_size))
  data.frame(start = unlist(block_starts),
             end = unlist(block_starts) + block_size - 1,
             lg = rep(names(block_starts), sapply(block_starts, length)))
}
get_genome_jackknife_indices <- function(lg, position, block_info){
  lapply(1:nrow(block_info), function(x) !(lg == block_info$lg[x] &
                                             position >= block_info$start[x] &
                                             position <= block_info$end[x]))
}
get_jackknife_sd <- function(FUN, input_dataframe, jackknife_indices){
  n_blocks <- length(jackknife_indices)
  overall_mean <- FUN(input_dataframe)
  sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - FUN(input_dataframe[jackknife_indices[[i]],])*(n_blocks-1)))
}

# Get the command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
allele_freqs_file_name <- args[1]
output_file_name <- args[2]
spc_p1 <- args[3]
spc_p2 <- args[4]
spc_p3 <- args[5]
spc_o <- args[6]
lg_length_table_file_name <- args[7]
output2_file_name <- args[8]

# Read the allele-frequencies table.
freq_table = read.table(allele_freqs_file_name, header=T, as.is=T)
#freq_table = read.table("All_as_species.tsv.gz", header=T, as.is=T)


## In case of NAN
freq_table<- as.data.frame(freq_table)
freq_table[is.na(freq_table)] <- 0.0

#install.packages("data.table")
library(data.table)
setnames(freq_table, old=c("X3M05","X3M16","X3M24","X3M35","X3M36","X3M37","X3M38","X3M41" ), new=c("3M05","3M16","3M24","3M35","3M36","3M37","3M38","3M41"))

# Open the output file.
output_file <- file(output_file_name, "w")
#output_file <- file("output.csv", "w")

# Output.
write(paste("Species 1: ", spc_p1, sep=""), output_file, append=T)
write(paste("Species 2: ", spc_p2, sep=""), output_file, append=T)
write(paste("Species 3: ", spc_p3, sep=""), output_file, append=T)
write(paste("Species O: ", spc_o, sep=""), output_file, append=T)
write("", output_file, append=T)
write(paste("Number of sites: ", nrow(freq_table), sep=""), output_file, append=T)

# Get the allele frequencies.
p1 = freq_table[,spc_p1]
p2 = freq_table[,spc_p2]
p3 = freq_table[,spc_p3]
p4 = freq_table[,spc_o]



#ITC0629B	malacensis
#ITC0266B    AA
#ITC0621A    banksii
#ITC0856B    sckizocarpa

#P1 <- "pop1"
#P2 <- "pop2"
#P3 <- "pop3"
#P4 <- "Outgroup"

#p1 = freq_table[,P1]
#p2 = freq_table[,P2]
#p3 = freq_table[,P3]
#p4 = freq_table[,P4]

# Calculate the number of abba and baba patterns.
BBAA = bbaa(p1, p2, p3, p4)
ABBA = abba(p1,	p2, p3,	p4)
BABA = baba(p1,	p2, p3,	p4)

# Calculate the d statistic.
ABBA_BABA_df = as.data.frame(cbind(ABBA,BABA))
D = D.stat(ABBA_BABA_df)
#D

# Calculate Simon Martin's fd statistic.
fd = fd.stat(p1, p2, p3, p4)
#fd

# Output.
write(paste("Number of BBAA sites: ", sum(BBAA), sep=""), output_file, append=T)
write(paste("Number of ABBA sites: ", sum(ABBA), sep=""), output_file, append=T)
write(paste("Number of BABA sites: ", sum(BABA), sep=""), output_file, append=T)
write(paste("D statistic: ", D, sep=""), output_file, append=T)
write(paste("fd statistic: ", fd, sep=""), output_file, append=T)
if( sum(ABBA) > sum(BBAA)){
  cat(paste("\nWARNING: The number of ABBA sites (", sum(ABBA) , ") is greater than the number of BBAA sites (", sum(BBAA) , "), indicating that ", spc_p2, " and ", spc_p3, " are more closely related than ", spc_p1, " and ", spc_p2, ". You should swap ", spc_p1, " and ", spc_p3, " to get the correct D-statistic.\n\n", sep=""))
}

# Read the lg length table
lg_table = read.table(lg_length_table_file_name)
lg_table = read.table("Chromosome_length_AA.txt")
lg_lengths = lg_table[,2]
names(lg_lengths) = lg_table[,1]

# Run jackknife iterations to estimate the standard deviation in d.
#blocks = get_genome_blocks(block_size=1e6, lg_lengths=lg_lengths)
blocks = get_genome_blocks(block_size=1e7, lg_lengths=lg_lengths)
n_blocks = nrow(blocks)
indices = get_genome_jackknife_indices(lg=freq_table$scaffold,position=freq_table$position,block_info=blocks)
D_sd = get_jackknife_sd(FUN=D.stat, input_dataframe=as.data.frame(cbind(ABBA,BABA)),jackknife_indices=indices)

# Calculate the p-value for d.
D_err <- D_sd/sqrt(n_blocks)
D_Z <- D / D_err
D_p <- 2*pnorm(-abs(D_Z))

# Feedback.
write(paste("D_sd: ", D_sd, sep=""), output_file, append=T)
write(paste("D_err: ", D_err, sep=""), output_file, append=T)
write(paste("D_Z: ", D_Z, sep=""), output_file, append=T)
write(paste("p-value: ", D_p, sep=""), output_file, append=T)

#cat(paste("Pop1:",spc_p1," ", "Pop2:",spc_p2," ", "Pop3:",spc_p3," ", "Outgroup:",spc_o," ", "D_sd: ", D_sd, " ", "D_err: ", D_err," ", "D_Z: ", D_Z," ", "p-value: ", D_p, sep="" ), output_file, append=T)



# Close the output file.
close(output_file)

### Export to a table

#if(!file.exists(output2_file_name)) {
#  file.create(output2_file_name)
#}

## window size option 
options(width = 200)

# Feedback.
#cat(paste("\nWrote results to file ", output_file_name, ".\n\n", sep=""))

Data <- cbind(spc_p1, spc_p2, spc_p3, spc_o, D, D_sd, D_err, D_Z, D_p, fd, sum(ABBA), sum(BABA), sum(BBAA))
df_total <- data.frame(Data, dec=4)

### Header name change 
#names(df_total) <- c("Pop1", "Pop2", "pop3", "Outgroup", "D", "D_sd", "D_err", "D_Z", "D_p", "fd", "ABBA", "BABA", "BBAA")
names(df_total) <- NULL
 

# Derivation des sorties �cran vers le fichier standard
# Append to the file
sink(file = output2_file_name, append=TRUE, split =TRUE) #permet d'ex�cuter un script contenu dans un fichier externe

df_total


# Close the file.
sink()


