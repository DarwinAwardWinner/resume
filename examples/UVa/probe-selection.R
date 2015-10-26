# We will use base 2 logarithms for now
log.base = 2

# Read the data in from the file
filename.data = "Processed_Data_Files/Normalized_Pair_Files/All_norm_pair.txt"
#filename.data = "Raw_Data_Files/Pair_Files/All_pair.txt"
intensity.lin <- read.table(file=filename.data, header=TRUE,row.names="PROBE_ID")

# Column names and numbers used to categorize the data
intensity.maincols = c("GENE_EXPR_OPTION", "SEQ_ID", "POSITION")
intensity.maincolnums = which(names(intensity.lin) %in% intensity.maincols)
# Column names and numbers containing intensity data
intensity.datacolnums = (1:dim(intensity.lin)[2])[-intensity.maincolnums] # i.e. "the rest"
intensity.datacols = names(intensity.lin)[intensity.datacolnums]

# Take the logarithm of the data
intensity.log <- data.frame(intensity.lin[intensity.maincols],log(intensity.lin[intensity.datacols], base = log.base))

# Separate into random (i.e. background) and data probes
# Split into "intensity$rand" and "intensity$data"
intensity <- split(x = intensity.log, f = (ifelse(intensity.log$GENE_EXPR_OPTION == "RANDOM","rand","data")))
rm(intensity.lin, intensity.log) # Discard unused stuff to free memory

# This gives a QQ plot of the data against the noise. I used it to select my cutoffs.
#qqplot( y=as.vector(as.matrix(intensity$data[sample(1:dim(intensity$data)[1],5000),intensity.datacols])), x=as.vector(as.matrix(intensity$rand[intensity.datacols])), ylab="Data", xlab="Random Control",main = "Log10 QQ plot of Specific Probe Intensities vs. Random Controls")

# Detection (low) threshold is 2 sd above mean random background
intensity.rand.vector = as.vector(as.matrix(intensity$rand[intensity.datacols]))
threshold.low = mean(intensity.rand.vector) + 2*sd(intensity.rand.vector)
rm(intensity.rand.vector)
# Saturation (high) threshold is 2-fold down from max possible
threshold.high = log(65535/2, base=log.base)
if (threshold.high <= threshold.low) print("Error: low threshold is too high")

# Compute needed row statistics
intensity$data$MAX <- apply(intensity$data[intensity.datacols],1,max)
# Actually, we only need the max, it seems. Uncomment these if needed.
#intensity$data$MIN <- apply(intensity$data[intensity.datacols],1,min)
#intensity$data$MEAN <- rowMeans(intensity$data[intensity.datacols]) # Don't need that one
#intensity$data$SD <- apply(intensity$data[intensity.datacols],1,sd) # Don't need that one

# Sort probes into three bins: absent, present, and saturated
# The integer value of BIN also serves as a rank:
# present < saturated < absent; lower is better
intensity$data$BIN <- factor(1 + (intensity$data$MAX > threshold.high) + (intensity$data$MAX < threshold.low) * 2, labels = c("present", "saturated", "absent"))

# Count how many probes from each probeset (SEQ_ID) went into each bin
# Also count total probes per set
# Counting is done by measuring length of data aggregated by SEQ_ID
num.probes <- lapply(c(list(total=intensity$data),split(x = intensity$data, f = intensity$data$BIN)), FUN = function (y) aggregate(x=rep(NA,dim(y)[1]),by=list(SEQ_ID=y$SEQ_ID),FUN=length))
if(mean(num.probes$present$x) < 3) { print("Error: Not enough present probes.") }


# Something below this needs updating



# 2 rankings: Nimblegen's and distance from probeset mean
# Nimblegen's rank is read from the design file
# Probeinfo file
# This is information parsed from the design file.
# The design file can't be used directly because info on probe
# selection is not in separate columns.
filename.probeinfo = "Design_Files/071031_U_Va_Tobacco_Expr.probeinfo"
probeinfo <- read.table(filename.probeinfo,header=TRUE,row.names="PROBE_ID",as.is="SEQ")
# Use probe names as row names for indexing
#row.names(probeinfo) <- as.character(probeinfo$PROBE_ID)

# Add Nimblegen rank to the main data frame
intensity$data[c("RANK","SEQ")] <- probeinfo[row.names(intensity$data),c("RANK","SEQ")]


# For present probes, rank is based on correlation to probeset mean

# Read probeset means from calls file
filename.calls = "Processed_Data_Files/Normalized_Calls_Files/All_norm_calls.txt"
calls.lin <- read.table(filename.calls,header=TRUE,row.names="SEQ_ID")
# This line simultaneously logs the data and gives the same column order as the intensity table
probeset.means = log(calls.lin[intensity.datacols], base=log.base)
rm(calls.lin)

# Collect the relevant data into matrices for efficiency
probes.present.data = t(intensity$data[intensity$data$BIN == "present",intensity.datacols])
probeset.means.data = t(probeset.means[intensity$data$SEQ_ID[intensity$data$BIN == "present"],])

# We invert the correlation so that lower is better
intensity$data$RANK[intensity$data$BIN == "present"] <- sapply(1:dim(probes.present.data)[2],function (x) { -cor(probes.present.data[,x],probeset.means.data[,x]) })

# Done with these
rm("probes.present.data","probeset.means.data")

# Sort by bin, then rank
intensity$data.ranked = intensity$data[order(intensity$data$BIN,intensity$data$RANK),]

# Split by SEQ_ID
probes.ranked <- split(x=row.names(intensity$data.ranked),
                       f=intensity$data.ranked$SEQ_ID,
                       drop=TRUE)

# Set the desired number of probes per sequence
num.probes.desired <- 3

# A function to make any vector have length n, by truncating longer ones and padding shorter ones with NA
firstN <- function (v,n) c(v,rep(NA,n))[1:n]
# A function to take a string and append P1, P2, P3, etc. up to PN.
probenamesN <- function (s,n) (paste(s,"P",1:n,sep=""))[1:n]

probes.selected <- c(sapply(probes.ranked,firstN,num.probes.desired))
names(probes.selected) <- c(sapply(names(probes.ranked),probenamesN,num.probes.desired))
# Filter NA
probes.selected <- probes.selected[!is.na(probes.selected)]

probes.selected.fasta <- paste(sep="",">",names(probes.selected),"\n",intensity$data[probes.selected,"SEQ"])

cat(sep="\n",probes.selected.fasta,file="selected_probes.fasta")

save.image(file="probesel.rda")
