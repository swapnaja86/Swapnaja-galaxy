setwd("~/Swapnaja/Correlation data")
setwd("~/Swapnaja/Correlation data")
input_file <-"correlation-input-table.tab"
otu_variables_start <- 7
signf_cutoff <- 0.05
includeTax <- 0
includeMeta <- 0
fill_NA <- 0
replace_zeros <- 1
prevalence_exclusion <- 0.3
min_pair_support <- 4
plot_pval_cutoff <- 0.05
plot_corr_cutoff <- 0.5
packages <-c("Hmisc","corrplot")
InsPack <- function(pack)
{
if ((pack %in% installed.packages()) == FALSE) {
install.packages(pack)
}
}
lapply(packages, InsPack)
lib <- lapply(packages, require, character.only = TRUE)
flag <- all(as.logical(lib))
my_data <-
read.table (
file = input_file,
check.names = FALSE,
header = TRUE,
dec = ".",
sep = "\t",
row.names = 1,
comment.char = ""
)
my_data <- my_data[!apply(is.na(my_data) | my_data=="",1,all),]
fill_NA.mean <- function(vec)
{
m <- mean(vec, na.rm = TRUE)
vec[is.na(vec)] <- m
return(vec)
}
log_ratio <- function(data)
{
log_data <- log(data)
gm <- exp(mean(log_data[is.finite(log_data)]))
log_gm <- log(gm)
data <- log_data - log_gm
return(data)
}
my_data <- as.data.frame(apply(my_data,2,as.numeric))
first_OTU <- colnames(my_data)[otu_variables_start]
my_meta_data <- my_data[1:otu_variables_start - 1]
my_otu_data <- my_data[otu_variables_start:dim(my_data)[2]]
if (fill_NA == 0) {
my_meta_fixed =  my_meta_data
}
if (fill_NA == 1) {
my_meta_fixed = apply(my_meta_data, 2, fill_NA.mean)
}
prevalence_cutoff <- dim(my_otu_data)[1] - (prevalence_exclusion*dim(my_otu_data)[1])
na_count <-sapply(my_otu_data, function(y) sum(length(which(is.na(y)))))
zero_count <-sapply(my_otu_data, function(y) sum(length(which(y==0))))
prevalence_count <- na_count + zero_count
my_otu_data <- my_otu_data[, prevalence_count <= prevalence_cutoff ]
if (replace_zeros == 1) {
my_otu_data[my_otu_data==0] <- NA
}
my_otu_data[my_otu_data==0] <- 0.0001
my_otu_fixed = apply(my_otu_data, 2, log_ratio)
transformed_data <- cbind(my_meta_fixed, my_otu_fixed)
my_scaled_data <- scale(transformed_data, center = TRUE, scale = TRUE)
my_rcorr <- rcorr(as.matrix(my_scaled_data, type = "pearson"))
var_names <- row.names(my_rcorr$r)
if(includeTax==1 & includeMeta==0){
row_names <- var_names[c(otu_variables_start:dim(my_rcorr$r)[1])]
col_names <- var_names
pairs <-expand.grid(row_names, col_names)
my_cor_matrix <- my_rcorr$r[c(otu_variables_start:dim(my_rcorr$r)[1]),]
my_pvl_matrix <-my_rcorr$P[c(otu_variables_start:dim(my_rcorr$P)[1]),]
my_num_matrix <- my_rcorr$n[c(otu_variables_start:dim(my_rcorr$n)[1]),]
diagonale=0
} else if(includeTax==1 & includeMeta==1){
row_names <-var_names
col_names <- var_names
pairs <-expand.grid(row_names, col_names)
my_cor_matrix <- my_rcorr$r
my_pvl_matrix <-my_rcorr$P
my_num_matrix <- my_rcorr$n
diagonale=0
} else if (includeTax==0 & includeMeta==1) {
row_names <-var_names[c(1:(otu_variables_start - 1))]
col_names <- var_names
pairs <-expand.grid(row_names, col_names)
my_cor_matrix <- my_rcorr$r[c(1:(otu_variables_start - 1)),]
my_pvl_matrix <-my_rcorr$P[c(1:(otu_variables_start - 1)),]
my_num_matrix <- my_rcorr$n[c(1:(otu_variables_start - 1)),]
diagonale=1
} else {
row_names <- var_names[c(1:(otu_variables_start - 1))]
col_names <- var_names[otu_variables_start:dim(my_rcorr$r)[1]]
pairs <-expand.grid(row_names, col_names)
my_cor_matrix <- my_rcorr$r[c(1:(otu_variables_start - 1)),c(otu_variables_start:dim(my_rcorr$r)[1])]
my_pvl_matrix <-my_rcorr$P[c(1:(otu_variables_start - 1)),c(otu_variables_start:dim(my_rcorr$P)[1])]
my_num_matrix <- my_rcorr$n[c(1:(otu_variables_start - 1)),c(otu_variables_start:dim(my_rcorr$n)[1])]
diagonale=1
}
p_vector <- as.vector(my_pvl_matrix)
c_vector <- as.vector(my_cor_matrix)
n_vector <- as.vector(my_num_matrix)
my_pairs <-
matrix(ncol = 5,
c(
as.character(pairs[, 2]),
as.character(pairs[, 1]),
c_vector,
p_vector,
n_vector
))
my_pairs <- subset(my_pairs, as.numeric(my_pairs[,5]) > min_pair_support)
pVal_BH <- round(p.adjust(my_pairs[,4], method = "BH"), 4)
my_pairs <- cbind(my_pairs,as.numeric(pVal_BH))
my_pairs <- my_pairs[!as.character(my_pairs[, 1]) == as.character(my_pairs[, 2]),]
my_pairs <- my_pairs[!duplicated(my_pairs[, 3]), ]
matrix_names <- list(c(rep("",times=dim(my_pairs)[1])),
c(
"variable1",
"variable2",
"correlation",
"pValue",
"support",
"Corrected"
))
dimnames(my_pairs) <- matrix_names
my_pairs_cutoff <- my_pairs[as.numeric(my_pairs[, 4]) <= signf_cutoff, ]
my_pairs_cutoff <- matrix(my_pairs_cutoff,ncol=6,dimnames = list(c(rep("",times=dim(my_pairs_cutoff)[1])),c("variable1","variable2","correlation","pvalue","support","corrected pvalue")))
my_pairs_cutoff_corr <- my_pairs_cutoff[abs(as.numeric(my_pairs_cutoff[, 3])) >= 0.5, ]
my_cor_matrix <- my_cor_matrix[, colSums(is.na(my_cor_matrix)) != nrow(my_cor_matrix)]
my_cor_matrix[is.na(my_cor_matrix)] <- 0
OriginalPath <- getwd()
prefix = paste(strsplit(input_file,"[.]")[[1]][1],sep="_")
newdir <- paste(prefix,Sys.Date(), sep = "_")
dir.create(newdir)
setwd(newdir)
pdf(file = "corrplot.pdf")
corrplot(
matrix(data=na.omit(my_cor_matrix),nrow=dim(as.data.frame(row_names))[1], ncol = dim(as.data.frame(col_names))[1],dimnames=list(row_names,col_names)),
tl.col = "black",
tl.srt = 65,
tl.cex = 0.6,
cl.cex = 0.5,
diag=diagonale
)
dev.off()
if (plot_pval_cutoff != signf_cutoff | plot_corr_cutoff != 0.5) {
my_pairs_cutoff <- my_pairs[as.numeric(my_pairs[, 4]) <= plot_pval_cutoff, ]
corr_pval_cutoff <- my_pairs_cutoff[abs(as.numeric(my_pairs_cutoff[, 3])) >= plot_corr_cutoff, ]
corr_pval_cutoff <- matrix(corr_pval_cutoff,ncol=6, dimnames=list(c(rep("",times=dim(corr_pval_cutoff)[1])),c("variable1","variable2","correlation","pvalue","support","corrected pvalue")))
} else  {
corr_pval_cutoff <- my_pairs_cutoff_corr
corr_pval_cutoff <- matrix(corr_pval_cutoff,ncol=6, dimnames=list(c(rep("",times=dim(corr_pval_cutoff)[1])),c("variable1","variable2","correlation","pvalue","support","corrected pvalue")))
}
pdf(file = "linear_sign_pairs.pdf",family="sans",fillOddEven=TRUE)
for (i in 1:dim(corr_pval_cutoff)[1]) {
x_df <- transformed_data[names(transformed_data) %in% corr_pval_cutoff[i, 2]]
x <- x_df[, 1]
y_df <- transformed_data[names(transformed_data) %in% corr_pval_cutoff[i, 1]]
y <- y_df[, 1]
clm <- lm(y ~ x, na.action = na.exclude)
steps <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/1000
newx <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), steps)
a <- predict(clm, newdata = data.frame(x = newx), interval = "confidence")
plot(
x,
y,
ylab="",
xlab="",
cex.axis = 0.75 ,
xaxt = 'n',
yaxt = 'n',
xlim=c(min(x,na.rm = TRUE),max(x, na.rm = TRUE)),
xaxs="i"
)
title(xlab = names(x_df),line=0.5,font.lab=2,cex.lab=1.4)
title(ylab = names(y_df),line=0.5,font.lab=2,cex.lab=1.4)
polygon(c(newx,rev(newx)),c(a[,2],rev(a[,3])),col="grey91",border=TRUE, lty="dashed")
points(x,y)
abline(clm, lwd=2)
pvalue_text <- paste("P-value:", round(as.numeric(corr_pval_cutoff[i, 4]), 4), sep = "")
pvalue_corr_text <- paste("Adj. p-value:", as.numeric(corr_pval_cutoff[i, 6]), sep = "")
corr_text <- paste("Pearson's r:", round(as.numeric(corr_pval_cutoff[i, 3]), 4), sep = "")
support_text <- paste("supported by ", round(as.numeric(corr_pval_cutoff[i, 5]), 4)," observations", sep = "")
mtext(corr_text, side = 3, line = 2)
mtext(pvalue_text, side = 3, line = 1)
mtext(pvalue_corr_text, side = 3, line = 0)
mtext(support_text, side = 1, line = 2)
abline(par("usr")[3],0)
segments(par("usr")[1],a[1,2],par("usr")[1],a[1,3])
abline(par("usr")[4],0)
}
dev.off()
OriginalPath <- getwd()
prefix = paste(strsplit(input_file,"[.]")[[1]][1],sep="_")
newdir <- paste(prefix,Sys.Date(), sep = "_")
dir.create(newdir)
setwd(newdir)
write.table(my_scaled_data,"transformed.tab",sep = "\t",col.names = NA,quote = FALSE)
write.table(my_cor_matrix,"correlation-table.tab",sep = "\t",col.names = NA,quote = FALSE)
write.table(my_pvl_matrix,"pval-table.tab",sep = "\t",col.names = NA,quote = FALSE)
write.table(my_num_matrix,"support-table.tab",sep = "\t",col.names = NA,quote = FALSE)
write.table(my_pairs_cutoff,"cutoff-pairs-corr-sign.tab",sep = "\t",col.names = NA,quote = FALSE)
write.table(corr_pval_cutoff,"plotted-pairs-stat.tab",sep = "\t",col.names = NA,quote = FALSE)
if(!flag) { stop("
It was not possible to install all required R libraries properly.
Please check the installation of all required libraries manually.\n
Required libaries:ade4, GUniFrac, phangorn, randomcoloR, Rcpp")
}
savehistory("~/Swapnaja/updated Correlation script.R")
