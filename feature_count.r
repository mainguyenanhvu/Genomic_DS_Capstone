sessionInfo()
f = read.csv("Week5QCdata.csv")
phenotype_table = f
rownames(phenotype_table) = phenotype_table[,1]
phenotype_table[,1] = NULL
#head(phenotype_table)

write.table(phenotype_table, file="phenotype.txt", col.names=TRUE, row.names=TRUE)

adult = f[1:3,]
fetal = f[4:6,]

# summary(adult)

# summary(fetal)

# t.test(fetal$alignment.rate, adult$alignment.rate)

# t.test(fetal$Average.Quality.per.read, adult$Average.Quality.per.read)

library('tidyverse',quietly=TRUE)
library(org.Hs.eg.db,quietly=TRUE)
library(annotate,quietly=TRUE)
# read feature count files
tabular_files = list.files(path = "./Data/FeatureCount", pattern = "tabular$", full.names = TRUE)
tabular_list = lapply(tabular_files, function(x) { read.table(x,stringsAsFactors=F)})

header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
tabular_list = lapply(tabular_list,header.true)

# merge the files by Geneid
feature_count_files = Reduce(function(x, y) merge(x, y, by="Geneid"), tabular_list)
#colnames(feature_count_files) = c("Geneid", "SRR1554536", "SRR1554538", "SRR1554556",  "SRR1554561", "SRR1554566", "SRR1554568")
#feature_count_files = feature_count_files[c("Geneid", "SRR1554536", "SRR1554538", "SRR1554556", "SRR1554561", "SRR1554566", "SRR1554568")]
# head(feature_count_files)
# print(typeof(feature_count_files))
# print(typeof(feature_count_files[1,1]))
#convert Geneid to gene_name
for (i in 1:nrow(feature_count_files)){
  #print(toString(lookUp(toString(feature_count_files[i,1]), 'org.Hs.eg', 'SYMBOL')))
  #print(typeof(toString(lookUp(toString(feature_count_files[i,1]), 'org.Hs.eg', 'SYMBOL'))))
  feature_count_files[i,1] = toString(lookUp(toString(feature_count_files[i,1]), 'org.Hs.eg', 'SYMBOL'))
  break
}
# head(feature_count_files)
#feature_count_files[,1] = getSYMBOL(feature_count_files[,1],data='org.Hs.eg')
rownames(feature_count_files) = make.names(feature_count_files[,1], unique=TRUE)
feature_count_files[,1] = NULL
feature_table = feature_count_files

head(feature_table)

write.table(feature_table, file="./feature_counts.txt", sep='\t', row.names=TRUE, col.names=TRUE)