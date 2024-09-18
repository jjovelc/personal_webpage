library (dplyr)
library (ggplot2)
library (wesanderson)

metaFile <- 'metadata.txt'
metadata <- read.table(metaFile, header = T, sep = '\t', row.names = 1)
counts_data <- read.table(file, header = T, sep = '\t', row.names = 1)
counts_data_df <- as.data.frame(counts_data)
head(counts_data_df)

# Import a list of transcripts IDs that were found differentially expressed as above
transc_ids <- read.table (ids_file, header = F, sep = "\t")

counter = 0

for (i in 1:nrow(transc_ids)) {
  
  counter = counter + 1
  trID <- transc_ids [i,]
  atch <- grep (trID, readLines(file), value = T)
  
  sel_row   <- counts_data_df[trID,]
  row_names <- colnames(counts_data_df)
  group     <- metadata$group
  counts <- as.numeric(sel_row[1:ncol(sel_row)])

  df = data.frame(row.names = row_names, 
                  condition = group,
                  counts = counts)
  
  gene <- row.names(sel_row)
  gene.name <- paste(gene, "Expression", sep = " ")
  fileName <- paste (gene, "pdf", sep = ".")
  message <- paste("Plot#:", counter, "Current plot:", fileName, sep = ' ')
  print (message )
  pdf(file = fileName)
  p <- ggplot(df, aes(condition, counts))  
    
  print (p + geom_boxplot(notch=F,  aes(fill = factor(group))) + theme_bw() +
    #scale_fill_brewer(palette="Dark2") +
      scale_fill_manual(values=wes_palette(n=2, name="Darjeeling2")) +
      ggtitle (gene.name) +
      xlab("Condition") +
      ylab("Abundance (TPM)") +
      #scale_y_log10 () +
      # Control the text features
      theme( axis.text.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
           axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),  
           axis.title.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
           axis.title.y = element_text(colour="black",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))
  )
  dev.off()
}
