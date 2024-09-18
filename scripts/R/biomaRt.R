# Load the biomaRt package
library(biomaRt)
listMarts()
setwd('/Users/juanjovel/jj/temp')
rm(list=ls())

ensembl=useMart("ensembl")
listDatasets(ensembl)

# Set the dataset to use
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

outfile <- "hsapiens_transcript_attributes_240908.tsv"

# Set the attributes to retrieve
# Make list of attributes to retrieve
attributes <- c("ensembl_transcript_id",
                "ensembl_gene_id",
                "external_gene_name",
                "entrezgene_id",
                "wikigene_description",
                "name_1006",
                "definition_1006",
                "namespace_1003")

# Get the transcripts and their attributes
transcripts_annotations <- getBM(attributes = attributes,
                     mart = mart)

# Write the results to a tsv file
write.table(transcripts_annotations, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)