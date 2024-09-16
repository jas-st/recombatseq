library(recombatseqv2)

count_red <- read.csv("tutorial/reduced_counts.csv", row.names = 1, header = TRUE)
count_red <- as.matrix(count_red)
count_batches <- read.csv("tutorial/batches.csv", row.names=1)
count_batches <- count_batches$batches

write.csv(count_recombseq, "count_recombseq.csv")


counts <- matrix(rnbinom(400, size=10, prob=0.1), nrow=50, ncol=8)
#colnames(counts) <- paste(rep("Sample", 8),1:8)
#rownames(counts) <- paste(rep("Gene", 50), 1:50)
batch2 <- c(rep(1, 4), rep(2, 4))

#adjusted2 <- ComBat_seq(cts_sub, batch=batch_sub, group=group_sub, shrink=FALSE)
ComBat_seq(counts, batch=batch2, shrink=FALSE)
combatseq_sub <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub, shrink=FALSE)

# how its supposed to work
combat_original <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub, shrink=FALSE)
cts_adjori_norm <- apply(combat_original, 2, function(x){x/sum(x)})

seobj_adjori <- SummarizedExperiment(assays=cts_adjori_norm, colData=col_data)
pca_obj_adjori <- plotPCA(DESeqTransform(seobj_adjori), intgroup=c("Batch", "Group"))
plt_adjori <- ggplot(pca_obj_adjori$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adjori$plot_env$percentVar[2])),
       title="Original ComBat")
plt_adjori


# test
combatseq_sub <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub, shrink=FALSE)
cts_adj_norm <- apply(combatseq_sub, 2, function(x){x/sum(x)})

seobj_adj <- SummarizedExperiment(assays=cts_adj_norm, colData=col_data)
pca_obj_adj <- plotPCA(DESeqTransform(seobj_adj), intgroup=c("Batch", "Group"))
plt_adj <- ggplot(pca_obj_adj$data, aes(x=PC1, y=PC2, color=Batch, shape=Group)) +
  geom_point() +
  labs(x=sprintf("PC1: %s Variance", percent(pca_obj_adj$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obj_adj$plot_env$percentVar[2])),
       title="ComBat-Seq")
plt_adj
