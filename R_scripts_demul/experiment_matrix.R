# Create a simple 4 (genes) x 3 (cells) matrix
dense_mat <- matrix(c(
  0, 5, 0,   # Gene1
  1, 0, 0,   # Gene2
  0, 2, 3,   # Gene3
  0, 0, 7    # Gene4
), nrow = 4, ncol = 3, byrow = TRUE)

rownames(dense_mat) <- c("Gene1", "Gene2", "Gene3", "Gene4")
colnames(dense_mat) <- c("Cell1", "Cell2", "Cell3")

dense_mat

library(Matrix)

sparse_mat <- as(dense_mat, "dgCMatrix")
sparse_mat
str(sparse_mat)

