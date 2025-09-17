library(susieR)
library(Matrix)

# load('output/coverages_mat_1trait_L1.RData')
# load('output/powers_mat_1trait_L1.RData')
# load('output/number_CSs_mat_1trait_L1.RData')
# load('output/number_SNPs_mat_1trait_L1.RData')

load('output/coverages_mat_1trait_L2.RData')
load('output/powers_mat_1trait_L2.RData')
load('output/number_CSs_mat_1trait_L2.RData')
load('output/number_SNPs_mat_1trait_L2.RData')


library(ggplot2)

list_mat = list(coverages, powers, no_CSs, no_SNPs)
list_name = c("Coverages", "Power", "Number of CSs", "Total number of SNPs")
plots = list()

num_rep = nrow(coverages)
rows_all_minus1 <- which(apply(no_CSs, 1, function(x) all(x == -1)))
rows_all_ok = setdiff(c(1:num_rep), rows_all_minus1)

library(boot)
for (i in 1:4){
  m = list_mat[[i]][rows_all_ok, c(1:5, 9, 10, 8)]
  m[m==-1] <- 0
  colnames(m) <- c('In-', 'Out-', 'L=1', 'reg.',
                   'P1(R)', '.9P1(R)+.1alphvv', '.5P1(R)+.5alphvv',
                   'alphavv+S')
  means <- colMeans(m)
  
  mean_fun <- function(data, idx) mean(data[idx])
  cis <- apply(m, 2, function(col) {
    if (length(unique(col)) == 1) {
      return(c(mean(col), mean(col)))
    } else {
      boot_res <- boot(col, mean_fun, R = 1000)
      return(boot.ci(boot_res, type = "perc")$percent[4:5])
    }
  })
  print(max(cis))
  df <- data.frame(
    variable = factor(colnames(m), levels = colnames(m)),
    mean = means,
    lower = cis[1, ],
    upper = cis[2, ]
  )
  
  plots[[i]] = ggplot(df, aes(x = variable, y = mean)) +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    labs(title = list_name[i],
         x = "Methods", y = "Mean with 95% CI")+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

library(patchwork)
wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "Performance when true L = 2"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )


