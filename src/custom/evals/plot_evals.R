library(dplyr)
library(stringr)
library(conflicted)
library(ggplot2)
library(data.table)
library(tidyr)

conflict_prefer("filter", "dplyr")

dir.create("data/figures/roc", recursive = TRUE)
dir.create("data/figures/pr", recursive = TRUE)

# import confusion matrices -----------------------------------------------------------------------
model_list <- lapply(list.files("data/evaluate_models/", full.names = TRUE), read.table, sep = "\t", header = TRUE)
model_df <- bind_rows(model_list) %>% mutate(sample = paste(condition, replicate, mark, sep = "_"))

print(head(model_df))

model_df$replicate <- as.factor(model_df$replicate)
model_df <- model_df %>% 
  mutate(f1 = 2 * (TP * FP) / (TP + FP)) %>%
  mutate(info = paste(method, condition, replicate, mark, sep = "_")) %>%
  mutate(info = as.factor(info))

# plot + export ROC -------------------------------------------------------------------------------
group_keys <- model_df %>% group_by(condition, mark) %>% group_keys()
print("Group Keys")
print(group_keys)

for (i in 1:nrow(group_keys)) {
  
  # define condition, mark
  tmp_condition <- as.character(group_keys[i,1])
  tmp_mark <- as.character(group_keys[i,2])
  
  # file I/O
  out_roc = paste0("data/figures/roc/", tmp_condition, "_", tmp_mark, ".roc.pdf") # \m/
  out_pr = paste0("data/figures/pr/", tmp_condition, "_", tmp_mark, ".pr.pdf")
  
  # filter the main DF with above groupings
  tmp_df <- model_df %>%
    filter(condition == tmp_condition) %>%
    filter(mark == tmp_mark)
  
  # plot ROC and PR
  tmp_roc <- ggplot(tmp_df, aes_string(x = "fpr", y = "recall", color = "method")) +
    geom_path(size = 2) +
    xlim(c(0,1)) +
    ylim(c(0,1)) +
    xlab("False Positive Rate") +
    ylab("Recall") +
    geom_abline(slope = 1, linetype = 2) +
    facet_wrap(sample ~ .) +
    scale_color_brewer(palette = "RdYlBu", direction = -1) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"))
  
  tmp_pr <- ggplot(tmp_df, aes_string(x = "recall", y = "precision", color = "method")) +
    geom_path(size = 2) +
    xlab("Recall") +
    ylab("Precision") +
    facet_wrap(sample ~ .) +
    scale_color_brewer(palette = "RdYlBu", direction = -1) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white")) +
    scale_x_continuous(breaks = c(seq(0.05, 0.95, 0.20)), limits = c(0.05, 0.95))
  
  ggsave(out_roc, tmp_roc, width = 16, height = 9, dpi = 600)
  ggsave(out_pr, tmp_pr, width = 16, height = 9, dpi = 600)
}

# AUC unit test -----------------------------------------------------------------------------------

# straight line is a triangle from x=0,y=0 to x=1,y=1.
# the area of the triangle should equal the sum of its trapezoids (AUC).
# area of triangle = 1/2 * (1) * (1) = 0.5
straight_line <- data.frame(x = seq(0, 1, 0.01), y = seq(0, 1, 0.01))
plot(straight_line$x, straight_line$y, type = "l")

trapezoid_area <- function(a, b, h) {
  return( 1/2 * (a + b) * h )
}

AUC <- function(roc_coordinates) {
  
  # input: roc coordinates in 2 cols: x and y values. Must be sorted by x-axis (ascending). Col order matters.
  # method: loop over each x value, define y(x), and x+i and y(x+i) where x+i is the next x-coord.
  # method: calculate trapezoidal area, sum all the trapezoidal areas.
  # output: area under the curve (float)
  
  auc = 0
  for (i in 1:nrow(roc_coordinates)) {
    
    # walk along the x-axis, define current point (x) and the next point (x_i, or x+i)
    # where i is diff btw two x coords. Keep track of the y values along the curve.
    x <- roc_coordinates[i, 1]
    x_i <- roc_coordinates[i + 1, 1]
    y_of_x <- roc_coordinates[i, 2]
    y_of_x_i <- roc_coordinates[i + 1, 2]
    i <- x_i - x

    # if at the end of the ROC curve, calculate the area one last time and add to total auc
    # then lastly stop the calculations and break the loop
    if (x_i == max(roc_coordinates[,1])) {
      
      last_area <- trapezoid_area(y_of_x, y_of_x_i, i)
      auc <- auc + last_area
      break
      
    }
    
    # Calculate trapezoid and add it to the total
    tmp_area <- trapezoid_area(y_of_x, y_of_x_i, i)
    auc <- auc + tmp_area
    
  }
  
  if (auc > 1) {
    stop("ERROR: AUC is greater than 1. Please review your data!")
  }
  
  return(auc)
  
}

if (AUC(straight_line) == 0.5) {
  print("Trapezoidal AUC Method worked for a right triangle!")
}

# AUC for ROCs ------------------------------------------------------------------------------------

# Calculate AUC for all method,condition,replicate,mark
list_of_auc <- lapply(split(model_df, model_df$info), function(x) {

  x %>% 
    select(fpr, recall) %>%
    arrange(fpr) %>%
    filter(fpr <= 1) %>%
    AUC()

})

auc_table <- gather(as.data.frame.list(list_of_auc)) %>%
  mutate(key = str_replace(key, "\\.", "-")) %>%
  separate(key, into = c("method", "condition", "replicate", "mark"), sep = "_") %>%
  rename("AUC" = value)

auc_graph <- auc_table %>% 
  ggplot(aes(x = replicate, y = AUC, fill = method)) +
    geom_col(position = "dodge") +
    facet_wrap(condition + mark ~ .) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    ggtitle("Area Under the ROC Curves by Condition and Mark") +
    geom_text(aes(label = round(AUC, 3), y = AUC - 0.05), position = position_dodge(0.9), vjust = 0, color = "white") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"))

out_auc = "data/figures/roc/auc.pdf"
ggsave(out_auc, auc_graph, width = 16, height = 9, dpi = 600)
