library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

df = read.table(snakemake@input[[1]], sep = "\t", header = TRUE)

p = df %>%
  mutate(label = paste0(intersections, "\n", "(", percent, ")")) %>%
  ggplot(aes(x = method, y = intersections, fill = method)) +
  geom_col() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.85), col = "white") +
  facet_wrap(condition + mark ~ ., scales = "free_y") +
  scale_fill_brewer(palette = "RdYlBu", direction = -1) +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white")) +
  labs(title = "Consensus Peak Intersections with Gold Standards")
p

ggsave(snakemake@output[[1]], p, width = 16, height = 9, dpi = 600)