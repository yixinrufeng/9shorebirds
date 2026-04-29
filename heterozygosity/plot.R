setwd("/path/angsd")

# install.packages("rcompanion")
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(rcompanion)

species_files <- list(
  bar = c(
    "NNU_ACE_0035.0", "NNU_ACE_0081.0", "NNU_ACE_0082.0", "NNU_ACE_0083.0",
    "NNU_ACE_0084.0", "NNU-ACE-0010.0", "SS_NNU_ACE_0199.0",
    "SS_NNU_ACE_0200.0", "SS_NNU_ACE_0201", "SS_NNU_ACE_0202.0"
  ),
  common = c(
    "NNU_ACE_0018.0", "NNU_ACE_0042.0", "NNU_ACE_0043.0", "NNU_ACE_0044.0",
    "NNU_ACE_0045.0", "NNU_ACE_0046.0", "NNU_ACE_0047.0", "NNU_ACE_0048.0",
    "NNU-ACE-0012.0", "SS_NNU_ACE_0049.0"
  ),
  eastern = c(
    "NNU_ACE_0020.0", "NNU_ACE_0021.0", "NNU_ACE_0022.0", "NNU_ACE_0023.0",
    "NNU_ACE_0024.0", "NNU_ACE_0068.0", "NNU_ACE_0069.0", "NNU_ACE_0070.0",
    "NNU-ACE-0005.0"
  ),
  great = c(
    "NNU_ACE_0025.0", "NNU_ACE_0026.0", "NNU_ACE_0027.0", "NNU_ACE_0032.0",
    "NNU_ACE_0071.0", "NNU-ACE-0004.0", "S-NNU-ACE-0028.0",
    "S-NNU-ACE-0029.0", "S-NNU-ACE-0030.0", "SS_NNU_ACE_0198.0"
  ),
  grey = c(
    "NNU_ACE_0085.0", "NNU_ACE_0086.0", "NNU_ACE_0087.0", "NNU_ACE_0088.0",
    "NNU_ACE_0089.0", "NNU_ACE_0090.0", "NNU_ACE_0091.0", "NNU_ACE_0092.0",
    "NNU_ACE_0093.0", "NNU-ACE-0011.0"
  ),
  red = c(
    "NNU_ACE_0072.0", "NNU_ACE_0073.0", "NNU_ACE_0074.0", "NNU_ACE_0076.0",
    "NNU_ACE_0078.0", "NNU_ACE_0079.0", "NNU_ACE_0080.0", "NNU-ACE-0009.0",
    "S-NNU-ACE-0075.0", "SS_NNU_ACE_0077.0"
  ),
  ruddy = c(
    "NNU_ACE_0050.0", "NNU_ACE_0051.0", "NNU_ACE_0052.0", "NNU_ACE_0053.0",
    "NNU_ACE_0054.0", "NNU_ACE_0055.0", "NNU_ACE_0056.0", "NNU_ACE_0057.0",
    "NNU_ACE_0058.0", "NNU-ACE-0007.0"
  ),
  spoon = c(
    "NNU_ACE_0002.0", "NNU_ACE_0003.0", "NNU-ACE-0001.0",
    "SAMN08968812", "SAMN08968815", "SAMN08968816", "SAMN08968817",
    "SAMN08968818", "SAMN08968819", "SAMN08968824"
  ),
  terek = c(
    "NNU_ACE_0059.0", "NNU_ACE_0060.0", "NNU_ACE_0061.0", "NNU_ACE_0062.0",
    "NNU_ACE_0063.0", "NNU_ACE_0064.0", "NNU_ACE_0065.0", "NNU_ACE_0066.0",
    "NNU_ACE_0067.0", "NNU-ACE-0008.0"
  )
)

species_labels <- c(
  bar = "Bar-tailed Godwit",
  common = "Common Snipe",
  eastern = "Far Eastern Curlew",
  great = "Great Knot",
  grey = "Grey-tailed Tattler",
  red = "Red-necked Stint",
  ruddy = "Ruddy Turnstone",
  spoon = "Spoon-billed Sandpiper",
  terek = "Terek Sandpiper"
)

species_cols <- c(
  bar = "#EBC15B",
  common = "#F9D71C",
  eastern = "#D1A146",
  great = "#ECF700",
  grey = "#AF780A",
  red = "#A6E800",
  ruddy = "#267300",
  spoon = "#5DB300",
  terek = "#C04800"
)

## 读入数据
dat <- read.csv("resuangsd.csv", header = FALSE, stringsAsFactors = FALSE)
colnames(dat) <- c("sample", "value")

## 去掉样本名末尾的 .ml
dat$sample <- sub("\\.ml$", "", dat$sample)

## 构建样本-物种映射表
sample_to_species <- unlist(
  lapply(names(species_files), function(sp) {
    x <- sub("\\.0$", "", species_files[[sp]])  # 去掉末尾 .0
    names(x) <- x
    rep(sp, length(x))
  })
)
names(sample_to_species) <- unlist(lapply(species_files, function(x) sub("\\.0$", "", x)))

## 给数据加上 species 列
dat$species <- sample_to_species[dat$sample]

## 检查未匹配样本
unmatched <- dat %>% filter(is.na(species))
print(unmatched)
write.csv(unmatched, "unmatched_samples.csv", row.names = FALSE)

## 只保留已匹配样本
matched <- dat %>% filter(!is.na(species))

## 设置物种顺序
species_order <- c("spoon", "great", "eastern", "red", "ruddy", "bar", "grey", "terek", "common")

matched$species <- factor(matched$species, levels = species_order)

## 总体检验
kw_res <- kruskal_test(matched, value ~ species)
print(kw_res)
write.csv(kw_res, "kruskal_test_result.csv", row.names = FALSE)

## 两两比较
dunn_res <- dunn_test(matched, value ~ species, p.adjust.method = "BH")
print(dunn_res)
write.csv(dunn_res, "pairwise_dunn_test_result.csv", row.names = FALSE)

## 生成显著性字母
letter_df <- dunn_res %>%
  mutate(comparison = paste(group1, group2, sep = "-")) %>%
  select(comparison, p.adj) %>%
  cldList(p.adj ~ comparison, data = ., threshold = 0.05)

print(letter_df)

## 给字母添加位置
label_df <- matched %>%
  group_by(species) %>%
  summarise(
    y = boxplot.stats(value)$stats[5],   # 上须位置，比 max(value) 更贴近箱线图
    .groups = "drop"
  ) %>%
  left_join(letter_df, by = c("species" = "Group"))

y_offset <- 0.02 * diff(range(matched$value, na.rm = TRUE))
label_df$y <- label_df$y + y_offset
## 保证 label_df 的物种顺序和 matched 一致
label_df$species <- factor(label_df$species, levels = species_order)

## 作图
y_min <- floor(min(matched$value, na.rm = TRUE) / 0.0005) * 0.0005
y_max <- ceiling(max(label_df$y, na.rm = TRUE) / 0.0005) * 0.0005

ggplot(matched, aes(x = species, y = value, fill = species)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85) +
  geom_text(
    data = label_df,
    aes(x = species, y = y, label = Letter),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = species_cols) +
  scale_x_discrete(labels = species_labels) +
  scale_y_continuous(
    breaks = seq(y_min, y_max, by = 0.0005),
    minor_breaks = seq(y_min, y_max, by = 0.00025)
  ) +
  theme_classic(base_size = 14) +
  labs(x = "Species", y = "Heterozygosity") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(y_min, y_max))
