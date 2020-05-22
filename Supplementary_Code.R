# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     Code for    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                                                     #
#                               Impact of age of onset of cannabis use                                #
#                                  on expression of psychopathology                                   #
#                             in the general population: a network view.                              #
#                                                                                                     #
#                                     developed by L. Betz                                            #
#                                                                                                     #
#                               - Analysis reported in supplement -                                   #
#                                                                                                     #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ---------------------------------- 0: Reproducibility  -----------------------------------

# for reproducibility, one can use the "checkpoint" package
# in a temporal directory, it will *install* those package versions used when the script was written
# these versions are then used to run the script
# to this end, a server with snapshot images of archived package versions needs to be contacted
# for more info visit: https://mran.microsoft.com/documents/rro/reproducibility

library(checkpoint)
checkpoint(
  snapshotDate = "2020-04-01",
  R.version = "3.6.3",
  checkpointLocation = tempdir()
)

# ---------------------------------- 1: Load packages & data -----------------------------------
library(haven)
library(qgraph)
library(bootnet)
library(mgm)
library(tidyverse)


# the data set is available at https://www.icpsr.umich.edu/icpsrweb/

# NCS-1
X06693_0001_Data <- read_sav("DS0001/06693-0001-Data.sav")


# ---------------------------------- 2: Data preparation & sample descriptives -----------------------------------

variable_names <- c(
  "age of onset",
  "cumulative use",
  "childhood abuse",
  "childhood neglect",
  "urbanicity",
  "panic",
  "anxious",
  "sad",
  "loss interest",
  "irritable",
  "manic",
  "SPYING/FOLLOWING YOU",
  "POISON/HURT YOU",
  "READING YOUR MIND",
  "HEAR YOUR THOUGHTS",
  "HEAR OTHERS THOUGHT",
  "CONTROLLED BY FORCE",
  "OTHERS STOLE THOUGHTS",
  "SPECIAL MESSAGES/TV",
  "hypnotized/MAGIC/FORCE",
  "SAW VISIONS",
  "HEARD NOISE/VOICE",
  "SMELLS/BODY ODORS",
  "FEELNGS IN/ON BODY"
) %>% tolower(.)


data <- X06693_0001_Data %>%
  filter(V4097 <= 40) %>% # age at time of assessment <= 40
  filter(V1907 == 1) %>% # lifetime cannabis use; N here 2,624
  transmute(
    CASEID,
    age_onset = V1908,
    # 1
    cumulative_use = V1909,
    # 2
    # traumatic threat experiences: check if they happened in childhood (< 18 y/o)
    molested = ifelse(V6126 == 5, 0, ifelse(V6127 < 18, 1, 0)),
    raped = ifelse(V6114 == 5, 0, ifelse(V6115 < 18, 1, 0)),
    # physical abuse as a child
    physical_abuse = ifelse(V6143 == 5, 0, 1),
    
    
    # combine threat experiences in childhood into abuse variable
    abuse = ifelse(molested == 1 | raped == 1 |
                     physical_abuse == 1,
                   #  physical_assault == 1,
                   1,
                   0),
    
    neglect = ifelse(V6144 == 5, 0, 1),
    # neglect in childhood
    # neglect
    
    urbanicity = ifelse(V7136 == 5 |
                          V7136 == 4, 1, 0),
    # urbanicity: city or suburb
    # urbanicity
    V301,
    # panic: Have you ever in your life had a spell or attack when all of a sudden you felt frightened, anxious or very uneasy in situations when most people would not be afraid or anxious?
    V302,
    # anxious: Have you ever had a period of one month or more when most of the time you felt worried or anxious?
    V308,
    # sad: In your lifetime, have you ever had two weeks or more when nearly every day you felt sad, blue, or depressed?
    V310,
    # loss interest: Has there ever been two weeks or more when you lost interest in most things like work, hobbies, or things you usually liked to do for fun?
    
    
    V313,
    # irritable: Has there ever been a period of several days when you were so irritable that you threw or broke things, started arguments, shouted at people, or hit someone?
    V312,
    #mania: Has there ever been a period of at least two days when you were so happy or excited that you got into trouble, or your family or friends worried about it, or a doctor said you were manic?
    
    V4101,
    # psychotic experiences
    # 3
    V4103,
    # 4
    V4105,
    # 5
    V4108,
    # 6
    V4110,
    # 7
    V4112,
    # 8
    V4114,
    # 9
    V4116,
    # 10
    V4118,
    # 11
    V4120,
    # 12
    V4122,
    # 13
    V4133,
    # 14
    V4135 # 15
    
  ) %>%
  select(-c(molested, raped, physical_abuse)) %>% # drop these variables
  na.omit() %>% # only complete cases, as network estimator does not allow missings
  mutate_all(as.numeric) %>%
  mutate_at(vars(matches("V41|V3")), ~ recode(.,
                                              `5` = 0))

data_network <- data %>% select(-CASEID)
colnames(data_network) <- as.character(1:24)

# -------------------------------- 3: Network Estimation ---------------------------------

# we use a mixed graphical model as implemented in mgm-package
graph_all <- estimateNetwork(
  data_network,
  default = "mgm",
  criterion = "EBIC",
  tuning = 0,
  rule = "OR"
)

# export for supplementary table 1
write.table(
  round(graph_all$graph, 2),
  file = "edge_weights_network.csv",
  sep = ",",
  row.names = F,
  col.names = F
)
# ---------------------------------- 4: Bootstrapping -----------------------------------
# we bootstrap our model => it takes rather long (~2h per analysis) to run this
edge_boot <- bootnet(graph_all, nBoots = 1000, nCores = 6)


case_boot <-
  bootnet(
    graph_all,
    statistics = "edge",
    # we are not interested in centrality measures
    nBoots = 1000,
    type = "case",
    caseN = 10,
    nCores = 6
  )

# -------------------------------- 5: Plot bootstrap results -------------------------
# ---------- Supplementary figure 1 --------------
# 95% CI for edge weights obtained from bootstrapping

# NOTE:
# the following code is an adapation of the "plot.bootnet" function in the bootnet package
# it allows to split the plot depicting bootstrap confidence intervals for edges
# into adjacent sub-plots (assuming split0=T)
# Particularly for graphs with many nodes, this ensures readability and visability

x <-
  edge_boot # an object of type bootnet, resulting from bootstrapping
n_blocks <-
  5 # number of sections into which the plot should be devided

# some other values, we leave them mostly at default
meanVar <- "mean"
minArea <- "q2.5"
maxArea <- "q97.5"

minArea <- paste0(minArea, "_non0")
maxArea <- paste0(maxArea, "_non0")
meanVar <- paste0(meanVar, "_non0")
statistics  <- "edge"

prop0_minAlpha = 0.50
prop0_cex = 0.75
prop0_alpha = 0.5

sumTable <-
  summary(x, statistics = "edge", rank = F)  %>% ungroup %>% dplyr::mutate_(type = ~
                                                                              factor(type, levels = statistics))

# Summarize first:

summary <-
  sumTable %>% dplyr::group_by_( ~ id) %>% dplyr::summarize_(sample = ~
                                                               sample[type == statistics[[1]]], mean = as.formula(paste0("~mean(", meanVar, ",na.rm=TRUE)")))

summary$order <- order(order(summary$sample, summary$mean))

sumTable <-
  sumTable %>% dplyr::left_join(dplyr::select_(summary,  ~ id,  ~ order), by = "id")

# Reorder:
sumTable <- sumTable %>%
  dplyr::arrange_( ~ dplyr::row_number(order))  %>%
  dplyr::mutate_(id = ~ gsub("^(E|N): ", "", as.character(id))) %>%
  dplyr::mutate_(id = ~ factor(id, levels = unique(id)))


# Some fancy transformation:
revTable <- function(x)
  x[nrow(x):1, ]

sumTable$lbound <- sumTable[[minArea]]
sumTable$ubound <- sumTable[[maxArea]]


sumTable2 <- dplyr::bind_rows(
  sumTable %>% select_(
    ~ type,
    ~ id,
    ~ node1,
    ~ node2,
    ~ sample,
    ci = ~ lbound,
    mean = meanVar,
    ~ prop0
  ),
  revTable(
    sumTable %>% select_(
      ~ type,
      ~ id,
      ~ node1,
      ~ node2,
      ~ sample,
      ci = ~ ubound,
      mean = meanVar,
      ~ prop0
    )
  )
)

sumTable <-
  sumTable %>% left_join(
    sumTable2 %>% arrange(id, ci) %>% mutate(type_ci = rep(
      c("lower", "upper"),  max(x$bootTable$rank_avg)
    )) %>%
      reshape2::dcast(., id ~ type_ci, value.var = "ci"),
    by = "id"
  )

gathered_sumTable <-
  tidyr::gather_(sumTable, "var", "value", c("sample", meanVar))


gathered_sumTable$alpha.a <- ifelse(
  gathered_sumTable$var == meanVar,
  prop0_minAlpha + (1 - prop0_minAlpha) * (1 -
                                             gathered_sumTable$prop0),
  1
)
gathered_sumTable$alpha.b <-
  prop0_minAlpha + (1 - prop0_minAlpha) * (1 - sumTable2$prop0)

pdf("edge_weights_bootstrapped_CI.pdf", height = 8.27, width =  11.69)
gathered_sumTable %>%
  mutate(split_plot = -ntile(mean, n_blocks)) %>%
  ggplot(.,
         aes_string(
           x = 'value',
           y = 'id',
           group = 'id',
           colour = "var"
         )) +
  facet_wrap( ~ split_plot, scales = "free_y", nrow = 1) +
  geom_point(aes_string(alpha = 'alpha.a')) +
  geom_segment(
    aes_string(
      x = "lbound",
      y = 'id',
      xend = "ubound",
      yend = 'id',
      alpha = 'alpha.b',
      group = 'id'
    ),
    colour = "black"
  ) +
  geom_label(
    aes(x = 0, label = format(round(prop0, 2), nsmall = 2)),
    cex = prop0_cex * 3.1,
    label.padding = unit(0.1, "lines"),
    label.size = 0.1,
    alpha = prop0_alpha,
    colour = "black"
  ) +
  theme_bw() +
  xlab("") +
  ylab("") +
  scale_color_manual(
    "",
    values = c("black", "red"),
    labels = c("Bootstrap mean", "Sample")
  ) +
  theme(legend.position = "top") + 
  scale_alpha(guide = "none") +
  theme(
    axis.text.y = element_text(size = 9),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
dev.off()

# ---------- Supplementary figure 2 --------------
# plot a network in which only stable connections (>50% of bootstraps) are contained

# first, filter out connections that are not stable by this criterion
only_stable_edges <- gathered_sumTable[c("id", "prop0")] %>%
  right_join(
    data.frame(
      from = main_network$Edgelist$from,
      to = main_network$Edgelist$to,
      weight = main_network$Edgelist$weight
    ) %>%
      mutate(id = paste0(from, "--", to)),
    by = "id"
  ) %>% filter(prop0 <= 0.5) %>%
  distinct() %>%
  select(from, to, weight)

# use this information to plot the network and save it as a pdf in wd ("main_network_only_stable_edges.pdf")
qgraph(
  only_stable_edges,
  directed = F,
  layout = lay * -1,
  # flip everything because it looks nicer
  fade = ifelse(only_stable_edges$from == 1, FALSE, TRUE),
  trans = TRUE,
  color = c("#FED439", "#F38D81", "#7FC8F8", "#cad3db"),
  groups = c(
    rep("Cannabis Use", 2),
    rep("Early Risk", 3),
    rep("Mood", 6),
    rep("Psychosis", 13)
  ),
  theme = "Borkulo",
  legend = T,
  #minimum = 0,
  #maximum = 0.2,
  cut = 0,
  labels = 1:ncol(data_network),
  GLratio = 1.75,
  label.cex = 1.9,
  vsize = 3,
  legend.cex = 0.52,
  filename = "main_network_only_stable_edges",
  filetype = "pdf",
  edge.width = 0.5,
  nodeNames = variable_names
)

# to keep track: which edges are NOT included in at least 50% of bootstraps?
gathered_sumTable[c("id", "prop0")] %>%
  right_join(
    data.frame(
      from = main_network$Edgelist$from,
      to = main_network$Edgelist$to,
      weight = main_network$Edgelist$weight
    ) %>%
      mutate(id = paste0(from, "--", to)),
    by = "id"
  ) %>%
  filter(prop0 > 0.5) %>%
  distinct() %>%
  select(from, to, weight)

# ---------- Supplementary figure 3 --------------
# edge weight case-dropping bootstrap
pdf("case_drop_bootstrap.pdf",
    width = 7,
    height = 5)
plot(case_boot, statistics = "edge")
dev.off()

corStability(case_boot) # 0.672 for edge
