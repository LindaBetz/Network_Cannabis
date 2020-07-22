# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     Code for    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                                                     #
#                               A network approach to relationships                                   #
#                      between cannabis use characteristics and psychopathology                       #
#                                       in the general population                                     #
#                                                                                                     #
#                                     developed by L. Betz                                            #
#                                                                                                     #
#                         - Analysis reported in Supplementary Material -                             #
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


# the data set is available for public use at https://www.icpsr.umich.edu/web/ICPSR/studies/6693

# NCS-1
X06693_0001_Data <- read_sav("DS0001/06693-0001-Data.sav")


# ---------------------------------- 2: Data preparation & sample descriptives -----------------------------------

variable_names <- c(
  "age of cannabis use initiation",
  "lifetime cumulative frequency",
  "childhood abuse",
  "childhood neglect",
  "urban upbringing",
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
  "FEELINGS IN/ON BODY"
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

data_network <-
  data %>% select(-CASEID) # drop participant ID for network analysis
colnames(data_network) <-
  as.character(1:24) # set colnames to numbers for nice plotting

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
# we bootstrap our model
# => !! it takes rather long (~2h per analysis on a PC with 16 GB RAM & 6 CPUs ~ 2.2 GHz) to run
edge_boot <-
  bootnet(graph_all, nBoots = 100, nCores = 6) # parallelization to multiple cores enabled


# => !! it takes rather long (~2h per analysis on a PC with 16 GB RAM & 6 CPUs ~ 2.2 GHz) to run
case_boot <-
  bootnet(
    graph_all,
    statistics = "edge",
    # we are not interested in centrality measures
    nBoots = 1000,
    type = "case",
    caseN = 10,
    nCores = 6
  ) # parallelization to multiple cores enabled

# -------------------------------- 5: Plot bootstrap results -------------------------
# ---------- Supplementary figure 1 --------------
# 95% CI for edge weights obtained from bootstrapping

# NOTE:
# the following code is an adapation of the "plot.bootnet" function in the bootnet package
# it allows to split the plot depicting bootstrap confidence intervals for edges
# into adjacent sub-plots (assuming split0=T)
# particularly for graphs with many nodes, this ensures readability and visability


x <- edge_boot # the bootnet object we want to plot
which_nodes <- 1:2 # the nodes, here given as integers
n_facets <-
  2 # how many "blocks" do we want the plot to be divided in

# some other settings, we leave them mostly at default
meanVar <- "mean"

minArea <- "q2.5"
maxArea <- "q97.5"

minArea <- paste0(minArea, "_non0")
maxArea <- paste0(maxArea, "_non0")
meanVar <- paste0(meanVar, "_non0")
statistics = "edge"

prop0_minAlpha = 0.25
prop0_cex = 0.75
prop0_alpha = 0.8
rank = F
meanColor <- bootColor <- "black"
sampleColor = "red"


# Compute summary stats:
sumTable <-
  summary(x, statistics = statistics, rank = rank)  %>% ungroup %>% dplyr::mutate_(type = ~
                                                                                     factor(type, levels = statistics))

summary <-
  sumTable %>% dplyr::group_by_(~ id) %>% dplyr::summarize_(sample = ~
                                                              sample[type == statistics[[1]]], mean = as.formula(paste0("~mean(", meanVar, ",na.rm=TRUE)")))

summary$order <- order(order(summary$sample, summary$mean))

sumTable <-
  sumTable %>% dplyr::left_join(dplyr::select_(summary,  ~ id,  ~ order), by = "id")

# Reorder:
sumTable <- sumTable %>%
  dplyr::arrange_(~ dplyr::row_number(order))  %>%
  dplyr::mutate_(id = ~ gsub("^(E|N): ", "", as.character(id))) %>%
  dplyr::mutate_(id = ~ factor(id, levels = unique(id)))


# Some fancy transformation:
revTable <- function(x)
  x[nrow(x):1,]

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

gathered_sumTable <-
  tidyr::gather_(sumTable, "var", "value", c("sample", meanVar))


gathered_sumTable$alpha <-
  ifelse(
    gathered_sumTable$var == meanVar,
    prop0_minAlpha + (1 - prop0_minAlpha) * (1 -
                                               gathered_sumTable$prop0),
    1
  )
sumTable2$alpha <-
  prop0_minAlpha + (1 - prop0_minAlpha) * (1 - sumTable2$prop0)



# for convenience, we determine the facet membership here
split_plot <- gathered_sumTable %>% arrange(node1, node2) %>%
  filter(node1 %in% which_nodes | node2 %in% which_nodes) %>%
  transmute(id, split_plot = -ntile(mean, n_facets)) %>%
  distinct(., id, .keep_all = TRUE)

# the acutal plot comes here
pdf("Supplementary_FigureS1.pdf",
    height = 8,
    width =  8)
gathered_sumTable %>% arrange(node1, node2) %>%
  filter(node1 %in% which_nodes | node2 %in% which_nodes) %>%
  left_join(split_plot, by = "id") %>%
  ggplot(., aes_string(
    x = 'value',
    y = 'id',
    group = 'id',
    colour = "var"
  )) +
  facet_wrap(~ split_plot, scales = "free_y", nrow = 1) +
  
  geom_point(aes_string(alpha = 'alpha'), size = 3) +
  geom_path(
    aes_string(
      y = 'id',
      group = 'id',
      x = 'ci',
      alpha = 'alpha'
    ),
    colour = bootColor,
    data = sumTable2 %>%  arrange(node1, node2) %>% filter(node1 %in% which_nodes |
                                                             node2 %in% which_nodes) %>% left_join(split_plot, by = "id")
  ) +
  geom_label(
    aes(x = 0, label = format(round(prop0, 2), nsmall = 2)),
    cex = prop0_cex * 3.5,
    data = sumTable %>% arrange(node1, node2) %>% filter(node1 %in% which_nodes |
                                                           node2 %in% which_nodes)  %>% left_join(split_plot, by = "id"),
    
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
    values = c(meanColor, sampleColor),
    labels = c("Bootstrap mean", "Sample")
  ) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  scale_alpha(guide = "none")
dev.off()


# ---------- Supplementary figure 2 --------------
# plot a network in which only stable connections (>50% of bootstraps) are contained
main_network <- qgraph(graph_all$graph, DoNotPlot = TRUE)

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


# use this information to plot the network and save it as a pdf in wd ("Supplementary_FigureS2.pdf")
# first, determine layout
all_lay <- qgraph(
  graph_all$graph[2:24, 2:24],
  layout = "spring",
  repulsion = 0.99,
  theme = "colorblind"
)$layout

# manually place age of onset & cumulative use in center of network
all_lay[1, ] <- c(0.1, -0.5)
lay <- rbind(c(0.1, -0.2), all_lay)


# we unfade edges connected to cannabis onset age (1) and cumulative use (2)
fade <- graph_all$graph < 1

fade[2:ncol(graph_all$graph),] <-
  fade[, 2:ncol(graph_all$graph)] <- TRUE
fade[1,] <- fade[, 1] <- fade[2,] <- fade[, 2]  <- FALSE

# acutal plotting & export comes here
pdf(
  colormodel = "cmyk",
  width = 7.0,
  height = 5,
  file = "Supplementary_FigureS2.pdf"
)
qgraph(
  only_stable_edges,
  directed = F,
  layout = lay * -1,
  # flip everything because it looks nicer
  fade = ifelse(only_stable_edges$from == 1, FALSE, TRUE),
  trans = TRUE,
  color =
    c("#BEDEC3", "#B9D1FF", "#FFF9F9", "#DADADA"),
  # color version
  # c("#b3b3b3", "#d9d9d9", "#f0f0f0", "#FFFFFF"), # greyscale version
  groups = c(
    rep("Cannabis Use Characteristics", 2),
    rep("Early Risk Factors", 3),
    rep("Mood", 6),
    rep("Psychosis", 13)
  ),
  theme = "colorblind",
  legend = T,
  cut = 0,
  labels = 1:ncol(data_network),
  GLratio = 1.75,
  label.cex = 1.9,
  vsize = 3,
  legend.cex = 0.35,
  edge.width = 0.5,
  nodeNames = variable_names,
  # edge.color = "black", # greyscale version only
  negDashed = TRUE,
  mar = c(1, 1, 1, 1),
  layoutScale = c(1.12, 1),
  layoutOffset = c(0.1, 0)
  
)
dev.off()

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
pdf("Supplementary_FigureS3.pdf",
    width = 7,
    height = 5)
plot(case_boot, statistics = "edge")
dev.off()

corStability(case_boot) # 0.672 for edge