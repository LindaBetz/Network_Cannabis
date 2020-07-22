# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     Code for    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                                                     #
#                               A network approach to relationships                                   #
#                      between cannabis use characteristics and psychopathology                       #
#                                       in the general population                                     #
#                                                                                                     #
#                                     developed by L. Betz                                            #
#                                                                                                     #
#                              - Analysis reported in main manuscript -                               #
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


# get descriptives for sample
data %>%  rename_at(vars(-CASEID), ~ paste0(variable_names)) %>%
  left_join(X06693_0001_Data[c("CASEID", "V12", "V13")] %>% mutate_all(as.numeric), by =
              "CASEID") %>%
  rename(age = V12,
         sex = V13) %>%
  mutate(sex = ifelse(sex == 1, 0, 1)) %>%
  summarise_all(c("mean", "sd")) %>%
  mutate_all( ~ round(., 3))

# how many % of participants were excluded due to missing values?
(2624 - nrow(data)) / 2624 # 0.0304878

data_network <-
  data %>% select(-CASEID) # drop participant ID for network analysis
colnames(data_network) <-
  as.character(1:24) # set colnames to numbers for nice plotting

# ---------------------------------- 3: Network Estimation -----------------------------------
# we use a mixed graphical model as implemented in mgm-package
graph_all <- estimateNetwork(
  data_network,
  default = "mgm",
  criterion = "EBIC",
  tuning = 0,
  rule = "OR"
)

# ---------------------------------- 4: Plotting the Network (Figure 1) -----------------------------------
# compute layout with all variables except age of onset
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


# here, we actually plot the network and save it as a pdf in wd ("Figure1.pdf")

pdf(
  colormodel = "cmyk",
  width = 7.0,
  height = 5,
  file = "Figure1.pdf"
)
main_network <- qgraph(
  graph_all$graph,
  layout = lay * -1,
  # flip everything because it looks nicer
  fade = fade,
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
  mar = c(4, 4, 4, 4),
  layoutScale = c(1.12, 1),
  layoutOffset = c(0, 0)
  
)
dev.off()