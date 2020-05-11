library(haven)
library(qgraph)
library(bootnet)
library(psychonetrics)
library(dplyr)


#X06693_0001_Data <- read_sav("DS0001/06693-0001-Data.sav")
#X06693_0002_Data <- read_sav("DS0002/06693-0002-Data.sav")

variable_names <- c(
  "age of onset",
  "cumulative use",
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
  "HYPNOTZ/MAGIC/FORCE",
  "SAW VISIONS",
  "HEARD NOISE/VOICE",
  "SMELLS/BODY ODORS",
  "FEELNG IN/ON BODY"
) %>% tolower(.)


data <- X06693_0001_Data %>%
  filter(V4097 <= 40) %>%
  transmute(
    #V1907,
    age_onset = V1908,
    # 1
    V1909,
    # 2
    
    mol = ifelse(V6126 == 5, 0, ifelse(V6127 < 18, 1, 0)),
    
    raped = ifelse(V6114 == 5, 0, ifelse(V6115 < 18, 1, 0)),
    physical_assault = ifelse(V6138 == 5, 0, ifelse(V6139 < 18, 1, 0)),
    sexual_abuse = ifelse(mol == 1 | raped == 1, 1, 0),
    
    
    V6143 = V6143,
    # sexual abuse
    
    abuse = ifelse(sexual_abuse == 1 |
                     V6143 == 1 | physical_assault == 1, 1, 0),
    
    V6144,
    # neglect
    
    V7136,
    # urbanicity
    V301,
    # panic: Have you ever in your life had a spell or attack when  all of a sudden you felt frightened, anxious or very uneasy in situations when most people would not be afraid or anxious?
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
  select(-c(mol, raped, sexual_abuse, V6143, physical_assault)) %>%
  filter(V4101 != 3) %>% # section K skipped entirely
  na.omit() %>%
  mutate_all(as.numeric) %>%
  mutate_at(vars(matches("V41|V3|V7403|V6|V1907")), ~ recode(.,
                                                             `5` = 0)) %>%
  mutate(V7136 = ifelse(V7136 == 5 | V7136 == 4, 1, 0))



graph_all <- estimateNetwork(
  data,
  # whole sample
  default = "mgm",
  criterion = "EBIC",
  tuning = 0,
  rule = "OR"
)




all_lay <- qgraph(
  graph_all$graph[2:24, 2:24],
  layout = "spring",
  repulsion = 0.9,
  theme = "colorblind"
)$layout


lay <- rbind(c(0, 0.4), all_lay)



fade <- graph_all$graph > 0

fade[2:ncol(graph_all$graph), ] <- TRUE
fade[, 2:ncol(graph_all$graph)] <- TRUE

fade[1, ] <- FALSE
fade[, 1] <- FALSE
fade[2, ] <- FALSE
fade[, 2] <- FALSE

#layout(t(matrix(c(1, 2,3), byrow = F)), widths = c(2, 2, 4))
all <- qgraph(
  graph_all$graph,
  layout = lay,
  fade = fade,
  trans = FALSE,
  color = c("#FED439", "#D2AF81", "#00A1D5FF", "#cad3db"),
  groups = c(
    rep("Cannabis Use", 2),
    rep("Early Risk", 3),
    rep("Mood", 6),
    rep("Psychosis", 13)
  ),
  theme = "colorblind",
  legend = T,
  minimum = 0.01,
  maximum = 0.8,
  cut = 0.6,
  #egend.mode = "groups",
  labels = 1:ncol(data),
  GLratio = 1.75,
  label.cex = 1.9,
  vsize = 3,
  legend.cex = 0.5,
  filename="out",
  filetype="pdf",
  edge.width = 1,
  nodeNames = variable_names
)

## bootstrapping procedures
set.seed(1)
case_boot_all <- bootnet(graph_all, nBoots = 1000, type="case", caseN = 10, nCores = 6)
corStability(case_boot_all)

set.seed(1)
edge_boot_all <- bootnet(graph_all, nBoots = 1000, nCores = 6)


compare_early_late <- NetworkComparisonTest::NCT(graph_early, graph_late, it=500)
