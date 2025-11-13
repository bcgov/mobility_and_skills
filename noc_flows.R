library(tidyverse)
library(here)
library(vroom)
library(janitor)
library(readxl)
library(transport)
library(ggalluvial)
library(conflicted)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
#constants-----------------------------
#set.seed(123)
n_sims <- 1000000
#functions--------------------------------
alluvial_plot <- function(tbbl, initial_age, subsequent_age, cost){
  tbbl|>
    ggplot(aes(axis1 = from, axis2 = to, y = mass)) +
    geom_alluvium(aes(fill = from)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = paste(after_stat(stratum)))) +
    scale_x_discrete(limits = str_replace_all(c(initial_age, subsequent_age), "_"," "))+
    theme_minimal()+
    labs(title=paste0("Lowest cost TEER flows between adjacent age brackets: min(cost)=", round(cost,1)),
         x="Age bracket",
         y=NULL,
         fill="Initial TEER"
    )+
    theme(
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
}
fix_labels <- function(tbbl, origin, destination){
  tbbl|>
    mutate(from=names(origin)[from],
           to=names(destination)[to]
    )
}
transitions_to_matrix <- function(transitions){
  transitions|>
    pivot_wider(names_from = to,  values_from = mass, values_fill = 0)|>
    arrange(from)|>
    column_to_rownames("from")|>
    as.matrix()
}
get_cost <- function(dist, trans){
  sum(dist*trans) #sum of the Hadamard product (element wise matrix multiplication)
}
char_to_int <- function(tbbl){
  tbbl |>
    mutate(from = as.integer(from),
           to   = as.integer(to)) |>
    as_tibble()
}
next_value <- function(from, split_transitions_stage){
  df <- split_transitions_stage[[as.character(from)]]
  if(nrow(df) == 1) return(df$to)
  sample(df$to, size = 1, prob = df$mass)
}

#skill distances----------------------------
skills <- read_excel(here("data","skills_data_for_career_profiles_2025-06-09.xlsx"))|>
  clean_names()|>
  mutate(score=sqrt(level_score*importance_score),
         noc2021=str_pad(noc2021, width = 5, pad="0"),
         teer=str_sub(noc2021, 2,2))|>
  select(noc_5=noc2021,
         teer,
         skill=skills_competencies,
         score)

skill_teer_dist <- skills|>
  group_by(teer, skill)|>
  summarise(score=mean(score))|>
  pivot_wider(names_from = skill, values_from = score)|>
  column_to_rownames("teer")|>
  dist()|>
  as.matrix()

#2011 census data-------------------------
census2011 <- vroom(here("data", "99-012-X2011033.csv"))|>
  pivot_longer(cols=-`Occupation - Na`, names_to = "age_2011", values_to = "count")|>
  rename(noc_2016=`Occupation - Na`)|>
  filter(str_detect(noc_2016, "\\b\\d{4}\\b"),
         !str_detect(noc_2016, "Total"))|>
  mutate(noc_2016=str_sub(noc_2016, 1, 4))
#need to map to noc_5------------------------------------
mapping <- vroom(here("data", "noc2016v1_3-noc2021v1_0-eng.csv" ))|>
  clean_names()|>
  group_by(noc_2016=noc_2016_v1_3_code)|>
  mutate(weight=1/n(), .after=noc_2016_v1_3_code)|> #multiple noc_2016 indicates a split: no info on weights, assume even.
  select(noc_2016, weight, noc_5=noc_2021_v1_0_code)

noc_5_2011 <- full_join(census2011, mapping, by="noc_2016")|>
  mutate(weighted_value=weight*count)|>
  group_by(noc_5, age_2011)|>
  summarize(count=sum(weighted_value))|>
  mutate(noc_5=if_else(noc_5 %in% c("00011","00012","00013","00014","00015"),"00018",noc_5),
         age_2011=str_replace_all(age_2011," ","_")
         )|>
  group_by(noc_5, age_2011)|>
  summarize(count=sum(count))

length(unique(noc_5_2011$noc_5)) #512

#aggregate to teer to begin with--------------------
teer_2011 <- noc_5_2011|>
  mutate(teer=str_sub(noc_5, 2,2))|>
  group_by(age_2011, teer)|>
  summarize(count=sum(count))|>
  group_by(age_2011)|>
  mutate(prop=count/sum(count))|>
  select(-count)|>
  group_by(age_2011)|>
  nest(props_2011=-age_2011)|>
  head(-1)

#2021 census--------------------------------
noc_5_2021<- vroom(here("data", "9810059301.csv"))|>
  clean_names()|>
  rename(noc_desc=contains("noc"),
         age_2021=contains("age"))|>
  mutate(noc_5=str_sub(noc_desc, 1, 5),
         desc=str_sub(noc_desc, 7),
         teer=str_sub(noc_5, 2,2),
         age_2021=str_replace_all(age_2021," ","_"))|>
  select(age_2021,
         noc_5,
         teer,
         desc,
         count=value)

length(unique(noc_5_2021$noc_5)) #512

#aggregate to teer to begin with---------------------------
teer_2021 <- noc_5_2021|>
  mutate(teer=str_sub(noc_5, 2,2))|>
  group_by(age_2021, teer)|>
  summarize(count=sum(count))|>
  group_by(age_2021)|>
  mutate(prop=count/sum(count))|>
  select(-count)|>
  group_by(age_2021)|>
  nest(props_2021=-age_2021)|>
  tail(-1)

#teer: minimize earth movers distance-----------------------
by_teer <- bind_cols(teer_2011, teer_2021)|>
  mutate(distance=list(skill_teer_dist),
         props_2011 = map(props_2011, deframe),
         props_2021 = map(props_2021, deframe),
         transitions=pmap(list(a=props_2011,
                               b=props_2021,
                               costm=distance),
                          transport,
                          method="shortsimplex"),
         transitions=pmap(list(transitions,
                               props_2011,
                               props_2021),
                          fix_labels),
         transition_matrix=map(transitions, transitions_to_matrix),
         cost=map2_dbl(distance,transition_matrix, get_cost),
         plot=pmap(list(transitions, age_2011, age_2021, cost), alluvial_plot)
  )

#save plots to pdf file--------------------------
pdf(here("out","alluvial.pdf"), onefile = TRUE, height=8.5, width=11)
by_teer%>%
  select(plot)%>%
  walk(print)
dev.off()

#heatmap not great for TEER
plt <- by_teer$transitions[by_teer$age_2011=="15_to_24_years"][[1]]|>
  mutate(`Origin TEER`=fct_rev(from),
         `Destination TEER`=fct_rev(to))|>
  ggplot(aes(`Origin TEER`, `Destination TEER`, fill=mass))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()

plotly::ggplotly(plt)

#simulation

#Pre-split each transition table by 'from' (instead of function next_value filtering repeatedly(slow)
split_transitions <- by_teer$transitions |>
  lapply(char_to_int) |>
  lapply(function(tbbl) split(tbbl, tbbl$from))

initial_probs <- enframe(by_teer$props_2011[by_teer$age_2011=="15_to_24_years"][[1]])|>
  mutate(name = as.integer(name))

stage_1 <- sample(initial_probs$name, size = n_sims, replace = TRUE, prob = initial_probs$value)
stage_2 <- vapply(stage_1, function(f) next_value(f, split_transitions[[1]]), integer(1))
stage_3 <- vapply(stage_2, function(f) next_value(f, split_transitions[[2]]), integer(1))
stage_4 <- vapply(stage_3, function(f) next_value(f, split_transitions[[3]]), integer(1))
stage_5 <- vapply(stage_4, function(f) next_value(f, split_transitions[[4]]), integer(1))

# --- 6. Combine into a tibble ---
sims <- tibble(sim = 1:n_sims,
               stage_1, stage_2, stage_3, stage_4, stage_5)

#  empirical proportions ---
empirical_probs <- sims %>%
  summarise(across(stage_1:stage_5, ~ list(prop.table(table(.)))))

#In simulation random sampling noise leads to over-representation in absorbing states.

print(empirical_probs$stage_1)
print(by_teer$props_2011[[1]])

print(empirical_probs$stage_2)
print(by_teer$props_2021[[1]])

print(empirical_probs$stage_3)
print(by_teer$props_2021[[2]])

print(empirical_probs$stage_4)
print(by_teer$props_2021[[3]])

print(empirical_probs$stage_5)
print(by_teer$props_2021[[4]])








