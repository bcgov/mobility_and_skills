library(tidyverse)
library(here)
library(vroom)
library(janitor)
library(readxl)
library(transport)
library(plotly)
library(scales)
library(grid)
library(ggalluvial)
library(conflicted)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
options(scipen = 999)
#functions------------------------------
fix_labels <- function(tbbl, origin, destination){
  tbbl|>
    mutate(from=names(origin)[from],
           to=names(destination)[to]
    )
}
read_data <- function(file_name){
  read_excel(here("data", "onet", file_name))%>%
    clean_names()%>%
    select(o_net_soc_code, element_name, scale_name, data_value)%>%
    pivot_wider(names_from = scale_name, values_from = data_value)%>%
    mutate(score=sqrt(Importance*Level), #geometric mean of importance and level
           #mutate(score=Level,
           category=(str_split(file_name,"\\.")[[1]][1]))%>%
    unite(element_name, category, element_name, sep=": ")%>%
    select(-Importance, -Level)
}
get_cost <- function(tbbl){
  tbbl|>
    left_join(census_2021$all_noc_names, by=c("from"="noc_5"))|>
    rename(from_name=desc)|>
    left_join(census_2021$all_noc_names, by=c("to"="noc_5"))|>
    rename(to_name=desc)|>
    mutate(cost=net_mass*distance)|>
    unite("to", "to", "to_name", sep=": ")|>
    unite("from", "from", "from_name", sep=": ")
}

make_segment_data <- function(transitions, coordinates){
    transitions|>
    left_join(coordinates, by = c("from" = "noc_5")) %>%
    rename(x_from = V1, y_from = V2) %>%
    left_join(coordinates, by = c("to" = "noc_5")) %>%
    rename(x_to = V1, y_to = V2) %>%
    filter(from != to) %>%
    filter(!(x_from == x_to & y_from == y_to))
}

net_segments <- function(transitions){
  transitions %>%
    filter(from != to)%>%
    mutate(
      pair1 = pmax(from, to),
      pair2 = pmin(from, to)
    )%>%
    group_by(pair1, pair2) %>%
    summarise(
      net_mass = abs(sum(if_else(from == pair1, mass, -mass))),
      .groups = "drop"
    )|>
    rename(from=pair1, to=pair2)
}

convert_to_teer <- function(transitions){
  transitions|>
    mutate(to=str_sub(to,2,2),
           from=str_sub(from,2,2))|>
    group_by(to, from)|>
    summarize(mass=sum(mass))
}

network_plot <- function(coordinates, segment_data, top_segments, age_2011, age_2021){
  ggplot()+
    geom_point(data=coordinates, mapping=aes(V1,V2), size=.25, alpha=.25)+
    geom_curve(data=segment_data, mapping=aes(x=x_from, y=y_from, xend=x_to, yend=y_to),
               linewidth=.2,
               alpha=.2,
               curvature = .15,
               colour = "steelblue")+
    geom_curve(data=top_segments, mapping=aes(x=x_from, y=y_from, xend=x_to, yend=y_to, colour=segment_name),
               curvature = .15,
               linewidth = .75,
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
    )+
    scale_color_brewer(palette = "Dark2")+
    theme_void()+
    labs(title="Solution to Optimal Transport Problem",
         subtitle=paste0("...for transition between ",
                         str_replace_all(age_2011, "_", " "),
                         " (2011) and ",
                         str_replace_all(age_2021, "_", " "),
                         " (2021)"),
         colour="Top 8 transitions")
}

source_dest_plot <- function(source, dest, age_2011, age_2021, cut_off){
  full_join(source, dest, by = join_by(noc_5, desc))|>
    mutate(diff=prop.y-prop.x,
           fill=if_else(diff<0, "Decreases", "Increases"))|>
    unite(noc_5, noc_5, desc, sep = ": ")|>
    filter(abs(diff)>cut_off)|>
    ggplot(aes(diff, fct_reorder(noc_5, diff, .desc = TRUE), fill=fill))+
    geom_col()+
    theme_minimal()+
    scale_fill_brewer(palette = "Dark2")+
    labs(title=paste0("Changes in Census NOC proportions larger than ",
                      cut_off),
         subtitle=paste0("... for transition between ",
                         str_replace_all(age_2011, "_", " "),
                         " (2011) and ",
                         str_replace_all(age_2021, "_", " "),
                         " (2021)"),
         x="Difference in Proportion",
         y=NULL,
         fill=NULL)+
    coord_cartesian(xlim = c(-.1,.1))
}

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


#2021 census--------------------------------

census_2021 <- list() #store all census 2021 stuff together

census_2021$implicit_missing<- vroom(here("data", "9810059301.csv"))|>
  clean_names()|>
  rename(noc_desc=contains("noc"),
         age_2021=contains("age"))|>
  mutate(noc_5=str_sub(noc_desc, 1, 5),
         desc=str_sub(noc_desc, 7),
         age_2021=str_replace_all(age_2021," ","_"))|>
  filter(age_2021!="15_to_24_years")|> #do not need in 2021
  select(age_2021,
         noc_5,
         desc,
         count=value)

#get the full grid of nocs and ages

census_2021$all_noc_names <- census_2021$implicit_missing|>
  select(noc_5, desc)|>
  unique() #should be length 512

census_2021$noc_5 <- unique(census_2021$implicit_missing$noc_5)
census_2021$age_2021 <- unique(census_2021$implicit_missing$age_2021)
census_2021$all_2021_combos <- crossing(noc_5=census_2021$noc_5, age_2021=census_2021$age_2021)

#join and replace nas with 0s

census_2021$census_2021 <- left_join(census_2021$all_2021_combos, census_2021$implicit_missing)|>
  mutate(count=if_else(is.na(count), 0 , count))|>
  group_by(age_2021)|>
  filter(age_2021!="15_to_24_years")|> #they were babies in 2011.
  mutate(prop=count/sum(count))|>
  select(-count)|>
  group_by(age_2021)|>
  nest(props_2021=-age_2021)

#2011 census-----------------------------------------------------

census_2011 <- list() #store all census 2021 stuff together

census_2011$census_2011_raw <- vroom(here("data", "99-012-X2011033.csv"))|>
  pivot_longer(cols=-`Occupation - Na`, names_to = "age_2011", values_to = "count")|>
  rename(noc_2016=`Occupation - Na`)|>
  filter(str_detect(noc_2016, "\\b\\d{4}\\b"),
         !str_detect(noc_2016, "Total"))|>
  mutate(noc_2016=str_sub(noc_2016, 1, 4))

#need to map to noc_5

census_2011$mapping <- vroom(here("data", "mapping","noc2016v1_3-noc2021v1_0-eng.csv" ))|>
  clean_names()|>
  group_by(noc_2016=noc_2016_v1_3_code)|>
  mutate(weight=1/n(), .after=noc_2016_v1_3_code)|> #multiple noc_2016 indicates a split: no info on weights, assume even.
  select(noc_2016, weight, noc_5=noc_2021_v1_0_code, desc=noc_2021_v1_0_title)

census_2011$census_2011 <- full_join(census_2011$census_2011_raw, census_2011$mapping, by="noc_2016")|>
  mutate(weighted_value=weight*count)|>
  group_by(noc_5, desc, age_2011)|>
  summarize(count=sum(weighted_value))|>
  mutate(noc_5=if_else(noc_5 %in% c("00011","00012","00013","00014","00015"),"00018",noc_5),
         desc=if_else(noc_5=="00018","Seniors managers - public and private sector" ,desc),
         age_2011=str_replace_all(age_2011," ","_"))|>
  group_by(noc_5, desc, age_2011)|>
  summarize(count=sum(count))|>
  filter(age_2011!="55_to_64_years")|>
  group_by(age_2011)|>
  mutate(prop=count/sum(count))|>
  select(-count)|>
  group_by(age_2011)|>
  nest(props_2011=-age_2011)

#get skills data-----------------------

skills <- list() #store all skill stuff together

skills$mapping <- read_excel(here("data","mapping", "onet2019_soc2018_noc2016_noc2021_crosswalk_consolodated.xlsx"))%>%
  mutate(noc_5=str_pad(noc2021, "left", pad="0", width=5))%>%
  select(noc_5, o_net_soc_code = onetsoc2019)%>%
  distinct()

#the onet data
skills$onet_raw <- tibble(file=c("Skills.xlsx", "Abilities.xlsx", "Knowledge.xlsx", "Work Activities.xlsx"))%>%
  mutate(data=map(file, read_data))%>%
  select(-file)%>%
  unnest(data)%>%
  pivot_wider(id_cols = o_net_soc_code, names_from = element_name, values_from = score)%>%
  inner_join(skills$mapping)%>%
  ungroup()%>%
  select(-o_net_soc_code)%>%
  select(noc_5, everything())%>%
  group_by(noc_5)%>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))%>% #mapping not one to one: mean gives one value per NOC
  mutate(across(where(is.numeric), ~ if_else(is.na(.), mean(., na.rm=TRUE), .))) #11 nocs 4 measures replace na with mean

skills$all_two_digit_skills <- skills$onet_raw|>
  mutate(noc_2=str_sub(noc_5,1,2))|>
  group_by(noc_2)|>
  summarise(across(contains(":"), ~mean(.x, na.rm = TRUE)))

skills$two_digit_skills <- anti_join(census_2021$all_noc_names, skills$onet_raw|>select(noc_5))|>
  select(noc_5)|>
  mutate(noc_2=str_sub(noc_5, 1,2))|>
  inner_join(skills$all_two_digit_skills)|>
  select(-noc_2)

skills$onet_full <- bind_rows(skills$onet_raw, skills$two_digit_skills)|>
  column_to_rownames("noc_5")

skills$onet_pca <- prcomp(skills$onet_full, center=TRUE, scale=TRUE)
skills$onet_scores <- skills$onet_pca$x[, 1:10]#keep first 10 components
skills$skills_noc_dist<- dist(skills$onet_scores, method = "euclidean")|>
  as.matrix()

skills$mds2 <- cmdscale(skills$skills_noc_dist, k = 2)|>
  as.data.frame()|>
  rownames_to_column("noc_5")|>
  left_join(census_2021$all_noc_names)|>
  unite(noc_5, noc_5, desc, sep=": ")

skills$distance_long <- skills$skills_noc_dist|>
  as.data.frame()|>
  rownames_to_column("from")|>
  pivot_longer(cols=-from, names_to = "to", values_to = "distance")

#minimize earth movers distance-----------------------

results <- bind_cols(census_2011$census_2011, census_2021$census_2021)|>
  mutate(facet_label=paste(age_2011,age_2021, sep=" -> "),
         distance=list(skills$skills_noc_dist),
         distance_long=list(skills$distance_long),
         props_2011_vec = map(props_2011, \(tbbl) tbbl |> select(-desc) |> deframe()),
         props_2021_vec = map(props_2021, \(tbbl) tbbl |> select(-desc) |> deframe()),
         transitions=pmap(list(a=props_2011_vec,
                               b=props_2021_vec,
                               costm=distance),
                          transport,
                          method="shortsimplex"),
         transitions=pmap(list(transitions,
                               props_2011_vec,
                               props_2021_vec),
                          fix_labels),
         teer_transitions=map(transitions, convert_to_teer),
         net_transitions=map(transitions, net_segments),
         net_transitions=map2(net_transitions, distance_long, left_join),
         net_transitions=map(net_transitions, get_cost),
         total_cost=map_dbl(net_transitions, ~ sum(.x$cost)),
         coordinates=list(skills[["mds2"]]),
         segment_data=map2(net_transitions, coordinates, make_segment_data),
         top_segments=map(segment_data, slice_max, order_by=net_mass, n=8),
         top_segments=map(top_segments, \(tbbl) tbbl |> unite(segment_name, from, to, sep = " -> ")),
         top_segments=map(top_segments, \(tbbl) tbbl |> mutate(segment_name=fct_reorder(segment_name, net_mass, .desc = TRUE))),
         source_destination_plot=pmap(list(props_2011, props_2021, age_2011, age_2021, cut_off=.005), source_dest_plot),
         transition_plot=pmap(list(coordinates, segment_data, top_segments, age_2011, age_2021), network_plot),
         teer_alluvium_plot=pmap(list(teer_transitions, age_2011, age_2021, total_cost), alluvial_plot)
         )



write_rds(skills, "skills.rds")
write_rds(results, "results.rds")





