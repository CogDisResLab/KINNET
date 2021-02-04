## code to prepare `stk_annotation` dataset goes here

library(tidyverse)
library(devtools)

# Loading Annotation ------------------------------------------------------

annotation <- read_csv("data-raw/data/RegPhos_Annotation_human.csv") %>%
  mutate(Kinase = str_to_upper(Kinase),
         Group = str_to_upper(Group),
         Family = str_to_upper(Family),
         Subfamily = str_to_upper(Subfamily),
         Gene_Symbol = str_to_upper(Gene_Symbol))

non_human <- c("CK2A1-RS", "SGK424", "STLK6-RS")

annotation_clean <- annotation %>%
  select(Kinase, Group, Family, Subfamily, Gene_Symbol) %>%
  filter(!Kinase %in% non_human) %>%
  mutate(Gene_Symbol = if_else(is.na(Gene_Symbol), Kinase, Gene_Symbol))


# Processing UKA Mapping --------------------------------------------------

uka_map <- read_csv("data-raw/data/uka_pep2kinase_STK.csv") %>%
  select(ID, Kinase_UniprotName) %>%
  rename(Kinase = Kinase_UniprotName) %>%
  filter(!is.na(Kinase)) %>%
  unique

uka_matched_by_kinase <- uka_map %>%
  left_join(annotation_clean, by = "Kinase") %>%
  filter(!is.na(Family), !is.na(Group))

uka_matched_by_family <- uka_map %>%
  rename(Family = Kinase) %>%
  left_join(annotation_clean, by = "Family") %>%
  filter(!is.na(Kinase), !is.na(Group))

uka_matched_by_subfamily <- uka_map %>%
  rename(Subfamily = Kinase) %>%
  left_join(annotation_clean, by = "Subfamily") %>%
  filter(!is.na(Kinase), !is.na(Group), !is.na(Family))

uka_matched_by_gene <- uka_map %>%
  rename(Gene_Symbol = Kinase) %>%
  left_join(annotation_clean, by = "Gene_Symbol") %>%
  filter(!is.na(Kinase), !is.na(Group), !is.na(Family))




# Processing KRSA Mapping -------------------------------------------------

krsa_map <- read_csv("data-raw/data/KRSA_Mapping_STK_Separated.csv")

krsa_matched_by_family <- krsa_map %>%
  left_join(annotation_clean, by = "Family") %>%
  filter(!is.na(Kinase), !is.na(Group))

krsa_matched_by_kinase <- krsa_map %>%
  rename(Kinase = Family) %>%
  left_join(annotation_clean, by = "Kinase") %>%
  filter(!is.na(Family), !is.na(Group))

krsa_matched_by_subfamily <- krsa_map %>%
  rename(Subfamily = Family) %>%
  left_join(annotation_clean, by = "Subfamily") %>%
  filter(!is.na(Kinase), !is.na(Group), !is.na(Family))

krsa_matched_by_gene <- krsa_map %>%
  rename(Gene_Symbol = Family) %>%
  left_join(annotation_clean, by = "Gene_Symbol") %>%
  filter(!is.na(Kinase), !is.na(Group), !is.na(Family))


# Processing Special Cases -----------------------------------------------


# ERK, PAKA, MTOR, PAKB are special cases from KRSA

krsa_mtor <- annotation %>%
  filter(str_detect(UniProt_ID, "MTOR")) %>%
  select(Kinase, Group, Family, Subfamily, Gene_Symbol) %>%
  mutate(Match = "MTOR") %>%
  left_join(krsa_map, by = c(Match = "Family")) %>%
  select(-Match)

krsa_paka <- annotation %>%
  filter(str_detect(Subfamily, "PAKA")) %>%
  select(Kinase, Group, Family, Subfamily, Gene_Symbol) %>%
  mutate(Match = "PAKA") %>%
  left_join(krsa_map, by = c(Match = "Family")) %>%
  select(-Match)

krsa_pakb <- annotation %>%
  filter(str_detect(Subfamily, "PAKB")) %>%
  select(Kinase, Group, Family, Subfamily, Gene_Symbol) %>%
  mutate(Match = "PAKB") %>%
  left_join(krsa_map, by = c(Match = "Family")) %>%
  select(-Match)

krsa_erks <- map_dfr(
  c(1:5,7),
  function(x) {
    annotation %>%
      filter(Kinase == str_c("ERK", 1, sep = "")) %>%
      select(Kinase, Group, Family, Subfamily, Gene_Symbol) %>%
      mutate(Match = "ERK") %>%
      left_join(krsa_map, by = c(Match = "Family")) %>%
      select(-Match)
  }
)


# Putting it all together -------------------------------------------------

stk_annotation <- bind_rows(
  krsa_matched_by_family,
  krsa_matched_by_kinase,
  krsa_matched_by_subfamily,
  krsa_matched_by_gene,
  krsa_erks,
  krsa_mtor,
  krsa_paka,
  krsa_pakb,
  uka_matched_by_family,
  uka_matched_by_gene,
  uka_matched_by_kinase,
  uka_matched_by_subfamily
) %>%
  select(ID, Group, Family, Subfamily, Kinase, Gene_Symbol) %>%
  unique()
