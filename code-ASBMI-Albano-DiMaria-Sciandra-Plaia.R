## Authors: Alessandro Albano, Chiara di Maria, Mariangela Sciandra and Antonella Plaia

#--Load libraries ------------------------------------------------
library(tidyverse)
library(tidytext)
library(tidylo)
library(tm)
library(grf)


#-- START MERGING MIMIC III --------------------------------------------
# In this phase two csv files from MIMIC III are required:
# 1) NOTEEVENTS to retrieve discharge summaries
# 2) DIAGNOSES_ICD to retrieve the final dignoses

# Read raw note events CSV and filter only "Report" descriptions
raw_notes <- read.csv2("NOTEEVENTS.csv", sep = ",", stringsAsFactors = FALSE)
note_reports <- raw_notes %>%
  filter(DESCRIPTION == "Report") %>%
  select(-DESCRIPTION) %>% 
  tibble()

# Load ICD-9 diagnosis codes
diagnoses_icd <- read_csv("DIAGNOSES_ICD.csv")

# Filter only unique discharge summaries per patient-admission
discharge_counts <- note_reports %>%
  filter(CATEGORY == "Discharge summary") %>%
  count(SUBJECT_ID, HADM_ID) %>%
  filter(n == 1) %>%
  select(SUBJECT_ID, HADM_ID)

unique_notes <- note_reports %>%
  filter(CATEGORY == "Discharge summary")  %>% 
  semi_join(discharge_counts, by = c("SUBJECT_ID", "HADM_ID"))

# Join notes with diagnosis records and clean up helper columns
notes_with_diag <- unique_notes %>%
  inner_join(diagnoses_icd, by = c("SUBJECT_ID", "HADM_ID")) %>%
  select(-ROW_ID.y, -SEQ_NUM) %>%
  rename(icd9_code = ICD9_CODE) %>%
  mutate(icd9_prefix = substr(icd9_code, 1, 3))

#-- Compute binary outcomes for each condition ---------------------------
# Prefix 250 = diabetes, 255 = adrenal disorders, 244 = hypothyroidism
disease_flags <- notes_with_diag %>%
  group_by(SUBJECT_ID, HADM_ID) %>%
  summarise(
    diab    = as.integer(any(icd9_prefix == "250")),
    surren     = as.integer(any(icd9_prefix == "255")),
    hypot = as.integer(any(icd9_prefix == "244")),
    .groups = "drop"
  )

# Build dataset: one row per patient-admission including flags
final_disease_data <- notes_with_diag %>%
  distinct(SUBJECT_ID, HADM_ID, .keep_all = TRUE) %>%
  left_join(disease_flags, by = c("SUBJECT_ID", "HADM_ID")) %>%
  select(-CHARTDATE, -CHARTTIME, -STORETIME, -CGID, -ISERROR) %>% 
  drop_na(diab, surren, hypot)


#-- Remove intermediate objects to clean environment ----------------------

rm(list = setdiff(ls(), c("final_disease_data")))

ls()


#-- END MERGING MIMIC III --------------------------------------------


# -- START PRE PROCESSING ------------------------------------------

#-- Text cleaning -----------------------------------------------
final_disease_data <- final_disease_data %>%
  mutate(
    Id     = row_number(),
    TEXT = TEXT %>%
      str_remove_all("[[:punct:]]") %>%
      str_remove_all("[[:digit:]]")
  )


#-- Document-Term Matrix preparation 

# Token extraction and DTM creation
tokens_dtm <- final_disease_data %>%
  unnest_tokens(word, TEXT) %>%
  anti_join(get_stopwords(), by = "word") %>%
  filter(str_length(word) > 2) %>% 
  count(Id, word) %>%
  cast_dtm(Id, word, n)

# Remove sparse terms
dtm_filtered <- removeSparseTerms(tokens_dtm, 0.95)
stop_words = stop_words[-938,]

# Filter tokens based on sparse-term vocabulary
filtered_tokens <- final_disease_data %>%
  unnest_tokens(word, TEXT) %>%
  anti_join(get_stopwords(), by = "word") %>%
  anti_join(stop_words)  %>%
  filter(str_length(word) > 2) %>%
  filter(word %in% Terms(dtm_filtered)) 


count_filtered_tokens <- filtered_tokens %>%  
  count(Id ,word, sort = TRUE)


DTM_byhand <- count_filtered_tokens |>
  spread(key= word, value=n, fill = 0)|>
  tibble()

# check inconsistencies
disease_data_reduced<- final_disease_data |>
  filter(Id %in% DTM_byhand$Id )



#-- Weighted-Log‑odds keyword analysis -------------------------------------------

#Diabetes
diab_tolog <- final_disease_data %>%
  unnest_tokens(word, TEXT) %>%
  anti_join(get_stopwords()) %>%
  filter(nchar(word)>2) %>%
  filter( word %in% colnames(dtm_filtered)) %>% 
  mutate(
    term = case_when(
      word %in% c("sugar",  "sugars")  ~ "sugar",
      word %in% c("obese",  "obesity") ~ "obesity",
      TRUE                             ~ word
    )
  ) %>%
  count(diab, term, sort = TRUE) %>%
  bind_log_odds(diab, term, n) %>%
  arrange(desc(log_odds_weighted)) %>%
  filter(diab == 1)


diabetes_terms <- c("insulin", "units", "sugar", "subcutaneous", 
                    "renal", "obesity", "dialysis", 
                    "cad", "foot", "hepatitis", "breast", 
                    "screening", "bilirubin", "head", "baby")

# Filter rows containing any diabetes keyword
diab_tolog %>%  filter(term %in% diabetes_terms)



# Hypothyroidism
hyp_tolog <- final_disease_data %>%
  unnest_tokens(word, TEXT) %>%
  anti_join(get_stopwords()) %>%
  filter(nchar(word)>2) %>%
  filter( word %in% colnames(dtm_filtered)) %>% 
  count(hypot, word, sort = TRUE) %>%
  anti_join(get_stopwords()) %>%
  filter(nchar(word)>2) %>%
  bind_log_odds(hypot, word, n) %>%
  arrange(desc(log_odds_weighted)) %>%
  filter(hypot == 1)


hypothyroid_terms <- c("tsh", "chf", "uti", "afib", 
                       "female", "woman", "man", "head", 
                       "cardiovascular", "male", "hepatitis", "baby")

hyp_tolog %>%  filter(word %in% hypothyroid_terms)


# Disorders of Adrenal Glands

sur_tolog <-final_disease_data %>%
  unnest_tokens(word, TEXT) %>%
  anti_join(get_stopwords()) %>%
  filter(nchar(word)>2) %>%
  filter( word %in% colnames(dtm_filtered)) %>% 
  mutate(
    term = case_when(
      word %in% c("steroid",  "steroids") ~ "steroid",
      TRUE                             ~ word
    )
  ) %>%
  count(surren, term, sort = TRUE) %>%
  bind_log_odds(surren, term, n) %>%
  arrange(desc(log_odds_weighted)) |> 
  filter(surren==1) 


adrenal_terms <- c("insufficiency", "steroid", "steroids", 
                   "prednisone", "hypotension", "stress", 
                   "vancomycin", "micu", "aortic", "postoperative", 
                   "valve", "coronary", "artery")


sur_tolog %>%  filter(term %in% adrenal_terms)




# -- END PRE PROCESSING ------------------------------------------


# -- START CAUSAL FOREST MODELS ------------------------------------------
set.seed(1995)
s = sample(1:nrow(DTM_byhand), 0.7*nrow(DTM_byhand))
DTM_byhand_train = DTM_byhand[s,]
DTM_byhand_test = DTM_byhand[-s,]

######################### Models for diabetes #########################


ordered_diab <- arrange(disease_data_reduced,Id)$diab

#1) insulin
mod_1_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$insulin,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="insulin")]
)
# ATE estimate
average_treatment_effect(mod_1_cf)

Y_diab_train <- (arrange(disease_data_reduced, Id)$diab)[s]
Y_diab_test <- (arrange(disease_data_reduced, Id)$diab)[-s]
insubin = 1*( DTM_byhand$insulin>0)


mod_train_cf_ins = causal_forest(
  Y = Y_diab_train,
  W = insubin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)=="insulin")]
)
tau.hat.eval.insul <- predict(mod_train_cf_ins, DTM_byhand_test[, -which(colnames(DTM_byhand_test)=="insulin")])$predictions
eval.forest.insul <- causal_forest(
  Y = Y_diab_test, 
  W = insubin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)=="insulin")] )

# RATE estimate
rate.cate.insul <- rank_average_treatment_effect(eval.forest.insul, tau.hat.eval.insul)

plot(rate.cate.insul, main = "TOC: By decreasing estimated CATE (insulin)")

# Repeat this process for each word to be tested


#2) units ####

mod_2_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$units,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="units")]
)

#ATE ESTIMATE
average_treatment_effect(mod_2_cf)  


uni_bin = 1*( DTM_byhand$units>0)

mod_train_cf_uni = causal_forest(
  Y = Y_diab_train,
  W = uni_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)=="units")]
)

tau.hat.eval.uni <- predict(mod_train_cf_uni, DTM_byhand_test[, -which(colnames(DTM_byhand_test)=="units")])$predictions

eval.forest.uni <- causal_forest(
  Y = Y_diab_test, 
  W = uni_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)=="units")] )

# RATE estimate
rate.cate.uni <- rank_average_treatment_effect(eval.forest.uni, tau.hat.eval.uni)


plot(rate.cate.uni, main = "TOC: By decreasing estimated CATE (units)")


# 3) sugars ###

mod_sug_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$sugars + DTM_byhand$sugar,
  X = DTM_byhand[, -which(colnames(DTM_byhand_train)%in%c("sugars", "sugar"))]
)

#ATE ESTIMATE
average_treatment_effect(mod_sug_cf) 

sug_bin = 1*( DTM_byhand$sugar + DTM_byhand$sugars  >0)

mod_train_cf_sug = causal_forest(
  Y = Y_diab_train,
  W = sug_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in%c("sugars", "sugar"))]
)
tau.hat.eval.sug <- predict(mod_train_cf_sug, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in%c("sugars", "sugar"))])$predictions
eval.forest.sug <- causal_forest(
  Y = Y_diab_test, 
  W = sug_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in%c("sugars", "sugar"))] )

rate.cate.sug <- rank_average_treatment_effect(eval.forest.sug, tau.hat.eval.sug)

plot(rate.cate.sug, main = "TOC: By decreasing estimated CATE (sugar)")


# 4) obesity

mod_obes_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$obese + DTM_byhand$obesity,
  X = DTM_byhand[, -which(colnames(DTM_byhand) %in% c("obese", "obesity"))]
)

# ATE estimate
average_treatment_effect(mod_obes_cf)


obes_bin = 1*( DTM_byhand$obese>0 |  DTM_byhand$obesity>0)
mod_train_cf_obes = causal_forest(
  Y = Y_diab_train,
  W = obes_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('obese', 'obesity'))]
)
tau.hat.eval.obes <- predict(mod_train_cf_obes, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('obese', 'obesity'))])$predictions
eval.forest.obes <- causal_forest(
  Y = Y_diab_test, 
  W = obes_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('obese', 'obesity'))] )

#RATE estimate
rate.cate.obes <- rank_average_treatment_effect(eval.forest.obes, tau.hat.eval.obes)
plot(rate.cate.obes, main = "TOC: By decreasing estimated CATE (obese/obesity)")



# 5) subcutaneous

mod_subc_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$subcutaneous,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="subcutaneous")]
)
# ATE estimate
average_treatment_effect(mod_subc_cf) 

subc_bin = 1*( DTM_byhand$subcutaneous>0)

mod_train_cf_subc = causal_forest(
  Y = Y_diab_train,
  W = subc_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('subcutaneous'))]
)

tau.hat.eval.subc <- predict(mod_train_cf_subc, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('subcutaneous'))])$predictions

eval.forest.subc <- causal_forest(
  Y = Y_diab_test, 
  W = subc_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('subcutaneous'))] )

#RATE estimate
rate.cate.subc <- rank_average_treatment_effect(eval.forest.subc, tau.hat.eval.subc)

plot(rate.cate.subc, main = "TOC: By decreasing estimated CATE (subcutaneous)")


#6) cad

mod_cad_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$cad,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="cad")] 
)
# ATE estimate
average_treatment_effect(mod_cad_cf) 
cad_bin = 1*( DTM_byhand$cad>0)

mod_train_cf_cad = causal_forest(
  Y = Y_diab_train,
  W = cad_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('cad'))]
)

tau.hat.eval.cad <- predict(mod_train_cf_cad, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('cad'))])$predictions
eval.forest.cad <- causal_forest(
  Y = Y_diab_test, 
  W = cad_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('cad'))] )

#RATE estimate
rate.cate.cad <- rank_average_treatment_effect(eval.forest.cad, tau.hat.eval.cad)

plot(rate.cate.cad, main = "TOC: By decreasing estimated CATE (cad)")


#7) dialysis

mod_dia_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$dialysis,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="dialysis")] 
)

#ATE estimate
average_treatment_effect(mod_dia_cf) 

dia_bin = 1*( DTM_byhand$dialysis>0)

mod_train_cf_dia = causal_forest(
  Y = Y_diab_train,
  W = dia_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('dialysis'))]
)
tau.hat.eval.dia <- predict(mod_train_cf_dia, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('dialysis'))])$predictions
eval.forest.dia <- causal_forest(
  Y = Y_diab_test, 
  W = dia_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('dialysis'))] )

#RATE estimate
rate.cate.dia <- rank_average_treatment_effect(eval.forest.dia, tau.hat.eval.dia)


plot(rate.cate.dia, main = "TOC: By decreasing estimated CATE (dialysis)")


#8) renal

mod_ren_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$renal,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="renal")] 
)

#ATE estimate
average_treatment_effect(mod_ren_cf) 

ren_bin = 1*( DTM_byhand$renal>0)

mod_train_cf_ren = causal_forest(
  Y = Y_diab_train,
  W = ren_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('renal'))]
)
tau.hat.eval.ren <- predict(mod_train_cf_ren, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('renal'))])$predictions
eval.forest.ren <- causal_forest(
  Y = Y_diab_test, 
  W = ren_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('renal'))] )

#RATE estimate
rate.cate.ren <- rank_average_treatment_effect(eval.forest.ren, tau.hat.eval.ren)


plot(rate.cate.ren, main = "TOC: By decreasing estimated CATE (renal)")


#9) foot

mod_foot_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$foot,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="foot")] 
)

#ATE estimate
average_treatment_effect(mod_foot_cf) 

foot_bin = 1*( DTM_byhand$foot>0)


mod_train_cf_foot = causal_forest(
  Y = Y_diab_train,
  W = foot_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('foot'))]
)

tau.hat.eval.foot <- predict(mod_train_cf_foot, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('foot'))])$predictions

eval.forest.foot <- causal_forest(
  Y = Y_diab_test, 
  W = foot_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('foot'))] )

# RATE estimate
rate.cate.foot <- rank_average_treatment_effect(eval.forest.foot, tau.hat.eval.foot)
plot(rate.cate.foot, main = "TOC: By decreasing estimated CATE (foot)")



#10) screening

mod_scr_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$screening,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="screening")] 
)
# ATE estimate

scr_bin = 1*( DTM_byhand$screening>0)

mod_train_cf_scr = causal_forest(
  Y = Y_diab_train,
  W = scr_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('screening'))]
)

tau.hat.eval.scr <- predict(mod_train_cf_scr, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('screening'))])$predictions
eval.forest.scr <- causal_forest(
  Y = Y_diab_test, 
  W = scr_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('screening'))] )

#RATE estimate
rate.cate.scr <- rank_average_treatment_effect(eval.forest.scr, tau.hat.eval.scr)

plot(rate.cate.scr, main = "TOC: By decreasing estimated CATE (screening)")



#11) breast

mod_bre_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$breast,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="breast")] 
)
#ATE estimate
average_treatment_effect(mod_bre_cf)

bre_bin = 1*( DTM_byhand$breast>0)

mod_train_cf_bre = causal_forest(
  Y = Y_diab_train,
  W = bre_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('breast'))]
)
tau.hat.eval.bre <- predict(mod_train_cf_bre, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('breast'))])$predictions
eval.forest.bre <- causal_forest(
  Y = Y_diab_test, 
  W = bre_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('breast'))] )

#RATE estimate
rate.cate.breast <- rank_average_treatment_effect(eval.forest.bre, tau.hat.eval.bre)

plot(rate.cate.breast, main = "TOC: By decreasing estimated CATE (breast)")



## 12) hepatitis

mod_hepa_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$hepatitis,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="hepatitis")] 
)
#ATE estimate
average_treatment_effect(mod_hepa_cf)  
hepa_bin = 1*( DTM_byhand$hepatitis>0)

mod_train_cf_hepa = causal_forest(
  Y = Y_diab_train,
  W = hepa_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('hepatitis'))]
)


tau.hat.eval.hepa <- predict(mod_train_cf_hepa, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('hepatitis'))])$predictions

eval.forest.hepa <- causal_forest(
  Y = Y_diab_test, 
  W = hepa_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('hepatitis'))] )

#RATE estimate
rate.cate.hepa <- rank_average_treatment_effect(eval.forest.hepa, tau.hat.eval.hepa)

plot(rate.cate.hepa, main = "TOC: By decreasing estimated CATE (hepatitis)")




# 13) bilirubin

mod_bili_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$bilirubin,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="bilirubin")] 
)
#ATE estimate
average_treatment_effect(mod_bili_cf) 
bili_bin = 1*( DTM_byhand$bilirubin>0)

mod_train_cf_bili = causal_forest(
  Y = Y_diab_train,
  W = bili_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('bilirubin'))]
)
tau.hat.eval.bili <- predict(mod_train_cf_bili, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('bilirubin'))])$predictions
eval.forest.bili <- causal_forest(
  Y = Y_diab_test, 
  W = bili_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('bilirubin'))] )

#RATE estimate
rate.cate.bili <- rank_average_treatment_effect(eval.forest.bili, tau.hat.eval.bili)

plot(rate.cate.bili, main = "TOC: By decreasing estimated CATE (bilirubin)")


# 14) job

mod_job_cf = causal_forest(
  Y = ordered_diab,
  W = DTM_byhand$job,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="job")] 
)
#ATE estimate
average_treatment_effect(mod_job_cf) 

job_bin = 1*( DTM_byhand$job>0)

mod_train_cf_job = causal_forest(
  Y = Y_diab_train,
  W = job_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train)%in% c('job'))]
)
tau.hat.eval.job <- predict(mod_train_cf_job, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('job'))])$predictions
eval.forest.job <- causal_forest(
  Y = Y_diab_test, 
  W = job_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c('job'))] )

#RATE estimate
rate.cate.job <- rank_average_treatment_effect(eval.forest.job, tau.hat.eval.job)

plot(rate.cate.job, main = "TOC: By decreasing estimated CATE (job)")




######################### Models for hypothiroidsm #########################


ordered_hyp <- arrange(disease_data_reduced,Id)$hypot

#1) female

mod_female_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$female,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="female")] 
)
#ATE estimate
average_treatment_effect(mod_female_hyp) 

Y_hyp_train <- (arrange(disease_data_reduced, Id)$hypot)[s]
Y_hyp_test <- (arrange(disease_data_reduced, Id)$hypot)[-s]


female_bin = 1*(DTM_byhand$female >0)



mod_train_cf_female= causal_forest(
  Y = Y_hyp_train,
  W = female_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("female"))]
)

tau.hat.eval.female <- predict(mod_train_cf_female, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("female"))])$predictions
eval.forest.female <- causal_forest(
  Y = Y_hyp_test, 
  W = female_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("female"))] )

#Rate estimate
rate.cate.female <- rank_average_treatment_effect(eval.forest.female, tau.hat.eval.female)

plot(rate.cate.female, 
     main = expression(bold("TOC: By decreasing estimated CATE ") * bolditalic("(female)")))


#2) male

mod_male_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$male,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="male")] 
)
#ATE estimate
average_treatment_effect(mod_male_hyp) 

male_bin = 1*(DTM_byhand$male >0)

mod_train_cf_male= causal_forest(
  Y = Y_hyp_train,
  W = male_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("male"))]
)

tau.hat.eval.male <- predict(mod_train_cf_male, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("male"))])$predictions
eval.forest.male <- causal_forest(
  Y = Y_hyp_test, 
  W = male_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("male"))] )
#RATE estimate
rate.cate.male <- rank_average_treatment_effect(eval.forest.male, tau.hat.eval.male)

plot(rate.cate.male, main = "TOC: By decreasing estimated CATE (male)")



#3) afib

mod_afib_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$afib,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="afib")] 
)
#ATE estimate
average_treatment_effect(mod_afib_hyp) 
afib_bin = 1*(DTM_byhand$afib >0)

mod_train_cf_afib= causal_forest(
  Y = Y_hyp_train,
  W = afib_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("afib"))]
)


tau.hat.eval.afib <- predict(mod_train_cf_afib, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("afib"))])$predictions
eval.forest.afib <- causal_forest(
  Y = Y_hyp_test, 
  W = afib_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("afib"))] )

#RATE estimate
rate.cate.afib <- rank_average_treatment_effect(eval.forest.afib, tau.hat.eval.afib)

plot(rate.cate.afib, main = "TOC: By decreasing estimated CATE (afib)")



#4) tsh
mod_tsh_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$tsh,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="tsh")] 
)
#ATE estimate
average_treatment_effect(mod_tsh_hyp) 

tsh_bin = 1*(DTM_byhand$tsh >0)
mod_train_cf_tsh= causal_forest(
  Y = Y_hyp_train,
  W = tsh_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("tsh"))]
)
tau.hat.eval.tsh <- predict(mod_train_cf_tsh, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("tsh"))])$predictions
eval.forest.tsh <- causal_forest(
  Y = Y_hyp_test, 
  W = tsh_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("tsh"))] )

#RATE estimate
rate.cate.tsh <- rank_average_treatment_effect(eval.forest.tsh, tau.hat.eval.tsh)

plot(rate.cate.tsh, 
     main = expression(bold("TOC: By decreasing estimated CATE ") * bolditalic("(tsh)")))


#5) head

mod_head_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$head,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="head")] 
)
#ATE estimate
average_treatment_effect(mod_head_hyp) 

head_bin = 1*(DTM_byhand$head >0)

mod_train_cf_head= causal_forest(
  Y = Y_hyp_train,
  W = head_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("head"))]
)

tau.hat.eval.head <- predict(mod_train_cf_head, DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("head"))])$predictions

eval.forest.head <- causal_forest(
  Y = Y_hyp_test, 
  W = head_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("head"))] )

#RATE estimate
rate.cate.head <- rank_average_treatment_effect(eval.forest.head, tau.hat.eval.head)
plot(rate.cate.head, main = "TOC: By decreasing estimated CATE (head)")


#6) chf

mod_chf_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$chf,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="chf")] 
)
#ATE estimate
average_treatment_effect(mod_chf_hyp) 

chf_bin = 1*(DTM_byhand$chf >0)

mod_train_cf_chf= causal_forest(
  Y = Y_hyp_train,
  W = chf_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("chf"))]
)

tau.hat.eval.chf <- predict(mod_train_cf_chf,
                            DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("chf"))])$predictions

eval.forest.chf<- causal_forest(
  Y = Y_hyp_test, 
  W = chf_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("chf"))] )

#RATE estimate
rate.cate.chf <- rank_average_treatment_effect(eval.forest.chf, tau.hat.eval.chf)
plot(rate.cate.chf, main = "TOC: By decreasing estimated CATE (chf)")


#7) uti

mod_uti_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$uti,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="uti")] 
)
#ATE estimate
average_treatment_effect(mod_uti_hyp) 

uti_bin = 1*(DTM_byhand$uti >0)

mod_train_cf_uti= causal_forest(
  Y = Y_hyp_train,
  W = uti_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("uti"))]
)

tau.hat.eval.uti <- predict(mod_train_cf_uti,
                            DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("uti"))])$predictions

eval.forest.uti <- causal_forest(
  Y = Y_hyp_test, 
  W = uti_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("uti"))] )

#RATE estimate
rate.cate.uti <- rank_average_treatment_effect(eval.forest.uti, tau.hat.eval.uti)
plot(rate.cate.uti, main = "TOC: By decreasing estimated CATE (uti)")


#8) woman

mod_woman_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$woman,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="woman")] 
)
#ATE estimate
average_treatment_effect(mod_woman_hyp) 

woman_bin = 1*(DTM_byhand$woman >0)

mod_train_cf_woman= causal_forest(
  Y = Y_hyp_train,
  W = woman_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("woman"))]
)

tau.hat.eval.woman <- predict(mod_train_cf_woman,
                              DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("woman"))])$predictions

eval.forest.woman <- causal_forest(
  Y = Y_hyp_test, 
  W = woman_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("woman"))] )

#RATE estimate
rate.cate.woman <- rank_average_treatment_effect(eval.forest.woman, tau.hat.eval.woman)
plot(rate.cate.woman, main = "TOC: By decreasing estimated CATE (woman)")


#9) cardiovascular 

mod_cardiovascular_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$cardiovascular ,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="cardiovascular")] 
)
#ATE estimate
average_treatment_effect(mod_cardiovascular_hyp) 

cardiovascular_bin = 1*(DTM_byhand$cardiovascular>0)

mod_train_cf_cardiovascular = causal_forest(
  Y = Y_hyp_train,
  W = cardiovascular_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("cardiovascular"))]
)

tau.hat.eval.cardiovascular  <- predict(mod_train_cf_cardiovascular ,
                                        DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in%
                                                                   c("cardiovascular"))])$predictions

eval.forest.cardiovascular <- causal_forest(
  Y = Y_hyp_test, 
  W = cardiovascular_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("cardiovascular"))] )

#RATE estimate
rate.cate.cardiovascular  <- rank_average_treatment_effect(eval.forest.cardiovascular , tau.hat.eval.cardiovascular )
plot(rate.cate.cardiovascular , main = "TOC: By decreasing estimated CATE (cardiovascular )")


#10) man 

mod_man_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$man ,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="man")] 
)
#ATE estimate
average_treatment_effect(mod_man_hyp) 

man_bin = 1*(DTM_byhand$man  >0)

mod_train_cf_man = causal_forest(
  Y = Y_hyp_train,
  W = man_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("man"))]
)

tau.hat.eval.man  <- predict(mod_train_cf_man ,
                             DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in%
                                                        c("man"))])$predictions

eval.forest.man <- causal_forest(
  Y = Y_hyp_test, 
  W = man_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("man"))] )

#RATE estimate
rate.cate.man  <- rank_average_treatment_effect(eval.forest.man , tau.hat.eval.man )
plot(rate.cate.man , main = "TOC: By decreasing estimated CATE (man)")


#12) hepatitis 

mod_hepatitis_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$hepatitis,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="hepatitis")] 
)
#ATE estimate
average_treatment_effect(mod_hepatitis_hyp) 

hepatitis_bin = 1*(DTM_byhand$hepatitis>0)

mod_train_cf_hepatitis = causal_forest(
  Y = Y_hyp_train,
  W = hepatitis_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("hepatitis"))]
)

tau.hat.eval.hepatitis  <- predict(mod_train_cf_hepatitis ,
                                   DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("hepatitis"))])$predictions

eval.forest.hepatitis <- causal_forest(
  Y = Y_hyp_test, 
  W = hepatitis_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("hepatitis"))] )

#RATE estimate
rate.cate.hepatitis  <- rank_average_treatment_effect(eval.forest.hepatitis , tau.hat.eval.hepatitis )
plot(rate.cate.hepatitis , main = "TOC: By decreasing estimated CATE (hepatitis )")


#13) baby 

mod_baby_hyp = causal_forest(
  Y = ordered_hyp,
  W = DTM_byhand$baby ,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="baby")] 
)
#ATE estimate
average_treatment_effect(mod_baby_hyp) 

baby_bin = 1*(DTM_byhand$baby>0)

mod_train_cf_baby = causal_forest(
  Y = Y_hyp_train,
  W = baby_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("baby"))]
)

tau.hat.eval.baby  <- predict(mod_train_cf_baby ,
                              DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("baby"))])$predictions

eval.forest.baby <- causal_forest(
  Y = Y_hyp_test, 
  W = baby_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("baby"))] )

#RATE estimate
rate.cate.baby<- rank_average_treatment_effect(eval.forest.baby , tau.hat.eval.baby )
plot(rate.cate.baby, main = "TOC: By decreasing estimated CATE (baby )")


######################### Models for surrenal gland disorders #########################


#------------------------------

#### surrenal  #######

#1) postoperative

ordered_surr <- arrange(disease_data_reduced,Id)$surren

mod_postoperative_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$postoperative ,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="postoperative")] 
)
#ATE estimate
average_treatment_effect(mod_postoperative_surr) 


Y_surr_train <- (arrange(disease_data_reduced, Id)$surren)[s]
Y_surr_test <- (arrange(disease_data_reduced, Id)$surren)[-s]


postoperative_bin = 1*(DTM_byhand$postoperative >0)

mod_train_cf_postoperative= causal_forest(
  Y = Y_surr_train,
  W = postoperative_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("postoperative"))]
)



tau.hat.eval.postoperative <- predict(mod_train_cf_postoperative, 
                                      DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                                 c("postoperative"))])$predictions

eval.forest.postoperative <- causal_forest(
  Y = Y_surr_test, 
  W = postoperative_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("postoperative"))] )

#RATE estimate
rate.cate.postoperative <- rank_average_treatment_effect(eval.forest.postoperative, tau.hat.eval.postoperative)
plot(rate.cate.postoperative, main = "TOC: By decreasing estimated CATE (postoperative)")



#2) insufficiency

mod_insufficiency_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$insufficiency ,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="insufficiency")] 
)
#ATE estimate
average_treatment_effect(mod_insufficiency_surr) 


insufficiency_bin = 1*(DTM_byhand$insufficiency >0)

mod_train_cf_insufficiency= causal_forest(
  Y = Y_surr_train,
  W = insufficiency_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("insufficiency"))]
)



tau.hat.eval.insufficiency <- predict(mod_train_cf_insufficiency, 
                                      DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                                 c("insufficiency"))])$predictions

eval.forest.insufficiency <- causal_forest(
  Y = Y_surr_test, 
  W = insufficiency_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("insufficiency"))] )

#RATE estimate
rate.cate.insufficiency <- rank_average_treatment_effect(eval.forest.insufficiency, tau.hat.eval.insufficiency)
plot(rate.cate.insufficiency, main = "TOC: By decreasing estimated CATE (insufficiency)")


#3) steroid(s)

mod_steroid_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$steroid + DTM_byhand$steroids,
  X = DTM_byhand[, -which(colnames(DTM_byhand_test)%in% c("steroid","steroids"))] 
)
#ATE estimate
average_treatment_effect(mod_steroid_surr) 


steroid_bin = 1*(DTM_byhand$steroid >0)

mod_train_cf_steroid= causal_forest(
  Y = Y_surr_train,
  W = steroid_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("steroid","steroids"))]
)



tau.hat.eval.steroid <- predict(mod_train_cf_steroid, 
                                DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                           c("steroid","steroids"))])$predictions

eval.forest.steroid <- causal_forest(
  Y = Y_surr_test, 
  W = steroid_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("steroid","steroids"))] )

#RATE estimate
rate.cate.steroid <- rank_average_treatment_effect(eval.forest.steroid, tau.hat.eval.steroid)
plot(rate.cate.steroid, main = "TOC: By decreasing estimated CATE (steroid)")




#4) prednisone

mod_prednisone_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$prednisone,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="prednisone")] 
)
#ATE estimate
average_treatment_effect(mod_prednisone_surr) 


prednisone_bin = 1*(DTM_byhand$prednisone>0)

mod_train_cf_prednisone= causal_forest(
  Y = Y_surr_train,
  W = prednisone_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("prednisone"))]
)



tau.hat.eval.prednisone <- predict(mod_train_cf_prednisone, 
                                   DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                              c("prednisone"))])$predictions

eval.forest.prednisone <- causal_forest(
  Y = Y_surr_test, 
  W = prednisone_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("prednisone"))] )

#RATE estimate
rate.cate.prednisone <- rank_average_treatment_effect(eval.forest.prednisone, tau.hat.eval.prednisone)
plot(rate.cate.prednisone, main = "TOC: By decreasing estimated CATE (prednisone)")



#5) hypotension

mod_hypotension_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$hypotension,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="hypotension")] 
)
#ATE estimate
average_treatment_effect(mod_hypotension_surr) 


hypotension_bin = 1*(DTM_byhand$hypotension>0)

mod_train_cf_hypotension= causal_forest(
  Y = Y_surr_train,
  W = hypotension_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("hypotension"))]
)



tau.hat.eval.hypotension <- predict(mod_train_cf_hypotension, 
                                    DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                               c("hypotension"))])$predictions

eval.forest.hypotension <- causal_forest(
  Y = Y_surr_test, 
  W = hypotension_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("hypotension"))] )

#RATE estimate
rate.cate.hypotension <- rank_average_treatment_effect(eval.forest.hypotension, tau.hat.eval.hypotension)
plot(rate.cate.hypotension, main = "TOC: By decreasing estimated CATE (hypotension)")



#6) stress

mod_stress_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$stress,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="stress")] 
)
#ATE estimate
average_treatment_effect(mod_stress_surr) 


stress_bin = 1*(DTM_byhand$stress>0)

mod_train_cf_stress= causal_forest(
  Y = Y_surr_train,
  W = stress_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("stress"))]
)



tau.hat.eval.stress <- predict(mod_train_cf_stress, 
                               DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                          c("stress"))])$predictions

eval.forest.stress <- causal_forest(
  Y = Y_surr_test, 
  W = stress_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("stress"))] )

#RATE estimate
rate.cate.stress <- rank_average_treatment_effect(eval.forest.stress, tau.hat.eval.stress)
plot(rate.cate.stress, main = "TOC: By decreasing estimated CATE (stress)")



#7) vancomycin

mod_vancomycin_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$vancomycin,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="vancomycin")] 
)
#ATE estimate
average_treatment_effect(mod_vancomycin_surr) 


vancomycin_bin = 1*(DTM_byhand$vancomycin>0)

mod_train_cf_vancomycin= causal_forest(
  Y = Y_surr_train,
  W = vancomycin_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("vancomycin"))]
)



tau.hat.eval.vancomycin <- predict(mod_train_cf_vancomycin, 
                                   DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                              c("vancomycin"))])$predictions

eval.forest.vancomycin <- causal_forest(
  Y = Y_surr_test, 
  W = vancomycin_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("vancomycin"))] )

#RATE estimate
rate.cate.vancomycin <- rank_average_treatment_effect(eval.forest.vancomycin, tau.hat.eval.vancomycin)
plot(rate.cate.vancomycin, main = "TOC: By decreasing estimated CATE (vancomycin)")


#8) micu

mod_micu_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$micu,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="micu")] 
)
#ATE estimate
average_treatment_effect(mod_micu_surr) 


micu_bin = 1*(DTM_byhand$micu>0)

mod_train_cf_micu= causal_forest(
  Y = Y_surr_train,
  W = micu_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("micu"))]
)



tau.hat.eval.micu <- predict(mod_train_cf_micu, 
                             DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                        c("micu"))])$predictions

eval.forest.micu <- causal_forest(
  Y = Y_surr_test, 
  W = micu_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("micu"))] )

#RATE estimate
rate.cate.micu <- rank_average_treatment_effect(eval.forest.micu, tau.hat.eval.micu)
plot(rate.cate.micu, main = "TOC: By decreasing estimated CATE (micu)")


#9) coronary

mod_coronary_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$coronary,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="coronary")] 
)
#ATE estimate
average_treatment_effect(mod_coronary_surr) 


coronary_bin = 1*(DTM_byhand$coronary>0)

mod_train_cf_coronary= causal_forest(
  Y = Y_surr_train,
  W = coronary_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("coronary"))]
)



tau.hat.eval.coronary <- predict(mod_train_cf_coronary, 
                                 DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                            c("coronary"))])$predictions

eval.forest.coronary <- causal_forest(
  Y = Y_surr_test, 
  W = coronary_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("coronary"))] )

#RATE estimate
rate.cate.coronary <- rank_average_treatment_effect(eval.forest.coronary, tau.hat.eval.coronary)
plot(rate.cate.coronary, main = "TOC: By decreasing estimated CATE (coronary)")


#10) artery

mod_artery_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$artery,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="artery")] 
)
#ATE estimate
average_treatment_effect(mod_artery_surr) 


artery_bin = 1*(DTM_byhand$artery>0)

mod_train_cf_artery= causal_forest(
  Y = Y_surr_train,
  W = artery_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("artery"))]
)



tau.hat.eval.artery <- predict(mod_train_cf_artery, 
                               DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                          c("artery"))])$predictions

eval.forest.artery <- causal_forest(
  Y = Y_surr_test, 
  W = artery_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("artery"))] )

#RATE estimate
rate.cate.artery <- rank_average_treatment_effect(eval.forest.artery, tau.hat.eval.artery)
plot(rate.cate.artery, main = "TOC: By decreasing estimated CATE (artery)")


#11) aortic

mod_aortic_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$aortic,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="aortic")] 
)
#ATE estimate
average_treatment_effect(mod_aortic_surr) 


aortic_bin = 1*(DTM_byhand$aortic>0)

mod_train_cf_aortic= causal_forest(
  Y = Y_surr_train,
  W = aortic_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("aortic"))]
)



tau.hat.eval.aortic <- predict(mod_train_cf_aortic, 
                               DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                          c("aortic"))])$predictions

eval.forest.aortic <- causal_forest(
  Y = Y_surr_test, 
  W = aortic_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("aortic"))] )

#RATE estimate
rate.cate.aortic <- rank_average_treatment_effect(eval.forest.aortic, tau.hat.eval.aortic)
plot(rate.cate.aortic, main = "TOC: By decreasing estimated CATE (aortic)")


#12) valve

mod_valve_surr = causal_forest(
  Y = ordered_surr,
  W = DTM_byhand$valve,
  X = DTM_byhand[, -which(colnames(DTM_byhand)=="valve")] 
)
#ATE estimate
average_treatment_effect(mod_valve_surr) 


valve_bin = 1*(DTM_byhand$valve>0)

mod_train_cf_valve= causal_forest(
  Y = Y_surr_train,
  W = valve_bin[s],
  X = DTM_byhand_train[, -which(colnames(DTM_byhand_train) %in% c("valve"))]
)



tau.hat.eval.valve <- predict(mod_train_cf_valve, 
                              DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% 
                                                         c("valve"))])$predictions

eval.forest.valve <- causal_forest(
  Y = Y_surr_test, 
  W = valve_bin[-s],
  X = DTM_byhand_test[, -which(colnames(DTM_byhand_test)%in% c("valve"))] )

#RATE estimate
rate.cate.valve <- rank_average_treatment_effect(eval.forest.valve, tau.hat.eval.valve)
plot(rate.cate.valve, main = "TOC: By decreasing estimated CATE (valve)")


# -- end CAUSAL FOREST MODELS ------------------------------------------

