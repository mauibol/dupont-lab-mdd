# Get file paths for raw data
# This might need to be updated as we get more data

library(glue)

get_raw_path_visium <- function(i) {
  base_path <- glue("{RAW_DATA_PATH}/250131_MAURA_MADELINE_7_HUMAN_LIBRARY_SELF_3B_PE150_X/")
  sub_dir <- glue("{i}/analysis/250131_MAURA_MADELINE_7_HUMAN_3B_PE150_X-{i}-spaceranger-count/")
  
  if (i %in% c('RM005', 'RM006', 'RM007', 'RM008', 'RM009')){
    base_path <- glue('{RAW_DATA_PATH}/250417_MAURA_RAKSHITHA_2_HUMAN_LIBRARY_SELF_3B_PE150_X/analysis/')
    sub_dir <- glue("250417_MAURA_RAKSHITHA_2_HUMAN_LIBRARY_SELF_3B_PE150_X-{i}-spa/")
  }
  
  if (i == 'RM010') {
    base_path <- glue('{RAW_DATA_PATH}/250417_MAURA_RAKSHITHA_2_HUMAN_LIBRARY_SELF_3B_PE150_X/analysis_REANALYSIS/')
    sub_dir <- glue("250417_MAURA_RAKSHITHA_2_HUMAN_LIBRARY_SELF_3B_PE150_X-{i}-spa/")
  }  
  
  else if (substr(i, 1,4)=='MM07') {
    base_path <- glue('{RAW_DATA_PATH}/230915_MAURA_MADELINE_SPATIAL/')
    sub_dir <- glue("{i}/outs/")
  }
  path <- paste0(base_path, sub_dir)
  return(path)
}

get_raw_path_multiome <- function(i) {
  if (substr(i, 1,2)=='CM') {
    dir <- "/ngs/release/singleCell/220502_MAURA_CHEICK_16_HUMAN_SEQUENCING/"
    sub_dir <- glue("{i}/analysis/220502_MAURA_CHEICK_16_HUMAN_SEQUENCING-{i}-cellranger-arc-count-default/{i}_cellranger_arc_count_outs")
  }
  
  else if (substr(i, 1,3)=='MM0' & grepl("^(0[1-9]|[12][0-9]|3[01])$", substr(i, 4,5))) {
    dir <- "/230817_MAURA_MADELINE_12_HUMAN_10X/"
    sub_dir <- glue('{i}/analysis/230817_MAURA_MADELINE_12_HUMAN_10X-{i}-cellranger-arc-count-default/{i}_cellranger_arc_count_outs')
  }
  
  else if (substr(i, 1,3)=='MM1' & grepl("^(4[5-9]|5[0-2])$", substr(i, 4,5))) {
    dir <- glue("/250131_MAURA_MADELINE_7_HUMAN_LIBRARY_SELF_3B_PE150_X/{i}/analysis/250131_MAURA_MADELINE_7_HUMAN_LIBRARY_SELF_3B_PE150_X-{i}-cellranger-arc-count-default/")
    sub_dir <- glue("{i}_cellranger_arc_count_outs")
  }
  else {
    dir <- "/240806_MAURA_MADELINE_6_HUMAN_LIBRARY_10X/"
    sub_dir <- glue('{i}/analysis/240806_MAURA_MADELINE_6_HUMAN_LIBRARY_10X-{i}-cellranger-arc-count-default/{i}_cellranger_arc_count_outs')
  }
  path <- paste0(RAW_DATA_PATH, dir, sub_dir)
  return(path)
}
