clean_dx <- function(vec) {
  dplyr::case_when(
    grepl("Control", vec, ignore.case = TRUE) ~ "Control",
    grepl("Kawasaki|KD", vec, ignore.case = TRUE) ~ "KD",
    grepl("Uncertain", vec, ignore.case = TRUE) ~ "Uncertain",
    TRUE ~ "Other"
  )
}