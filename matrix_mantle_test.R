library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(vegan)
library(iNEXT)
library(bipartiteD3)
library(bipartite)
library(betapart)
library(usedist)
library(readxl)
library(writexl)
library(openxlsx)

# loading data and prepare data -all parasitoids and caterpillars without singletons and doubletons
main_table <- readxl::read_excel("data/input/Beta_August22.xlsx")
main_table <- dplyr::rename(main_table, par_morph = NEW_M_code, cat_morph = cat_whole_name, plant_sp = sample_code)

# Main dataset - reduced
testik <- main_table %>%
  dplyr::filter(Par_remove == "0", Remove_cats == "0", singlecat == "0", singlepara == "0") %>%
  dplyr::select(LOC2, par_morph, cat_morph, plant_sp)

# distance matrix
Distance <-
  read.csv2(
    here::here("data/input/Distance.csv"),
    row.names = 1
  )

#######   Food-web interactions between caterpillars and parasitoids

# Preparing contingency table - caterpillars vs parasitoids

count_indiv <- function(data_source,
                        sel_locality) {
  data_filter <-
    data_source %>%
    dplyr::filter(LOC2 == sel_locality) %>%
    dplyr::select(cat_morph, par_morph)

  result <- as.matrix(table(data_filter$cat_morph, data_filter$par_morph))

  return(result)
}

count_indiv(
  data_source = testik,
  sel_locality = "Elem"
)

loc_list <-
  testik$LOC2 %>%
  unique() %>%
  rlang::set_names() %>%
  purrr::map(
    .x = .,
    .f = ~ count_indiv(
      data_source = testik,
      sel_locality = .x
    )
  )

str(loc_list, 2)
loc_list["Elem"]


################################# run the analyses

FW_PC_res <-
  bipartite::betalinkr_multi(
    bipartite::webs2array(loc_list),
    partitioning = "commondenom",
    partition.st = TRUE
  ) %>%
  tibble::add_column(dataset = "Cat_Par") # adding colums for excel comparison


######### Prepare tables for Mantle test

# function, which prepare the matrix suitable for the Mantle test

make_a_matrix <- function(data_source,
                          var_name) {
  xmnames <- unique(c(data_source$i, data_source$j))
  xmnames_length <- length(xmnames)
  xm <- matrix(0, nrow = xmnames_length, ncol = xmnames_length)

  colnames(xm) <- xmnames
  rownames(xm) <- xmnames

  xm[lower.tri(xm, diag = FALSE)] <- purrr::pluck(FW_PC_res, var_name)

  xm[upper.tri(xm, diag = FALSE)] <- purrr::pluck(FW_PC_res, var_name)

  return(xm)
}

make_a_matrix(
  data_source = FW_PC_res,
  var_name = "S"
)

data_for_mantel <-
  FW_PC_res %>%
  dplyr::select(-c(i, j, dataset)) %>%
  names() %>%
  rlang::set_names() %>%
  purrr::map(
    .f = ~ make_a_matrix(
      data_source = FW_PC_res,
      var_name = .x
    )
  )

str(data_for_mantel, 2)

result_mantel <-
  purrr::map(
    .x = data_for_mantel,
    .f = ~ vegan::mantel(
      xdis = .x,
      ydis = Distance,
      method = "spear"
    )
  )

str(result_mantel, 2)

result_table <-
  tibble::tibble(
    names(result_mantel)
  ) %>%
  dplyr::mutate(
    mantel_stat = purrr::map_dbl(
      .x = result_mantel,
      .f = ~ .x$statistic
    ),
    p_value = purrr::map_dbl(
      .x = result_mantel,
      .f = ~ .x$signif
    ),
  )

readr::write_csv(
  result_table,
  here::here("data/output/mantel_result.csv")
)
