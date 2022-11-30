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
main_table <- readxl::read_excel("data/input/Beta_August22.xlsx") %>%
  dplyr::rename(par_morph = NEW_M_code, cat_morph = cat_whole_name, plant_sp = sample_code)

# Main dataset - reduced
testik <-
  main_table %>%
  dplyr::filter(Par_remove == "0", Remove_cats == "0", singlecat == "0", singlepara == "0") %>%
  dplyr::select(LOC2, par_morph, cat_morph, plant_sp)

# distance matrix
Distance <-
  read.csv2(
    here::here("data/input/Distance.csv"),
    row.names = 1
  )


# helper fnc
count_indiv <- function(data_source,
                        sel_locality,
                        group_a,
                        group_b) {
  data_filter <-
    data_source %>%
    dplyr::filter(LOC2 == sel_locality) %>%
    dplyr::select(
      dplyr::all_of(
        c(
          group_a,
          group_b
        )
      )
    )

  result <- as.matrix(
    table(
      data_filter %>%
        purrr::pluck(group_a),
      data_filter %>%
        purrr::pluck(group_b)
    )
  )

  return(result)
}

make_a_matrix <- function(data_source,
                          var_name) {
  xmnames <- unique(c(data_source$i, data_source$j))
  xmnames_length <- length(xmnames)
  xm <- matrix(0, nrow = xmnames_length, ncol = xmnames_length)

  colnames(xm) <- xmnames
  rownames(xm) <- xmnames

  xm[lower.tri(xm, diag = FALSE)] <- purrr::pluck(data_source, var_name)

  xm[upper.tri(xm, diag = FALSE)] <- purrr::pluck(data_source, var_name)

  return(xm)
}

#######   Food-web interactions between caterpillars and parasitoids

# Preparing contingency table - caterpillars vs parasitoids

mantel_for_group <- function(data_source,
                             data_distance,
                             group_a,
                             group_b) {
  sel_group <-
    paste(
      group_a,
      group_b,
      sep = "_"
    )

  message(
    paste(
      "groups:",
      sel_group
    )
  )

  message(" - count indiv")

  loc_list <-
    data_source$LOC2 %>%
    unique() %>%
    rlang::set_names() %>%
    purrr::map(
      .x = .,
      .f = ~ count_indiv(
        data_source = testik,
        sel_locality = .x,
        group_a = group_a,
        group_b = group_b
      )
    )

  # run the analyses

  message(" - betalinkr_multi")

  FW_PC_res <-
    bipartite::betalinkr_multi(
      bipartite::webs2array(loc_list),
      partitioning = "commondenom",
      partition.st = TRUE
    ) %>%
    # adding colums for excel comparison
    tibble::add_column(dataset = sel_group)

  ######### Prepare tables for Mantle test

  # function, which prepare the matrix suitable for the Mantle test
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

  message(" - mantel test")

  result_mantel <-
    purrr::map(
      .x = data_for_mantel,
      .f = ~ vegan::mantel(
        xdis = .x,
        ydis = Distance,
        method = "spear"
      )
    )

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
      dataset = sel_group
    )

  return(result_table)
}

all_groups_result <-
  tibble::tribble(
    ~group_a, ~group_b,
    "par_morph", "cat_morph",
    "cat_morph", "plant_sp",
    "par_morph", "plant_sp"
  ) %>%
  dplyr::mutate(
    res_table = purrr::map2(
      .x = group_a,
      .y = group_b,
      .f = ~ mantel_for_group(
        data_source = testik,
        data_distance = Distance,
        group_a = .x,
        group_b = .y
      )
    )
  )

readr::write_csv(
  all_groups_result %>%
    tidyr::unnest(res_table),
  here::here("data/output/mantel_result.csv")
)
