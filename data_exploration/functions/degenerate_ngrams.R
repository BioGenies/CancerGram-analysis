
#' Create all combinations of traits
#'
#' Creates all combinations of input traits.
#' @param vtraits a vector of trait indices in the expanded aaindex table
#' create by \code{\link{choose_properties}}.
#'
#' @return a list of traits combinations of given length (from single trait to
#' all traits).


create_traits_combination <- function(ftraits) {
  #vector of traits
  vtraits <- unlist(ftraits)
  
  #all combinations of traits
  #number of combinations: sum(sapply(all_traits_combn_list, nrow))
  pblapply(1L:length(vtraits), function(i)
    t(combn(vtraits, i)))
}



#' Create alphabets
#'
#' Creates alphabets (3-6 groups long) from list of traits.
#' @param ftraits a vector of trait indices in the expanded aaindex table
#' create by \code{\link{choose_properties}}.
#' @param list_duplicates if \code{TRUE} returns also a list of duplicates.
#' @return a named vector of alphabets (for example 
#' \code{iknty_degpqrs_acfhlmvw})
create_alphabets <- function(ftraits, list_duplicates = FALSE) {
  
  paste_enc <- function(x)
    paste0(sapply(x, paste0, collapse = ""), collapse = "_")
  
  grouping_properties <- aaprop[ftraits, ]
  all_traits_combn_list <- create_traits_combination(ftraits)
  
  #create alphabets
  all_aa_groups <- {
    res <- unlist(lapply(all_traits_combn_list, function(all_traits_combn)
      vapply(1L:nrow(all_traits_combn), function(single_trait_combn) {
        cl <- t(aaprop[unlist(all_traits_combn[single_trait_combn, , drop = FALSE]), , drop = FALSE]) %>%
          dist %>%
          hclust(method = "ward.D2")
        #cl <- hclust(dist(t(aaprop[unlist(all_traits_combn[single_trait_combn, , drop = FALSE]), , drop = FALSE])))
        gr <- cutree(cl, k = 6)
        names(gr) <- tolower(names(gr))
        agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
        #inside alphabets, amino acids are ordered alphabetically
        agg_gr <- lapply(agg_gr, sort)
        #groups are sorted by their length
        paste_enc(agg_gr[order(lengths(agg_gr))])
      }, "a")))
    names(res) <- paste0("ID", 1L:length(res), "K", 6)
    res
  }
  
  #get indices of unique alphabets
  aa_id <- lapply(all_aa_groups, function(i) !duplicated(i))
  
  aa_duplicates <- unlist(lapply(1L:length(aa_id), function(i) 
    lapply(all_aa_groups[[i]][aa_id[[i]]], function(j)
      names(which(j == all_aa_groups[[i]])))
  ), recursive = FALSE)
  #aa_duplicates <- aa_duplicates[lengths(aa_duplicates) > 1]
  
  #remove from aa_groups redundant alphabets
  aa_groups <- unlist(lapply(1L:length(aa_id), function(i) {
    all_aa_groups[[i]][aa_id[[i]]]
  }), recursive = FALSE)
  
  #add as a benchmark two alphabets from the literature
  aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"), 
             `2` = c("k", "r", "h"), 
             `3` = c("d", "e"), 
             `4` = c("f", "w", "y", "s", "t", "c", "n", "q"))
  
  aa2 = list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
             `2` = c("k", "r", "h"), 
             `3` = c("d", "e"), 
             `4` = c("s", "t", "c", "n", "q", "y", "w"))
  
  if(list_duplicates) {
    list(aagroups = c(aa1 = paste_enc(aa1), aa2 = paste_enc(aa2), aa_groups),
         aa_duplicates = aa_duplicates,
         features = ftraits)
  } else {
    list(aagroups = c(aa1 = paste_enc(aa1), aa2 = paste_enc(aa2), aa_groups),
         features = ftraits)
  }
}

#create_alphabets()




