library(biogram)

tar <- sample(0L:1, 100, replace = TRUE)
feat <- sample(0L:1, 100, replace = TRUE)
calc_ig(feat, tar, 100, sum(tar))

tar_names <- c("apc", "amp", "neg")
tar <- unlist(lapply(tar_names, rep, 50))
feat <- unlist(lapply(c(0.8, 0.1, 0.1), function(ith_prob) {
  sample(0L:1, 50, replace = TRUE, prob = c(ith_prob, 1 - ith_prob)) 
}))

do.call(rbind, lapply(combn(tar_names, 2, simplify = FALSE), function(ith_cmbn) {
  data.frame(tar1 = ith_cmbn[1], tar2 = ith_cmbn[2],
             ig = calc_ig(feature = feat[tar %in% ith_cmbn], 
                          target = tar[tar %in% ith_cmbn] == ith_cmbn[1],
                          len_target = length(tar[tar %in% ith_cmbn] == ith_cmbn[1]),
                          pos_target = sum(tar[tar %in% ith_cmbn] == ith_cmbn[1])))
}))
