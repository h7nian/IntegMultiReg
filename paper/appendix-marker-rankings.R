# Appendix G selection rankings only; posterior-mode coefficients are separate.
appendix_marker_rankings <- function(fit, original, out_dir) {
  stopifnot(identical(fit$model$platform_names, names(original$platforms)))
  counts <- c(mrna = 10L, mirna = 6L, methylation = 10L)
  equations <- c('111' = 'E1', '011' = 'E2', '101' = 'E3', '001' = 'E4')
  all <- list()
  for (p in seq_along(fit$model$platform_names)) {
    platform <- fit$model$platform_names[p]
    feature <- fit$model$feature_names[[p]]
    raw <- original$paper_alignment$original_feature_names[[platform]]
    stopifnot(platform %in% names(counts), length(raw) == length(feature),
      identical(feature, names(original$platforms[[platform]])[-1L]))
    prob <- fit$posterior$inclusion_probabilities[[p]]
    subgroup <- fit$model$subgroup_names[fit$model$platform_subgroups[[p]]]
    stopifnot(nrow(prob) == length(subgroup), ncol(prob) == length(feature),
      all(is.finite(prob)), all(prob >= 0 & prob <= 1), all(subgroup %in% names(equations)))
    for (s in seq_along(subgroup)) {
      idx <- order(-prob[s, ], seq_along(feature))
      n <- min(counts[[platform]], length(feature))
      cutoff <- prob[s, idx[n]]
      all[[length(all) + 1L]] <- data.frame(platform, subgroup = subgroup[s],
        equation = unname(equations[subgroup[s]]), rank = seq_along(idx),
        feature_index = idx, feature = feature[idx], original_name = raw[idx],
        mpip = prob[s, idx], in_top_list = seq_along(idx) <= n,
        tied_at_cutoff = prob[s, idx] == cutoff & sum(prob[s, ] == cutoff) > 1L)
    }
  }
  full <- do.call(rbind, all)
  key <- paste(full$platform, full$feature_index, sep = ':')
  count <- tapply(as.integer(full$in_top_list), key, sum)
  full$n_top_subgroups <- as.integer(count[key])
  top <- full[full$in_top_list, ]
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(full, file.path(out_dir, 'all-marker-rankings.csv'), row.names = FALSE)
  write.csv(top, file.path(out_dir, 'top-marker-rankings.csv'), row.names = FALSE)
  writeLines(c('Rank within each platform and availability subgroup by descending mPIP.',
    'Top counts follow Appendix G: 10 genes, 6 miRNAs, 10 methylation probes.',
    'Ties are ordered by archived feature-column position; cutoff ties are flagged.',
    'Original names are mapped by validated feature-column order, not approximate name matching.',
    'n_top_subgroups counts membership in these top lists; it is not a threshold-based selection count.',
    'No posterior-mode coefficient is supplied by this exporter. coef(imr) is mPIP, not a regression slope.',
    'Historical coefficient conditioning/model choices and IPA remain separate unresolved work.'),
    file.path(out_dir, 'INTERPRETATION.txt'))
  invisible(list(all = full, top = top))
}
