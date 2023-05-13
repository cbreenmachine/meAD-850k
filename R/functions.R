


create_id_mapping <- function(master.file, array.file){

  master.df <- read_csv(master.file, show_col_types = F)
  array.df <- read_csv(array.file, show_col_types = F)

  inner_join(master.df, array.df, by = c("study_id" = "adrcnum...15")) %>%
    dplyr::transmute(sample_id, study_id, array_id = str_split(Basename, "/")) %>%
    unnest(array_id) %>%
    dplyr::filter(str_detect(array_id, "_"))
}

load_targets_from_csv <- function(file){
  read_csv(file, show_col_types=F)
}


load_from_color_channel_data <- function(idir, targets.df){
  my.files <- list.files(path = idir, full.names = T, recursive = T, pattern = "Grn.idat")
  my.files.stripped <- stringr::str_remove(basename(my.files), "_Grn.idat")

  valid.IDs <- intersect(my.files.stripped, targets.df$ID)
  print(valid.IDs)

  # Filter file paths and samplesheet (targets)
  my.files <- my.files[my.files.stripped %in% valid.IDs]
  targets.df <- targets.df %>% dplyr::filter(ID %in% valid.IDs)

  # Heavy lifting--read RGset here...
  RGset <- minfi::read.metharray(my.files, verbose = TRUE)
  return(RGset)
}


get_detection_pvals <- function(RGset){
  minfi::detectionP(RGset)
}


filter_out_bad_probes <- function(m.set, detection.pvals){
  # Check that probes are in correct order
  detection.pvals <- detection.pvals[match(featureNames(m.set), rownames(detection.pvals)), ]

  if (!all(rownames(detection.pvals) == featureNames(m.set))){
    error("Probes are not ordered correctly")
  }

  # First, based on detection p-values
  keep.ix <- rowSums(detection.pvals < 0.01) == ncol(m.set)
  m.set.2 <- m.set[keep.ix,]

  # Map and drop SNPs
  m.set.3 <- mapToGenome(m.set.2) %>%
    dropLociWithSnps(maf = 0.05) %>%
    dropMethylationLoci(dropRS = T, dropCH = T)

  return(m.set.3)
}

get_methylation_estimates <- function(m.set.filtered, ids.df){

  loci <- getLocations(m.set.filtered) %>%
    as.data.frame() %>%
    dplyr::select(seqnames, start) %>%
    dplyr::mutate(end = start + 2)

  methylation <- getBeta(m.set.filtered) %>%
    as.data.frame()

  out <- cbind(loci, methylation)

  # Now we renamse the columns to RA id instead of beadchip (array)
  ix <- match(colnames(out)[-3:-1], ids.df$array_id)

  #Check
  all(ids.df$array_id[ix] == colnames(out)[-3:-1])

  colnames(out)[-3:-1] <- ids.df$sample_id[ix]
  out
}


download_chain_from_ucsc <- function(dir, file){
  # Assign output name
  ofile.name <- file.path(dir, basename(file))

  # Download
  download.file(file, destfile = ofile.name)

  # Decompress
  R.utils::gunzip(ofile.name, remove = F, overwrite = T)

  # Return a chain object
  rtracklayer::import.chain(stringr::str_remove(ofile.name, "\\.gz"))
}



liftover_wrapper <- function(gr, chain){
  seqlevelsStyle(gr) <- "UCSC"

  unlist(rtracklayer::liftOver(gr, chain))
}


lift_m_estimates_to_hg38 <- function(m.estimates.df, chain){
  m.gr <- makeGRangesFromDataFrame(m.estimates.df,
                                   starts.in.df.are.0based = T,
                                   keep.extra.columns = T)

  liftover_wrapper(m.gr, chain)
}



write_bed <- function(data.gr, file){

  names(data.gr) <- NULL

  out <- data.gr %>%
    as.data.frame() %>%
    dplyr::select(-width) %>%
    dplyr::rename(chr = seqnames)

  colnames(out)[5:ncol(out)] <- colnames(mcols(data.gr))

  fwrite(out, file, sep = "\t")
  return(file)
}

# filter_out_bad_probes_and_samples <- function(RGset, targets){
#   det.p <- minfi::detectionP(RGset)
#   mean.by.sample <- colMeans(det.p)
#
#   # Samples
#   # Keep if the mean p-value is below 0.05
#   samples.ix <- which(mean.by.sample < 0.05)
#   probes.ix <-
#
#   # Filter
#   good.samples <- names(which(mean.p <= 0.05))
#   RGset <- RGset[ ,good.samples]
#   targets <- dplyr::filter(targets, ID %in% good.samples)
#   return(list("RGset" = RGset, "targets" = targets, "det.p" = det.p))
# }
#
# estimate_cell_comp <- function(RGset){
#   cellCounts <- FlowSorted.Blood.EPIC::estimateCellCounts2(RGset, compositeCellType="Blood", referencePlatform = "IlluminaHumanMethylationEPIC")
#   cc.df <- as.data.frame(cellCounts)  %>%
#     rownames_to_column("ID") %>%
#     rename_with( function(x) tolower(str_remove(x, "counts.")), starts_with("counts"))
#   return(cc.df)
# }
#
#
# preprocess <- function(RGset, targets, method){
#
#   if (!all(colnames(RGset) == targets$ID)){
#     error("Column names in RGset and targets$ID don't match--remake RGset with wrapper")
#   } else {
#     print("Names match!")
#   }
#
#   # First drop samples with uniformly bad detection
#   out <- drop_bad_detection_samples(RGset, targets)
#   RGset <- out$RGset
#   targets <- out$targets
#   det.p <- out$det.p
#
#   # SWAN or NOOB preprocessing
#   mSetProcessed <- normalize(out$RGset, method)
#   mSetFiltered <- map_and_filter(mSetProcessed)
#   targets <- out$targets
#   pSex.df <- out$pSex.df
#
#   # Sex prediction
#   out <- predict_sex(mSetFiltered)
#   mSetFiltered <- out$mSetFiltered
#   targets <- out$targets
#   pSex.df <- out$pSex.df
#
#   # One more round of detection p-value culling
#   mSetFiltered <- drop_bad_detection_probes(SetFiltered)
#
#   out <- get_methylation_estimates(mSetFiltered)
#   M <- out$M
#   beta <- out$M
#
#   # Parallel filtering
#   RGset <- RGset[, colnames(mSetFiltered)]
#   cellcomp.df <- estimate_cell_comp(RGset)
#
#   return(list("M" = M, "beta" = beta, "targets" = targets, "cell.comp" = cellcomp.df, "sex.pred" = pSex.df))
# }
