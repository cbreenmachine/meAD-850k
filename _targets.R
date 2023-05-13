
library(targets)

tar_source("R/")


# Where to store reference data (e.g. chain files)
DATA.REFERENCE.DIR <- "DataRef"


# Set target options:
tar_option_set(
  packages = c(
    "data.table",
    "tidyverse",
    "minfi",
    "IlluminaHumanMethylationEPICmanifest",
    "data.table",
    "FlowSorted.Blood.EPIC",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
    ),
  format = "rds" # default storage format
)

# Replace the target list below with your own:
list(
  tar_target(chain.19to38,
              download_chain_from_ucsc(
                DATA.REFERENCE.DIR,
                "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz")),
  tar_target(ids.df,
             create_id_mapping(
               "DataRaw/masterSamplesheet.csv",
               "DataRaw/EPIC-850k-run-info.csv"
               )),
  tar_target(targets.df, load_targets_from_csv("DataRaw/targets.csv")),
  tar_target(rg.set, load_from_color_channel_data("DataRaw/EPICArrays/", targets.df)),
  tar_target(detection.pvals, get_detection_pvals(rg.set)),
  # If not using 84, need to filter samples here
  tar_target(m.set, preprocessSWAN(rg.set)),
  tar_target(m.set.filtered, filter_out_bad_probes(m.set, detection.pvals)),
  tar_target(m.estimates.df, get_methylation_estimates(m.set.filtered, ids.df)),
  tar_target(m.hg38.gr, lift_m_estimates_to_hg38(m.estimates.df, chain.19to38)),
  tar_target(m.hg38.bed, write_bed(m.hg38.gr, "DataDerived/array.M.hg38.bed"), format = "file")
)
