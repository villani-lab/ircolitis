cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9",
  "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7", "#000000"
)

my_colors <- list(

  class = c(
    "CPI colitis"    = "#E31A1C",
    "CPI no colitis" = "#A6CEE3",
    "Normal"         = cbPalette[4]
  ),

  class_short = c(
    "C"  = "#E31A1C",
    "NC" = "#A6CEE3",
    "N"  = cbPalette[4]
  ),

  case = c(
    "Case"  = "#E31A1C",
    "Control" = "#1F78B4"
  ),

  chemistry = c(
    "3p"   = cbPalette[2],
    "3pV2" = cbPalette[2],
    "3pV3" = cbPalette[2],
    "5pV2" = cbPalette[3],
    "5p"   = cbPalette[3]
  )
)
