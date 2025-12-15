suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(openxlsx))

# Args for testing
args <- list()
args$input <- "output/pf=pf8-version/all_hspf_analyses_summary.tsv"
args$output <- "output/pf=pf8-version/all_hspf_analyses_summary.xlsx"

# Read TSV
result <- read_tsv(args$input)

# Prepare 'out' with correct column names and empty column
out <- result %>%
  transmute(
    `Reported in` = Reported,
    `Cell type` = celltype,
    `Cell size (degrees)` = cellsize,
    `HbS map: r` = HbSr0,
    `HbS map: sigma` = HbSsigma0,
    allele = allele,
    area = area,
    Countries = countries,
    `Min distance (km) to HbS survey point` = min_km_to_survey_pt,
    `Number of Pf data points` = n_data_points,
    `Covariate in HbS-Pf model` = covariate,
    `Δf+\n(mean)` = round(delta_mean, 3),
    `Δf+\n(median)` = round(delta_median, 3),
    `Δf+\n(2.5% quantile)` = round(delta_q2.5, 3),
    `Δf+\n(97.5% quantile)` = round(delta_q97.5, 3)
  ) %>%
  mutate(blank = NA) %>%  # empty column after Δf+ columns
  bind_cols(
    result %>% transmute(
      `5%` = round(pf_at_0.05, 3),
      `10%` = round(pf_at_0.1, 3),
      `15%` = round(pf_at_0.15, 3),
      `20%` = round(pf_at_0.2, 3),
      `25%` = round(pf_at_0.25, 3),
      `30%` = round(pf_at_0.3, 3)
    )
  )

# Create workbook
wb <- createWorkbook()
addWorksheet(wb, "ST3 results summary")

title <- "ST3 Summary of geographical regression analysis results"
writeData(wb, "ST3 results summary", x = title, startCol = 1, startRow = 1)
mergeCells(wb, "ST3 results summary", cols = 1:8, rows = 1)
title_style <- createStyle(fontSize = 11, textDecoration = "bold", halign = "left")
addStyle(wb, "ST3 results summary", title_style, rows = 1, cols = 1)

# Create header style
header_style <- createStyle(
#  fontSize = 11,
  textDecoration = "bold",
  halign = "center",      # Horizontal alignment
  valign = "bottom",      # Vertical alignment (bottom)
  fgFill = "gray95",     # Light grey background
#  border = "TopBottomLeftRight",  # Add borders for better visibility,
  wrapText = TRUE 
)

# Add merged row with "Slope estimates"
writeData(wb, "ST3 results summary", 
          x = "Slope estimates", 
          startRow = 3, 
          startCol = 12)

# Merge cells for the slope estimates header
mergeCells(wb, "ST3 results summary", cols = 12:15, rows = 3)

# Apply the same header style to this merged cell
addStyle(wb, "ST3 results summary", header_style, rows = 3, cols = 12:15, gridExpand = TRUE)

writeData(wb, "ST3 results summary", 
          x = "Modelled Pfsa+ frequencies at HbAS/SS frequency given by header", 
          startRow = 3, 
          startCol = 17)

# Merge cells for the slope estimates header
mergeCells(wb, "ST3 results summary", cols = 17:ncol(out), rows = 3)

# Apply the same header style to this merged cell
addStyle(wb, "ST3 results summary", header_style, rows = 3, cols = 17:ncol(out), gridExpand = TRUE)


# Write column names at row 4 with styling
writeData(wb, "ST3 results summary", 
          x = as.data.frame(t(names(out))),  # Transpose column names to write horizontally
          startRow = 4, 
          startCol = 1,
          colNames = FALSE)

# Apply header style to the column names row
addStyle(wb, "ST3 results summary", header_style, rows = 4, cols = 1:15, gridExpand = TRUE)
addStyle(wb, "ST3 results summary", header_style, rows = 4, cols = 17:ncol(out), gridExpand = TRUE)
# Merge cells vertically to create 3-row vertical space for headers
for (col in 1:15) {
  mergeCells(wb, "ST3 results summary", cols = col, rows = 4:6)
}
for (col in 17:ncol(out)) {
  mergeCells(wb, "ST3 results summary", cols = col, rows = 4:6)
}
# Set column width for column 16 to make it narrower
setColWidths(wb, "ST3 results summary", cols = 16, widths = 3)  # Adjust width as needed

# When writing column names, set column 16 name to empty string
col_names <- names(out)
col_names[16] <- ""  # Remove column name for column 16

# Use the modified column names when writing headers
writeData(wb, "ST3 results summary", 
          x = as.data.frame(t(col_names)),  # Use modified column names
          startRow = 4, 
          startCol = 1,
          colNames = FALSE)

# Write actual data starting at row 7 (after the 3-row header space)
writeData(wb, "ST3 results summary", out, startRow = 7, colNames = FALSE)
# # Create a blank style for column 16
blank_style <- createStyle(
  fgFill = "#FFFFFF",  # White background (or use NULL for no fill)
)

# Apply blank style to column 16
addStyle(wb, "ST3 results summary", blank_style,rows = 1:nrow(out)+6, cols = 16, gridExpand = TRUE)
# Save workbook
saveWorkbook(wb, args$output, overwrite = TRUE)
