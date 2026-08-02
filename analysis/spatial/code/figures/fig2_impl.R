area_mapping <- function() {
	tibble::tribble(
		~area,                                          ~Region,                                        ~order, ~include, ~parent,
		"global",                                       "Global",                                            1,        1, "Global",
		"africa",                                       "Africa",                                            1,        1, "Global",
		"waf",                                          "West Africa, Cameroon and Gabon",                   2,        1, "Africa",
		"wwaf",                                         "Western populations\n(Western area)",               3,        1, "West Africa",
		"ewaf",                                         "Western populations\n(Eastern area)",               3,        1, "West Africa",
		"gambia+senegal",                               "Gambia\n& Senegal",                                 4,        0, "West Africa",
		"mali",                                         "Mali",                                              4,        0, "West Africa",
		"ghana",                                        "Ghana",                                             4,        0, "West Africa",
		"ghana+burkina+togo",                           "Ghana, Burkina Faso\nand Togo",                     4,        0, "West Africa",
		"ghana+burkina+togo+benin+ivorycoast",          "Ghana, Burkina Faso,\nTogo, Benin and Ivory Coast", 4,        0, "West Africa",
    "nigeria",                                      "Nigeria",                                           4,        0, "West Africa",
		"caf",                                          "Central Africa",                                    2,        0, "Africa",
		"DRC+eaf",                                      "East Africa and DRC",                               2,        1, "Africa",
		"DRC",                                          "DRC",                                               4,        1, "Central Africa",
		"eaf",                                          "East Africa",                                       4,        1, "Africa",
		"tanzania+kenya+uganda+rwanda",                 "Tanzania, Kenya,\nUganda and Rwanda",               4,        0, "East Africa",
		"uganda",                                       "Uganda",                                            4,        0, "East Africa",
		"tanzania",                                     "Tanzania",                                          4,        0, "East Africa",
    "mozambique",                                   "Mozambique",                                        4,        0, "Eest Africa"
	)
}

# Generalised link function
gl <- function(v, parameters) {
  x <- parameters[["intercept"]] + parameters[["beta"]] * v
  nu <- exp(parameters[["log_nu"]])
  1 / (1 + exp(-x))^(1 / nu)
}

calc_slope <- function(intercept, beta, log_nu) {
  gl(0.2, list(intercept = intercept, beta = beta, log_nu = log_nu)) -
    gl(0.1, list(intercept = intercept, beta = beta, log_nu = log_nu))
}


make_region_labels <- function(region_order,region_label) {
  tibble::tibble(
    Region = region_order,
    RegionLabel = region_label
  )
}

make_summary <- function(raw, region_order, region_labels_df) {
  area_meta <- raw %>%
    distinct(locus, area, Region, N, `Pfsa+`)

  region_meta <- area_meta %>%
    group_by(locus, Region) %>%
    summarise(
      N = sum(N, na.rm = TRUE),
      Pfsa_plus = sum(`Pfsa+`, na.rm = TRUE),
      .groups = "drop"
    )

  region_draws <- raw %>%
    filter(Region %in% region_order) %>%
    mutate(Region = factor(Region, levels = region_order))

  res_sum <- region_draws %>%
    group_by(locus, Region) %>%
    summarise(
      estimate       = median(slope, na.rm = TRUE),
      lower          = quantile(slope, 0.025, na.rm = TRUE),
      upper          = quantile(slope, 0.975, na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    left_join(region_meta, by = c("locus", "Region")) %>%
    left_join(region_labels_df, by = "Region") %>%
    mutate(
      estimate_pct   = 100 * estimate,
      lower_pct      = 100 * lower,
      upper_pct      = 100 * upper,
      N_lab          = scales::comma(N),
      df_lab         = sprintf("%.2f (%.2f-%.2f)", estimate, lower, upper),
      freq           = Pfsa_plus / N,
      freq_lab       = scales::percent(freq, accuracy = 0.1)
    )

  res_sum
}
