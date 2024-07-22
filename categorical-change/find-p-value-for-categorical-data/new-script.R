# Provided data as character vectors
EM_data <- "
Nonpolar->Nonpolar 349
Nonpolar->Polar 459
Nonpolar->Negatively_charged 123
Nonpolar->Positively_charged 133
Polar->Nonpolar 261
Polar->Polar 189
Polar->Negatively_charged 432
Polar->Positively_charged 467
Negatively_charged->Nonpolar 101
Negatively_charged->Polar 92
Negatively_charged->Negatively_charged 42
Negatively_charged->Positively_charged 93
Positively_charged->Nonpolar 394
Positively_charged->Polar 143
Positively_charged->Negatively_charged 35
Positively_charged->Positively_charged 49
"

nonEM_data <- "
Nonpolar->Nonpolar 7427
Nonpolar->Polar 3926
Nonpolar->Negatively_charged 1607
Nonpolar->Positively_charged 2944
Polar->Nonpolar 2765
Polar->Polar 2750
Polar->Negatively_charged 433
Polar->Positively_charged 1796
Negatively_charged->Nonpolar 1058
Negatively_charged->Polar 977
Negatively_charged->Negatively_charged 289
Negatively_charged->Positively_charged 1115
Positively_charged->Nonpolar 1791
Positively_charged->Polar 2405
Positively_charged->Negatively_charged 394
Positively_charged->Positively_charged 1046
"

contingency_table <- data.frame(
  "EM genes" = c(432, 45),
  "non-EM genes" = c(433, 797),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )
