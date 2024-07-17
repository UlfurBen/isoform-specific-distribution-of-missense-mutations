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

# Function to parse the data into a named vector
parse_data <- function(data) {
  lines <- strsplit(data, "\n")[[1]]
  lines <- lines[lines != ""] # Remove empty lines
  data_list <- lapply(lines, function(line) {
    parts <- strsplit(line, " ")[[1]]
    value <- as.numeric(parts[length(parts)])
    name <- paste(parts[1:(length(parts) - 1)], collapse = " ")
    return(c(name, value))
  })
  data_vector <- as.numeric(sapply(data_list, function(x) x[2]))
  names(data_vector) <- sapply(data_list, function(x) x[1])
  return(data_vector)
}

# Parsing EM and nonEM data
EM_vector <- parse_data(EM_data)
nonEM_vector <- parse_data(nonEM_data)

# Function to create the contingency table and perform Fisher's exact test
fisher_test <- function(EM_vector, nonEM_vector, change) {
  EM_count <- EM_vector[change]
  nonEM_count <- nonEM_vector[change]
  
  EM_total <- sum(EM_vector)
  nonEM_total <- sum(nonEM_vector)
  
  EM_other <- EM_total - EM_count
  nonEM_other <- nonEM_total - nonEM_count
  
  table <- matrix(c(EM_count, nonEM_count, EM_other, nonEM_other), nrow = 2)
  
  p_value <- fisher.test(table)$p.value
  
  return(list(table = table, p_value = p_value))
}

# List of categorical changes
changes <- names(EM_vector)

# Creating tables and performing Fisher's exact test
results <- lapply(changes, function(change) {
  fisher_test(EM_vector, nonEM_vector, change)
})

# Naming the results list
names(results) <- changes

# Creating a data frame to display the results
results_df <- data.frame(
  Change = changes,
  EM_Count = sapply(results, function(x) x$table[1, 1]),
  nonEM_Count = sapply(results, function(x) x$table[1, 2]),
  EM_Other = sapply(results, function(x) x$table[2, 1]),
  nonEM_Other = sapply(results, function(x) x$table[2, 2]),
  P_Value = sapply(results, function(x) format(x$p_value, digits = 6, scientific = FALSE))
)

# Print the results without scientific notation and limit to 6 digits
options(scipen = 999)  # Turn off scientific notation
print(results_df)
