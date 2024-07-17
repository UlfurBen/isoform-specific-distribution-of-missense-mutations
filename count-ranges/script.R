# Load necessary library
library(dplyr)

# Read the data
file_path <- "enrichment-2.txt"
data <- read.csv(file_path, header = TRUE)

# Filter out rows where Variation_count is 0
filtered_data <- data %>% filter(Variation_count != 0)

# Define the ranges
ranges <- list(
  "1-5" = c(1, 5),
  "6-10" = c(6, 10),
  "11-20" = c(11, 20),
  "21-30" = c(21, 30),
  "31-40" = c(31, 40),
  "41-50" = c(41, 50),
  "51-100" = c(51, 100),
  "101-150" = c(101, 150),
  "151-200" = c(151, 200),
  "201-1000" = c(201, 1000),
  ">1000" = c(1001, Inf)
)

# Function to count the number of rows in each range
count_ranges <- function(data, ranges) {
  counts <- sapply(ranges, function(range) {
    sum(data$Variation_count >= range[1] & data$Variation_count <= range[2])
  })
  return(counts)
}

# Get counts for each range
range_counts <- count_ranges(filtered_data, ranges)

# Print the results
print(range_counts)

# Optional: Write the results to a file
write.csv(data.frame(Range = names(range_counts), Count = range_counts), "range_counts.csv", row.names = FALSE)
