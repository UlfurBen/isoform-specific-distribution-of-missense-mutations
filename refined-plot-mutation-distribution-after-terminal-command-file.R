# Load necessary library
library(ggplot2)

# Define the path to your .txt file
file_path <- "/Users/ulfurfjolnisson/Downloads/output-3.txt"

# Read the data from the text file
data <- read.table(file_path, sep=",", header=FALSE)

# Extract the values from the second column
first_value <- data[1, 2]
second_value <- data[2, 2]
third_value <- data[3, 2]

# Perform the subtraction operation
modified_value <- first_value - second_value - third_value

# Save the modified value back to the data frame
data[1, 2] <- modified_value

# Optional: Write the modified data back to the text file
write.table(data, file_path, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)


# Read the data
data <- read.csv(file_path, header=FALSE)
new_column_headers <- c("Name", "Count")
colnames(data) <- new_column_headers

# Print the data to verify it's loaded correctly
print(data)

# Create a bar plot
p <- ggplot(data, aes(x = Name, y = Count)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Bar Plot of Counts by Name", x = "Name", y = "Count")

# Display the plot
print(p)
