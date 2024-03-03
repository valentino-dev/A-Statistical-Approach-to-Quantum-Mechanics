print("Start.")


file <- "data/data_fig5_check.csv"
data <- read.csv(file, sep=";")
data_flat = as.vector(t(data[100:5000,]))
hist(data_flat)

print("Done.")