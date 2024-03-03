print("Start.")

file <- "data/data_fig5.csv"
data <- read.csv(file, sep=";")
print(dim(data))
flat_data <- array(data, dim=c(prod(dim(data))))
x = array(c(c(1, 2), c(3, 4)), dim=c(4))
print(prod(dim(data)))
print(length(data.flat))
print(class(data.flat))
#print(data.flat)
print(class(flat_data))
# hist(data.flat)

print("Done.")