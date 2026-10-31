library(data.table)
print("Start.")

m0 = 0.5;
a = 0.5;
mu = 2.0;
lambda = 0.0;

file <- "data/data_fig5_check.csv";
data <- as.matrix(read.csv(file, sep=";"));
print(class(data));

shift(data.table(data), 1);
data_shift <- t(apply(data, 1, function(x) shift(as.data.table(x), n = 1)$V1))
print(dim(as.matrix(data_shift)));

x_action <- m0*0.5*(data-data_shift)/a+a*(mu*data**2*0.5+lambda*data**4);
action = sum(x_action);



data_flat = as.vector(data);



# hist(data_flat[data_flat < 200 & data_flat > -200]);

print("Done.")