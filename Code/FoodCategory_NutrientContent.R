#############################
#######AMSA Assignment#######
#############################

#By Antoine Ruzette - r0829308
#Applied Multivariate Statistical Anaysis [I0P16a]
#MSc in Bioinformatics, KU Leuven

#Use of data from the Australian Food Composition Database 
#to investigate the type of food and its content in nutrients






#load packages
library(naniar)
library(readxl)
library(dplyr)
library(cluster)
library(factoextra)
library(stringr)
library(tree)
library(ISLR)
library(dplyr)
library(ggplot2)
library(NbClust)
library(miscFuncs)
library(readr)
library(reshape2)
library(devtools)
#install_github("vqv/ggbiplot")
library('ggbiplot')
library(Rtsne)
library(randomForest)
library(rpart)

#set your specific environment
setwd("C:/Users/antoi/OneDrive/Documentos/Bioinformatics/M2/M2_Q1/AMSA/")

##1. Extraction, Transformation and Loading
#load the data from excel
nutrient.rawdata = read_excel("Release 1 - Food nutrient database.xlsx", 
                               ginsheet = 'All solids & liquids per 100g')
rawdata = as.data.frame(nutrient.rawdata)
nrow(rawdata)#1535 food items
ncol(rawdata)#252 nutrients information

#load the food categories (keys)
foodCategories = read_delim("FoodCategoriesConversion.txt",
                            delim = "\t", escape_double = FALSE, 
                            col_names = FALSE)

#remove all columns have more than 90% missing observations
missing_col = sapply(rawdata, function(x) sum(is.na(x))/1535)
missing_col = as.data.frame(missing_col)
colnames(missing_col)[1] = 'MissingPercentage'
#subset only the column having less than 10% missing data
nonmissing_col = subset(missing_col, missing_col[, 1] < 0.1)
nonmissing_col['colName'] = row.names(nonmissing_col)
nrow(nonmissing_col)
#63 on 252 columns that have less than 10% missing

#remove the highly missing variables
data_sub = rawdata[, nonmissing_col$colName]

#subsetting the variable space to the variables that we are interested in 
all.names = c('Public Food Key', 'Classification', 'Food Name', 'Energy, without dietary fibre', 'Moisture (water)', 'Protein', 'Nitrogen', 'Total Fat', 'Ash', 'Total dietary fibre', 
              'Alcohol', 'Total sugars', 'Starch', 'Iodine (I)', 'Iron (Fe)', 'Magnesium (Mg)', 'Phosphorus (P)', 'Potassium (K)', 'Vitamin A retinol equivalents', 
              'Sodium (Na)', 'Zinc (Zn)', 'Dietary folate equivalents', 'Vitamin E', 'Calcium (Ca)', 'Selenium (Se)')
features.names = c('Energy, without dietary fibre', 'Moisture (water)', 'Protein', 'Nitrogen', 'Total Fat', 'Ash', 'Total dietary fibre', 
                   'Alcohol', 'Total sugars', 'Starch', 'Iodine (I)', 'Iron (Fe)', 'Magnesium (Mg)', 'Phosphorus (P)', 'Potassium (K)', 'Vitamin A retinol equivalents', 
                   'Sodium (Na)', 'Zinc (Zn)', 'Dietary folate equivalents', 'Vitamin E', 'Calcium (Ca)', 'Selenium (Se)')
target.name = c('Classification')
data = data_sub[all.names]

#removing remaining missing data, if there still exist some
data = na.omit(data)
sum(is.na(data))#0 missing values
nrow(data)#1534 rows remaining
ncol(data)

data.x = data
X = select(data.x, - c("Public Food Key", "Classification", "Food Name"))
X = lapply(X,  as.numeric)
X = as.data.frame(X)

###2. Exploratory Analysis (correlations, k-means clustering, PCA and tSNE)
##correlation between variables
cormat = round(cor(X), 2)
cormat
melted_cormat <- melt(cormat)
#plot heatmap
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + 
  scale_x_discrete(guide = guide_axis(n.dodge=2))


#scaling data
X = as.data.frame(X)
X.scaled = scale(as.data.frame(X))
X.scaled = as.data.frame(X.scaled)

data.scaled = cbind(data$'Public Food Key', data$'Classification', data$'Food Name', X.scaled)
data = rename(data.scaled, replace = c(`data$"Public Food Key"` = "FoodKey", 
                                       `data$Classification` = 'Classification', 
                                       `data$"Food Name"` = 'FoodName'))
#create the 22 food categories
data$Classification = str_extract(as.character(data$Classification), "\\d{2}")
data$Classification = as.factor(data$Classification)
length(unique(data$Classification))#22 food categories

#distribution of the observations among classes
distr = aggregate(FoodKey ~ Classification,                                            # Count rows of all groups
          data = data,
          FUN = length)

##Principal Components Analysis
nutrient.pca <- prcomp(data[, c(4:25)], center = TRUE, scale = FALSE)
summary(nutrient.pca)#2 PCs explain 33% of the variance, 13 PCs explain 90%

#plot PC1 vs PC2
ggbiplot(nutrient.pca, ellipse = FALSE , groups = data$Classification, 
         var.axes = TRUE, varname.size = 3, alpha = 0.2, circle = TRUE) + 
  theme_classic() + 
  labs(color="Food type") + 
  ggtitle('Biplot of PCA')

##tSNE 
colors = rainbow(length(unique(data$Classification)))
names(colors) = unique(data$Classification)

tsne <- Rtsne(data[, c(4:25)], check_duplicates = FALSE, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
tsne.Y = as.data.frame(tsne$Y)
#plot t-SNE1 vs t-SNE2
ggplot(tsne.Y, aes(V1, V2, colour = data$Classification)) + 
  geom_point() + 
  theme_classic() + 
  labs(color = 'Food type') +
  ggtitle('tSNE exploration')
  

##k-means clustering

#explore the data scatterplot
set.seed(2)
best_clust = NbClust(data = data[, c(4:25)], distance = "euclidean", min.nc = 20, max.nc = 25, method = 'kmeans')
fviz_nbclust(best_clust, ggtheme = theme_minimal())
km = kmeans(data[, c(4:25)], centers = 22, nstart = 25, iter.max = 1000)
#plot clusters (with ellipse) over the data scatterplot
fviz_cluster(object=km, data = data[, c(4:25)], ellipse.type = 'norm', geom = 'point', main = 'K-Means clustering', ggtheme = theme_minimal())

#table to compare the clustering patterns and the real food categories
data$km_cluster = as.numeric(km$cluster)
assignment_table = table(data$km_cluster, data$Classification)
latextable(foodCategories)
aggregate(cbind(Protein, Nitrogen, Ash, Alcohol) ~ km_cluster, data = data, FUN=mean)



##DECISION TREE MODELLING
features.names = c("Classification", "Energy..without.dietary.fibre", "Moisture..water.", "Protein", "Nitrogen", "Total.Fat", 
                   "Ash", "Total.dietary.fibre", "Alcohol", "Total.sugars", "Starch", "Iodine..I.", "Iron..Fe.", 
                   "Magnesium..Mg.", "Phosphorus..P.", "Potassium..K.", "Vitamin.A.retinol.equivalents", "Sodium..Na.", 
                   "Zinc..Zn.", "Dietary.folate.equivalents", "Vitamin.E", "Calcium..Ca.", "Selenium..Se.")

##3. Tree based Modeling
#split in train and test sets (approx. 75/25 split)
set.seed(3)
train = data %>% sample_n(1200)
X_train = train[features.names]
y_train = train[target.name]

test = data %>% setdiff(train)
X_test = test[features.names]
y_test = test[target.name]

#SINGLE TREE 
tree_nutrient = tree(Classification ~ . , X_train)
summary(tree_nutrient)

single_pred = predict(tree_nutrient, X_test, type = "class")
mean(single_pred == X_test$Classification)
#train accuracy = 1 - 0.2439 = 0.7561 = 75.6%
#number of terminal nodes = 20 (compared to 22 food categories in the beginning)

pdf("SingleTree.pdf", width = 8, height = 7, 
    bg = "white", colormodel = "cmyk", 
    paper = "A4")  
plot(tree_nutrient)
text(tree_nutrient, pretty = 2, cex = 0.45)
dev.off() 

#test accuracy
tree_pred = predict(tree_nutrient, X_test, type = "class")
mean(tree_pred == X_test$Classification)
#75.4% of test accuracy

#SINGLE PRUNED TREE
set.seed(4)
cv_nutrient = cv.tree(tree_nutrient, FUN = prune.misclass)
plot(cv_nutrient$size, cv_nutrient$dev, type = "b")
prune_nutrient = prune.misclass(tree_nutrient, best = 7)
summary(prune_nutrient)
#train accuracy = 0.635 
#reduced accuracy for the stake of interpretability


pdf("PrunedTree.pdf", width = 8, height = 7, 
    bg = "white", colormodel = "cmyk", 
    paper = "A4")  
plot(prune_nutrient)
text(prune_nutrient, pretty = 0, cex = 0.6)
dev.off() 

prune_pred = predict(prune_nutrient, X_test, type = "class")
mean(prune_pred == X_test$Classification)
#65% of test accuracy

#RANDOM FOREST
set.seed(5)
rf_nutrient = randomForest(Classification~., 
                           data = X_train, 
                           importance = TRUE)
#print summary 
rf_nutrient
#OOB accuracy = 90.75% 

rf_estimate = predict(rf_nutrient, newdata = X_test)
mean(rf_estimate == X_test$Classification)
#test accuracy = 94%

ggplot() + 
  geom_point(aes(x = X_test$Classification, y = rf_estimate)) +
  geom_abline()

#calculate importance of the used variables
importance(rf_nutrient)
#plot importance measure of the used variables ()
varImpPlot(rf_nutrient)


