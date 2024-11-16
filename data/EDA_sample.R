rm(list=ls())

pkgs <- c("dplyr", "ggplot2", "tidyr","gridExtra")
load_pkgs <- sapply(pkgs, require, character.only = TRUE)

# specify the working dir
setwd("path/to/github/repo") 

# Load data
load("data/Input.RData")

set.seed(1)
X = as.matrix(X)
n_total = nrow(X)
p = ncol(X)

n = 200
subsample = sample(n_total, 2*n)
train_idx = subsample[1:n]
test_idx = subsample[(n+1):(2*n)]

X_train = X[train_idx,]
Y_train = Y[train_idx]

df <- data.frame(X_train,Y=Y_train)

table(df$bedroom_category)
table(df$view)
table(df$waterfront)
table(df$zone)
table(df$condition)
## EDA

# violin plot of Y ~ basement

basement_plot <- ggplot(df, aes(factor(basement_1), Y)) +
  #geom_violin(trim=FALSE, fill="gray") +
  geom_violin(trim=FALSE) +
  labs(x = "Basement", y = "House price") +
  #geom_boxplot(width=0.05)+
  theme_minimal()

# violin plot of Y ~ bedroom

df$bedroom_category <- ifelse(df$bedrooms_12 == 1, "1-2",
                              ifelse(df$bedrooms_35 == 1, "3-5", ">=6"))

df$bedroom_category <- factor(df$bedroom_category,levels = c("1-2","3-5",">=6"))

bedroom_plot <- ggplot(df, aes(x = bedroom_category, y = Y)) +
  geom_violin() +
  labs(x = "Number of Bedrooms", y = "House price") +
  theme_minimal()


# violin plot of y~condition
df$condition <- ifelse(df$condition_average==1,'average',
                       ifelse(df$condition_good==1,'good','poor'))

df$condition <-factor(df$condition,levels = c('poor','average','good'))

condition_plot <- ggplot(df, aes(x = condition, y = Y)) +
  geom_violin() +
  labs(x = "House condition", y = "House price") +
  theme_minimal() +
  scale_x_discrete(drop = FALSE)


# violin plot of Y~year built 
df <- df %>%
  mutate(year_built_category = case_when(
    yr_built_1920 == 1 ~ "1920", 
    yr_built_1940 == 1 ~ "1940",
    yr_built_1960 == 1 ~ "1960",
    yr_built_1980 == 1 ~ "1980",
    yr_built_2000 == 1 ~ "2000",
    TRUE ~ "2015"
  ))

df$year_built_category <- factor(df$year_built_category, levels = c("1920", "1940", "1960", "1980", "2000", "2015"))


yr_built_plot <- ggplot(df, aes(x = year_built_category, y = Y)) +
  geom_violin() +
  labs(x = "Year built", y = "House price") +
  theme_minimal()

# violin plot of Y~zone 
df <- df%>%
  mutate(zone = case_when(
    zone_1==1 ~ 'zone_1',
    zone_2==1 ~ 'zone_2',
    TRUE ~ 'zone_3'))

df$zone <- factor(df$zone,levels = c('zone_1','zone_2','zone_3'))

zone_plot <- ggplot(df,aes(x = zone, y =Y)) +
  geom_violin() +
  labs(x = "Zone", y = "House price") +
  theme_minimal()


# violin plot of Y~waterfront
df <- df%>%
  mutate(waterfront = case_when(
    waterfront1==1 ~ 'yes',
    TRUE ~ 'no'))

df$waterfront <- factor(df$waterfront,levels = c('no','yes'))

waterfront_plot <- ggplot(df,aes(x = waterfront, y =Y)) +
  geom_violin() +
  labs(x = "Waterfront", y = "House price") +
  theme_minimal()

# scatter plot of Y~view 

## local reference line cannot be added 

#ggplot(df, aes(x = view, y = Y)) +
#  geom_point() +  
#  geom_smooth(method = "loess", se = FALSE, color = '#619CFF',linetype='dashed',size=1) +
#  labs(x = "View", y = "House price") +
#  theme_minimal()

view_scatter_plot <- ggplot(df, aes(x = view, y = Y)) +
  geom_point(alpha=0.25) +
  geom_smooth(method = "lm", se = FALSE, color = '#619CFF',linetype='dashed',size=1) + 
  labs(x = "View", y = "House Price") +
  theme_minimal()

# violin plot of Y ~ view 
view_violin_plot <- ggplot(df, aes(x = as.factor(view), y = Y)) +
  geom_violin() +
  labs(x = "View", y = "House price") +
  theme_minimal()


# scatter plot of Y ~ bathrooms 

bathroom_plot <- ggplot(df, aes(x = bathrooms, y = Y)) +
  geom_point(alpha=0.25) +
  geom_smooth(method = "loess", se = FALSE, color = '#619CFF',linetype='dashed',size=1) + 
  labs(x = "Number of Bathrooms", y = "House price") +
  theme_minimal()

# scatter plot of Y ~ grade 
grade_plot <-ggplot(df, aes(x = grade, y = Y)) +
  geom_point(alpha=0.25) +  # Scatter plot
  geom_smooth(method = "loess", se = FALSE, color = '#619CFF',linetype='dashed',size=1) + 
  labs(x = "Grade", y = "House price") +
  theme_minimal()

# scatter plot of Y ~ grade2
grade2_plot <- ggplot(df, aes(x = grade2, y = Y)) +
  geom_point(alpha=0.25) +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = '#619CFF',linetype='dashed',size=1) + 
  labs(x = "Grade Squared", y = "House price") +
  theme_minimal()

# scatter plot of Y~sqrtSqft_living15
living15_plot <- ggplot(df, aes(x = sqrtSqft_living15, y = Y)) +
  geom_point(alpha=0.25) +  
  geom_smooth(method = "lm", se = FALSE, color = '#619CFF',linetype='dashed',size=1) + 
  labs(x = "sqrtSqft_living15", y = "House price") +
  theme_minimal()

# scatter plot of Y~sqrtSqft_living
living_plot <- ggplot(df, aes(x = sqrtSqft_living, y = Y)) +
  geom_point(alpha=0.25) +  #
  geom_smooth(method = "lm", se = FALSE, color = '#619CFF',linetype='dashed',size=1) +   
  labs(x = "sqrtSqft_living", y = "House price") +
  theme_minimal()



plot_list_v1 <- list(basement_plot,bedroom_plot,condition_plot,
                     yr_built_plot,zone_plot,waterfront_plot,
                     view_violin_plot,bathroom_plot,grade_plot,
                     grade2_plot,living15_plot,living_plot)

grid_arranged_v1 <- grid.arrange(grobs = plot_list_v1, ncol = 3, nrow = 4)

# treat view as numeric 

plot_list_v2 <- list(basement_plot,bedroom_plot,condition_plot,
                     yr_built_plot,zone_plot,waterfront_plot,
                     view_scatter_plot,bathroom_plot,grade_plot,
                     grade2_plot,living15_plot,living_plot)

# save as 12 * 12 final EDA plot
grid_arranged_v2 <- grid.arrange(grobs = plot_list_v2, ncol = 3, nrow = 4)
