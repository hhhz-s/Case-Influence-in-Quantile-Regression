pkgs <- c("dplyr", "tidyr", "ggplot2", "rgdal")
load_pkgs <- sapply(pkgs, require, character.only = TRUE)

wd <- "~/Documents/osu/Research/Case_weight_adjusted/Codes/Real_data_analysis/Housing_data" 
setwd(wd)
source("utilities_wls.R")
source("utilities_rda.R")

house_data <- read.csv("data/kc_house_data.csv")

## basic transformation
house_data <- house_data %>%
  drop_na() %>%
  mutate(date = substr(date, 1, 8),
         date = as.Date(date, format = "%Y%m%d"),
         sold_year = factor(format(date, "%Y")),
         basement = factor(ifelse(sqft_basement > 0, 1, 0)), 
         waterfront = factor(waterfront, levels = c(0,1)),
         condition = factor(case_when(
           condition <= 2 ~ "poor",
           condition == 3 ~ "average",
           TRUE ~ "good"
         ), levels = c("poor", "average", "good")),
         yr_built = factor(case_when(
           yr_built <= 1920 ~ "1900 - 1920",
           yr_built %in% 1921:1940 ~ "1921 - 1940",
           yr_built %in% 1941:1960 ~ "1941 - 1960",
           yr_built %in% 1961:1980 ~ "1961 - 1980",
           yr_built %in% 1981:2000 ~ "1981 - 2000",
           TRUE ~ "2001 - 2015"))) %>%
  select(!c(id, date, yr_renovated, sqft_lot15, sqft_above, sqft_basement)) %>%
  filter(bedrooms <= 10)

house2014 <- filter(house_data, sold_year == "2014") %>% select(-sold_year)
house2015 <- filter(house_data, sold_year == "2015") %>% select(-sold_year)

## location
zips <- readOGR(file.path(wd, "king_county_zip", "Zip_Codes.SHP"))
zips <- spTransform(zips, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
zips_df <- broom::tidy(zips) %>% 
  mutate(id = as.numeric(id)) %>%
  left_join(zips@data, by = c("id" = "OBJECTID")) %>%
  filter(ZIP %in% house_data$zipcode)

# add NICHE grade 
# https://www.niche.com/places-to-live/search/best-zip-codes-to-live/c/king-county-wa/
niche <- read.csv(file.path(wd, "niche_grade.csv"))

zips_df2 <- zips_df %>%
  left_join(niche, by = c("ZIP" = "zipcode")) %>%
  mutate(niche_overall_grade = case_when(
    niche_overall_grade %in% c("A+","A") ~ "zone 1",
    niche_overall_grade %in% c("A-", "B+") ~ "zone 2",
    TRUE ~ "zone 3"))

# boxplot of price vs. zone
house2014 %>%
  left_join(niche, by = "zipcode") %>%
  mutate(niche_overall_grade = case_when(
    niche_overall_grade %in% c("A+","A") ~ "zone 1",
    niche_overall_grade %in% c("A-", "B+") ~ "zone 2",
    TRUE ~ "zone 3")) %>%
  ggplot() +
  geom_boxplot(aes(x = niche_overall_grade, y = sqrt(price))) +
  theme_bw()

# Transformation
house2014_transf1 <- house2014 %>%
  mutate(sqrtPrice = sqrt(price),
         bedrooms = factor(case_when(
           bedrooms %in% c(1,2) ~ "1-2",
           bedrooms %in% 3:5 ~ "3-5",
           TRUE ~ "6 or more"
         )),
         sqrtSqft_living = sqrt(sqft_living),
         sqrtSqft_living15 = sqrt(sqft_living15)) %>%
  select(!c(sqft_living, sqft_living15, 
            sqft_lot, price, floors, lat, long))

relevant_features <- c("basement", "bathrooms", "bedrooms", "condition", 
                       "grade", "grade2", "sqrtSqft_living15", "sqrtSqft_living", 
                       "view", "waterfront", "yr_built", "zone")

# Final data 
house2014_transf <- house2014_transf1 %>%
  mutate(grade2 = grade^2,
         zone = factor(case_when(
           zipcode == 98039 ~ "zone 1",
           zipcode == 98004 ~ "zone 2",
           TRUE ~ "zone 3"))) %>%
  select(all_of(relevant_features), sqrtPrice)

xnames <- setdiff(colnames(house2014_transf), c("sqrtPrice"))
num_vars <- xnames[sapply(house2014_transf[,xnames], class) == "numeric"]
int_vars <- xnames[sapply(house2014_transf[,xnames], class) == "integer"]
cat_vars <- xnames[sapply(house2014_transf[,xnames], class) == "factor"]

# Process categorical variables
final_features <- c("basement_1", "bathrooms", "bedrooms_12", "bedrooms_35", "condition_good",
                    "condition_average", "waterfront1", "yr_built_2000", "yr_built_1980", "yr_built_1960",
                    "yr_built_1940", "yr_built_1920", "zone_2", "zone_1" ,"grade", "grade2", "sqrtSqft_living15", 
                    "sqrtSqft_living", "view")

Data <- house2014_transf %>%
  mutate(basement_1 = case_when(
           basement == 1 ~ 1,
           TRUE ~ 0),
         zone_1 = case_when(
           zone == "zone 1" ~ 1,
           TRUE ~ 0),
         zone_2 = case_when(
           zone == "zone 2" ~ 1,
           TRUE ~ 0),
         yr_built_1920 = case_when(
           yr_built == "1900 - 1920" ~ 1,
           TRUE ~ 0),
         yr_built_1940 = case_when(
           yr_built == "1921 - 1940" ~ 1,
           TRUE ~ 0),
         yr_built_1960 = case_when(
           yr_built == "1941 - 1960" ~ 1,
           TRUE ~ 0),
         yr_built_1980 = case_when(
           yr_built == "1961 - 1980" ~ 1,
           TRUE ~ 0),
         yr_built_2000 = case_when(
           yr_built == "1981 - 2000" ~ 1,
           TRUE ~ 0),
         waterfront1 = case_when(
           waterfront == 1 ~ 1,
           TRUE ~ 0),
         condition_good = case_when(
           condition == "good" ~ 1,
           TRUE ~ 0),
         condition_average = case_when(
           condition == "average" ~ 1,
           TRUE ~ 0),
         bedrooms_12 = case_when(
           bedrooms == "1-2" ~ 1,
           TRUE ~ 0),
         bedrooms_35 = case_when(
           bedrooms == "3-5" ~ 1,
           TRUE ~ 0)) %>%
  select(all_of(final_features), sqrtPrice)

Y = Data$sqrtPrice
X = Data[, 1:(ncol(Data)-1)]
save(Y, X, file = "Input.RData")
