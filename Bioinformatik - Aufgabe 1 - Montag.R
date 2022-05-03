#section 1 assessment 

#Q1: what are the two solutions to 2x^2 -x -4 = 0? Use the quadratic equasion

# assigning values to variables
a < -2
b <- -1
c <- -4

# solving the quadratic equation
(-b + sqrt(b^2 - 4*a*c))/(2*a)
(-b - sqrt(b^2 - 4*a*c))/(2*a)

#Q2: Use R to compute log base 4 of 1024. You can use the help() function to learn how to use arguments to change the base of the log() function.

log(1024, 4)

#Q3
library(dslabs)
data(movielens)
#Q3a: How many rows are in the dataset?
str(movielens)

#Q3b: s.Q3a

#Q4: We already know we can use the levels() function to determine the levels of a factor. 
#A different function, nlevels(), may be used to determine the number of levels of a factor.´
#Use this function to determine how many levels are in the factor genres in the movielens data frame.

nlevels(movielens$genres)

#section 2 assessment 

name <- c("Mandi", "Amy", "Nicole", "Olivia")
distance <- c(0.8, 3.1, 2.8, 4.0)
time <- c(10, 30, 40, 50)

#Write a line of code to convert time to hours. Remember there are 60 minutes in an hour. 

time_h <- time/60

#Then write a line of code to calculate the speed of each runner in miles per hour. Speed is distance divided by time.

speed <- distance/time_h

#section 3 assessment 

#load packages and prepare the data - heights dataset
library(tidyverse)
library(dslabs)
library(data.table)
data(heights)
options(digits = 3)
#above average height 
average_height <- mean(heights$height)
ind <- heights$height > mean(average_height)
sum(ind)
#above average height and female 
female_sample <- heights$sex == "Female"
taller_female <- heights$height > mean(average_height)
indf <- female_sample & taller_female
sum(indf)
#what proportion of individuals in the dataset are females?
mean(heights$sex == "Female")
#minimum height
min(heights$height)
#Use the match() function to determine the index of the first individual with the minimum height.
match(50, heights$height)
#Subset the sex column of the dataset by the index in 5 to determine the individual’s sex.
heights$sex[match(50,heights$height)]
#max height
max(heights$height)
#Write code to create a vector x that includes the integers between the minimum and maximum heights in this dataset (as numbers).
max(heights$height)
min(heights$height)
x <- 50:83
#How many of the integers in x are NOT heights in the dataset?
sum(!(x %in% heights$height))
#Using the heights dataset, create a new column of heights in centimeters named ht_cm. Recall that 1 inch = 2.54 centimeters. Save the resulting dataset as heights2.
heights2 <- mutate(heights, ht_cm = height*2.54)
head(heights2)
#What is the height in centimeters of the 18th individual (index 18)?
heights2$ht_cm[18]
#What is the mean height in centimeters?
mean(heights2$ht_cm)
#How many females are in the heights2 dataset?
sum(female_sample2 <- heights2$sex == "Female")
#or
female_sample2 <- filter(heights2, sex == "Female")
nrow(female_sample2)
#What is the mean height of the females in centimeters?
mean(female_sample2$ht_cm)

#Q8
library(dslabs)
data(olive)
head(olive)

#Plot the percent palmitic acid versus palmitoleic acid in a scatterplot. What relationship do you see?
palmitic_acid <- olive$palmitic
palmitoleic_acid <- olive$palmitoleic
plot(palmitic_acid,palmitoleic_acid)

#Create a histogram of the percentage of eicosenoic acid in olive.
hist(olive$eicosenoic)

#Make a boxplot of palmitic acid percentage in olive with separate distributions for each region.
boxplot(palmitic ~ region, data = olive)


#DataCamp
#section 1 DataCamp

# Here is how you compute the sum for the first 20 integers
20*(20+1)/2 

# However, we can define a variable to use the formula for other values of n
n <- 20
n*(n+1)/2

n <- 25
n*(n+1)/2

# Below, write code to calculate the sum of the first 100 integers
n <- 100
n*(n+1)/2

# Below, write code to calculate the sum of the first 1000 integers 
n <- 1000
n*(n+1)/2

# log to the base 2 
log2(16)

# sqrt of the log to the base 2 of 16:
sqrt(log2(16))

# Compute log to the base 10 (log10) of the sqrt of 100. Do not use variables.
sqrt(log10(100))

#section 2 DataCamp

# Load package and data

library(dslabs)
data(murders)

# Use the function names to extract the variable names 
names(murders)

# To access the population variable from the murders dataset use this code:
p <- murders$population 

# To determine the class of object `p` we use this code:
class(p)

# Use the accessor to extract state abbreviations and assign it to a

a <- murders$abb

# Determine the class of a
class(a)

# We extract the population like this:
p <- murders$population

# This is how we do the same with the square brackets:
o <- murders[["population"]] 

# We can confirm these two are the same
identical(o, p)

# Use square brackets to extract `abb` from `murders` and assign it to b
b <- murders[["abb"]]
# Check if `a` and `b` are identical 
identical(a, b)

# We can see the class of the region variable using class
class(murders$region)

# Determine the number of regions included in this variable 
length(murders$region)
levels(murders$region)

length(levels(murders$region))

# The "c" in `c()` is short for "concatenate," which is the action of connecting items into a chain
# The function `c()` connects all of the strings within it into a single vector, which we can assign to `x`
x <- c("a", "a", "b", "b", "b", "c")
# Here is an example of what the table function does
table(x)

#section 3 DataCamp 

# Here is an example creating a numeric vector named cost
cost <- c(50, 75, 90, 100, 150)

# Create a numeric vector to store the temperatures listed in the instructions into a vector named temp
# Make sure to follow the same order in the instructions
temp <- c(35, 88, 42, 84, 81, 30)

# here is an example of how to create a character vector
food <- c("pizza", "burgers", "salads", "cheese", "pasta")

# Create a character vector called city to store the city names
# Make sure to follow the same order as in the instructions
city <- c("Beijing", "Lagos", "Paris", "Rio de Janeiro", "San Juan", "Toronto")

# Associate the cost values with its corresponding food item
cost <- c(50, 75, 90, 100, 150)
food <- c("pizza", "burgers", "salads", "cheese", "pasta")
names(cost) <- food

# You already wrote this code
temp <- c(35, 88, 42, 84, 81, 30)
city <- c("Beijing", "Lagos", "Paris", "Rio de Janeiro", "San Juan", "Toronto")

# Associate the temperature values with its corresponding city
names(temp) <- city

# cost of the last 3 items in our food list:
cost[3:5]

# temperatures of the first three cities in the list:
temp[1:3]

# Access the cost of pizza and pasta from our food list 
cost[c(1,5)]

# Define temp
temp <- c(35, 88, 42, 84, 81, 30)
city <- c("Beijing", "Lagos", "Paris", "Rio de Janeiro", "San Juan", "Toronto")
names(temp) <- city

# Access the temperatures of Paris and San Juan
temp[c(3,5)]

# Create a vector m of integers that starts at 32 and ends at 99.
m <- 32:99

# Determine the length of object m.
length(m)

# Create a vector x of integers that starts at 12 and ends at 73.
x <- 12:73
# Determine the length of object x.
length(x)

# Create a vector with the multiples of 7, smaller than 50.
seq(7, 49, 7) 

# Create a vector containing all the positive odd numbers smaller than 100.
# The numbers should be in ascending order
seq(1, 100, 2)

# We can create a vector with the multiples of 7, smaller than 50 like this 
seq(7, 49, 7) 

# But note that the second argument does not need to be the last number
# It simply determines the maximum value permitted
# so the following line of code produces the same vector as seq(7, 49, 7)
seq(7, 50, 7)

# Create a sequence of numbers from 6 to 55, with 4/7 increments and determine its length
length(seq(6, 55, 4/7))

# Store the sequence in the object a
a <- seq(1, 10, length.out = 100)

# Determine the class of a
class(a)

# Store the sequence in the object a
a <- seq(1, 10)

# Determine the class of a
class(a)

# Check the class of 1, assigned to the object a
class(1)

# Confirm the class of 1L is integer
class(1L)

# Define the vector x
x <- c(1, 3, 5,"a")

# Note that the x is character vector
x 

# Redefine `x` to typecast it to get an integer vector using `as.numeric`.
# You will get a warning but that is okay
x <- as.numeric(x)

#section 4 DataCamp

# Access the `state` variable and store it in an object 
states <- murders$state 

# Sort the object alphabetically and redefine the object 
states <- sort(states) 

# Report the first alphabetical value  
states[1]

# Access population values from the dataset and store it in pop
pop <- murders$population
# Sort the object and save it in the same object 
pop <- sort(pop)
# Report the smallest population size 
pop[1]

# Access population from the dataset and store it in pop
pop <- murders$population
# Use the command order to find the vector of indexes that order pop and store in object ord
ord <- order(pop)
# Find the index number of the entry with the smallest population size
ord

# Find the index of the smallest value for variable total 
which.min(murders$total)

# Find the index of the smallest value for population
which.min(murders$population)

# Define the variable i to be the index of the smallest state
i <- which.min(murders$population)

# Define variable states to hold the states
states <- murders$state
# Use the index you just defined to find the state with the smallest population
states[i]

# Store temperatures in an object 
temp <- c(35, 88, 42, 84, 81, 30)

# Store city names in an object 
city <- c("Beijing", "Lagos", "Paris", "Rio de Janeiro", "San Juan", "Toronto")

# Create data frame with city names and temperature 
city_temps <- data.frame(name = city, temperature = temp)

# Define a variable states to be the state names 
states <- murders$state

# Define a variable ranks to determine the population size ranks 
ranks <- rank(murders$population)

# Create a data frame my_df with the state name and its rank
my_df <- data.frame(name = states, population = ranks)

my_df

# Define a variable states to be the state names from the murders data frame
states <- murders$state

# Define a variable ranks to determine the population size ranks 
ranks <- rank(murders$population)

# Define a variable ind to store the indexes needed to order the population values
ind <- order(ranks)

# Create a data frame my_df with the state name and its rank and ordered from least populous to most 

my_df <- data.frame(name = states[ind], ranks = ind[ranks])
my_df

# Using new dataset 
library(dslabs)
data(na_example)

# Checking the structure 
str(na_example)

# Find out the mean of the entire dataset 
mean(na_example)

# Use is.na to create a logical index ind that tells which entries are NA
ind <- is.na(na_example)
# Determine how many NA ind has using the sum function
sum(ind)

# Note what we can do with the ! operator
x <- c(1, 2, 3)
ind <- c(FALSE, TRUE, FALSE)
x[!ind]

# Create the ind vector
library(dslabs)
data(na_example)
ind <- is.na(na_example)

# We saw that this gives an NA
mean(na_example)

# Compute the average, for entries of na_example that are not NA 
mean(na_example[!ind])

#section 5 DataCamp

# Assign city names to `city` 
city <- c("Beijing", "Lagos", "Paris", "Rio de Janeiro", "San Juan", "Toronto")

# Store temperature values in `temp`
temp <- c(35, 88, 42, 84, 81, 30)

# Convert temperature into Celsius and overwrite the original values of 'temp' with these Celsius values
temp <- 5/9*(temp-32)
# Create a data frame `city_temps` 
city_temps <- data.frame(name = city, temperature = temp)

# Define an object `x` with the numbers 1 through 100
x <- seq(1, 100)
# Compute the sum 
sum(1/x^2)

# Load the data
library(dslabs)
data(murders)

# Store the per 100,000 murder rate for each state in murder_rate
murder_rate <- murders$total/murders$population*100000
# Calculate the average murder rate in the US 
mean(murder_rate)

#section 6 assessment DataCamp

# Store the murder rate per 100,000 for each state, in `murder_rate`
murder_rate <- murders$total / murders$population * 100000

# Store the `murder_rate < 1` in `low` 
low <- murder_rate < 1

# Store the murder rate per 100,000 for each state, in murder_rate
murder_rate <- murders$total/murders$population*100000

# Store the murder_rate < 1 in low 
low <- murder_rate < 1

# Names of states with murder rates lower than 1
murders$state[low]

# Store the murder rate per 100,000 for each state, in `murder_rate`
murder_rate <- murders$total/murders$population*100000

# Store the `murder_rate < 1` in `low` 
low <- murder_rate < 1

# Create a vector ind for states in the Northeast and with murder rates lower than 1. 
ind <- low & murders$region=="Northeast"
# Names of states in `ind` 
murders$state[ind]

# Store the murder rate per 100,000 for each state, in murder_rate
murder_rate <- murders$total/murders$population*100000


# Compute the average murder rate using `mean` and store it in object named `avg`
avg <- mean(murder_rate)

# How many states have murder rates below avg ? Check using sum 
sum(murder_rate < avg)

# Store the 3 abbreviations in a vector called `abbs` (remember that they are character vectors and need quotes)
abbs <- c("AK", "MI", "IA")
# Match the abbs to the murders$abb and store in ind
ind <- match(abbs, murders$abb)
# Print state names from ind
murders$state[ind]

# Store the 5 abbreviations in `abbs`. (remember that they are character vectors)
abbs <- c("MA", "ME", "MI", "MO", "MU")

# Use the %in% command to check if the entries of abbs are abbreviations in the the murders data frame
abbs%in%murders$abb

# Store the 5 abbreviations in abbs. (remember that they are character vectors)
abbs <- c("MA", "ME", "MI", "MO", "MU") 

# Use the `which` command and `!` operator to find out which index abbreviations are not actually part of the dataset and store in `ind`
ind <- which(!abbs%in%murders$abb)
# Names of abbreviations in `ind`
abbs[ind]

#section 7 DataCamp 

# Loading data
library(dslabs)
data(murders)

# Loading dplyr
library(dplyr)

# Redefine murders so that it includes a column named rate with the per 100,000 murder rates
murders <- mutate(murders, rate = total / population * 100000)
murders

# Note that if you want ranks from highest to lowest you can take the negative and then compute the ranks 
x <- c(88, 100, 83, 92, 94)
rank(-x)

# Defining rate
rate <-  murders$total/ murders$population * 100000

# Redefine murders to include a column named rank
# with the ranks of rate from highest to lowest
murders <- mutate(murders, rank(-rate))
murders

# Load dplyr
library(dplyr)

# Use select to only show state names and abbreviations from murders
select(murders, state, abb)

# Add the necessary columns
murders <- mutate(murders, rate = total/population * 100000, rank = rank(-rate))

# Filter to show the top 5 states with the highest murder rates
filter(murders, rank <=5)

# Use filter to create a new data frame no_south
no_south <- filter(murders, region != "South")
# Use nrow() to calculate the number of rows
nrow(no_south)

# Create a new data frame called murders_nw with only the states from the northeast and the west
murders_nw <- filter(murders, region %in% c("Northeast", "West"))
# Number of states (rows) in this category 
nrow(murders_nw)

# add the rate column
murders <- mutate(murders, rate =  total / population * 100000, rank = rank(-rate))

# Create a table, call it my_states, that satisfies both the conditions 
my_states <- filter(murders, region %in% c("Northeast", "West") & rate < 1)
# Use select to show only the state name, the murder rate and the rank
select(my_states, state, rate, rank)

# Load library
library(dplyr)

## Define the rate column
murders <- mutate(murders, rate =  total / population * 100000, rank = rank(-rate))

# show the result and only include the state, rate, and rank columns, all in one line, in that order
filter(murders, region %in% c("Northeast", "West") & rate < 1) %>% select(state, rate, rank)

# Loading the libraries
library(dplyr)
data(murders)

# Create new data frame called my_states (with specifications in the instructions)
my_states <- murders %>%
  mutate(rate = total / population * 100000, rank = rank(-rate)) %>%
  filter(region %in% c("Northeast", "West") & rate < 1) %>%
  select(state, rate, rank)
my_states

#section 8 DataCamp 

# Load the datasets and define some variables
library(dslabs)
data(murders)

population_in_millions <- murders$population/10^6
total_gun_murders <- murders$total

plot(population_in_millions, total_gun_murders)

# Transform population (not population in millions) using the log10 transformation and save to object log10_population
log10_population <- log(murders$population, 10)

# Transform total gun murders using log10 transformation and save to object log10_total_gun_murders
log10_total_gun_murders <- log(murders$total, 10)

# Create a scatterplot with the log scale transformed population and murders 
plot(log10_population, log10_total_gun_murders)

# Store the population in millions and save to population_in_millions 
population_in_millions <- murders$population/10^6


# Create a histogram of this variable
hist(population_in_millions)

# Create a boxplot of state populations by region for the murders dataset
boxplot(population~region, data = murders)

#section 9 DataCamp 

# Assign the state abbreviation when the state name is longer than 8 characters 

new_names <- ifelse(nchar(murders$state)>8, murders$abb, murders$state)
head(new_names)

# Create function called `sum_n`
sum_n <- function(n){
  x <- 1:n
  sum(x)
}
# Use the function to determine the sum of integers from 1 to 5000
sum_n(5000)

# Create `altman_plot` 
altman_plot <- function(x, y){
  plot(x + y, y - x)
}

# Run this code 
x <- 3
my_func <- function(y){
  x <- 5
  y+5
}

# Print the value of x 
my_func(x)

# Here is an example of a function that adds numbers from 1 to n
example_func <- function(n){
  x <- 1:n
  sum(x)
}

# Here is the sum of the first 100 numbers
example_func(100)

# Write a function compute_s_n with argument n that for any given n computes the sum of 1 + 2^2 + ...+ n^2
compute_s_n <- function(n){
  x <- 1:n
  sum(x^2)
}

# Report the value of the sum when n=10
compute_s_n(10)

# Define a function and store it in `compute_s_n`
compute_s_n <- function(n){
  x <- 1:n
  sum(x^2)
}

# Create a vector for storing results
s_n <- vector("numeric", 25)

# write a for-loop to store the results in s_n
for(i in 1:25){
  s_n[i] <- compute_s_n(i)
}

# Define the function
compute_s_n <- function(n){
  x <- 1:n
  sum(x^2)
}

# Define the vector of n
n <- 1:25

# Define the vector to store data
s_n <- vector("numeric", 25)
for(i in n){
  s_n[i] <- compute_s_n(i)
}

#  Create the plot 
plot(n, s_n)

# Define the function
compute_s_n <- function(n){
  x <- 1:n
  sum(x^2)
}

# Define the vector of n
n <- 1:25

# Define the vector to store data
s_n <- vector("numeric", 25)
for(i in n){
  s_n[i] <- compute_s_n(i)
}

# Check that s_n is identical to the formula given in the instructions.
identical(s_n, (n*(n+1)*(2*n+1))/6)

