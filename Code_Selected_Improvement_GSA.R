#############
## DOCE  ####
#############
# Selective improvement of global datasets for the computation of locally relevant environmental indicators:
# a method based on global sensitivity analysis

# method					: GSA by Squared standardized regression coefficients
# Uncertainty propagation 	: Monte Carlo Simulation
# Last update				: May 24 2017
# Toolbox					: Statistics analysis

####################
## Import Library ##
####################

library("MASS") 			# For inverted matrix
library("ggplot2")			# Plotting
library("dplyr")			# Data management
library ("sensitivity")		# Sensitivity analysis
library ("boot") 			# Bootstrap techniques
library ("fitdistrplus")  	# Fitting distribution 
library ("triangle")		# Triangular distribution for Emission Factors
library ("reshape2")		# Graphics with ggplot2

#1. WORKING DIRECTORY
#setwd("C:/Users/...")

#######################
# Import input files ##
#######################



#################
## FUNCTIONS ####
#################


######################
## MODEL START HERE ##
######################


n = 5000

# Monte Carlo simulation for input parameters
set.seed(20)
p1 = data.frame(matrix(rnorm(n, 100, 10), nrow=n, ncol=1)) # For normal distribution

set.seed(20)
p2 = data.frame(matrix(runif(n, 5, 30), nrow=n, ncol=1)) # For uniform distribution

set.seed(20)
p3 = data.frame(matrix(log10(rltriangle(n, 10^0.003, 10^0.03, 10^0.01)), nrow=n, ncol=1)) #0.01 direct volatilization for grass (fert, crop residues)


# Add here your own model
# eg

model <- function (a1, a2, a3){
	out = a1 + a2 * a3 
	return (out)
	} 
####################
# GSA: SRC method ##
####################

# Make a matrix of input parameters n-column = parameters n-rows = Monte Carlo runs
# If some parameters are equals, please add a parameter once
# Parameters with 0 values are excluded for GSA

h1 = data.frame (p1, p2, p3)  # h1 is a data frame with input parameters


y = with (h1, (model (a1 = p1, a2 = p2, a3 = p3))) 	# Use sensitivity package to estimate the regression coefficient

x = src (X=h1, y=y,  nboot = 1000)     				# Estimate of standardized regression coefficient through bootstrap technique

print (x)											# Visualization of SRC

plot(x)

#
# Calculate the squared SRC
sensitivity_table = data.frame(as.matrix (x$SRC))
sensitivity_table$variable = row.names(sensitivity_table)
row.names(sensitivity_table) = NULL
new_data = data.frame(sensitivity_table$variable, sensitivity_table$original^2)  #Squared Regression coefficienct
names(new_data) = c("fact","coer_nue")
new_data = new_data[order(new_data[,2], decreasing = TRUE),]

windows (width=4, height=4)
barplot(new_data[c(1:20),2], xlab= "", ylab = "SRC", ylim = ,names = new_data[c(1:20),1],
		main= "Sensitivity indices")
sum(new_data[,2])

#