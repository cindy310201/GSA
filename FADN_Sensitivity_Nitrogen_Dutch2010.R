#################################################
## SENSITIVITY ANALYSIS  DUTCH DAIRY SYSTEM  ####
#################################################

#http://www.stat.umn.edu/geyer/old/5101/rlook.html
#http://tagteam.harvard.edu/hub_feeds/1981/feed_items/202072

##############
## NITROGEN ##
##############

#####################
##    Library   #####
#####################
library("Matrix")                       # Matrix calculations
library("MASS")                         # Inverted matrix
library("ggplot2")                      # Plot and graphs
library("dplyr")
library("grid")                         # 
library ("sensitivity")                 # Sensitivity analysis package
library ("boot")                        # bootstrap
library ("fitdistrplus")                # Fitting distribution 
library ("triangle")					# Triangular distribution for EF

#1. WORKING DIRECTORY
setwd("C:/Users/uwizeye/Documents/LEAP_PhD_Working_Documents/LEAP_working_documents/NUE_paper_framework/Life_cycle_NUE/3.Third paper")

#######################
# Import input files ##
#######################
# Data from FADN
feed   = read.table("1- Data/2- Dutch data/2. Data file/Feed_production.csv", header = T, sep = ",", row.names = NULL)  #Feed production
animal = read.table("1- Data/2- Dutch data/2. Data file/Dairy_production.csv", header = T, sep = ",", row.names = NULL)     #Dairy production

animal = animal[-c(250:268),]


#################
## FUNCTIONS ####
#################

# Fit function for all types of distributions
fit_function = function (data){
	windows(width=4, height=4)
	plotdist(data, histo=TRUE, demp=TRUE)
	windows(width=4, height=4)
	descdist(data, boot=1000)

	#output = data.frame(0, 16, 4)

	fn_rr = fitdist (data, "norm")
	summary(fn_rr)
	fln_rr = fitdist(data, "lnorm")
	summary(fln_rr)
	fg_rr = fitdist(data, "gamma")
	summary(fg_rr)
	fu_rr = fitdist(data, "unif")
	summary(fu_rr)

	
	windows(width=4, height=4)
	par(mfrow=c(2, 2))
	plot.legend =c("norm", "lognormal")
	denscomp(list(fn_rr, fln_rr), legendtext = plot.legend)
	qqcomp(list(fn_rr, fln_rr), legendtext = plot.legend)
	cdfcomp(list(fn_rr, fln_rr), legendtext = plot.legend)
	ppcomp(list(fn_rr, fln_rr), legendtext = plot.legend)

	windows(width=4, height=4)
	par(mfrow=c(2, 2))
	plot.legend =c("gamma", "unif")
	denscomp(list(fg_rr, fu_rr), legendtext = plot.legend)
	qqcomp(list(fg_rr, fu_rr), legendtext = plot.legend)
	cdfcomp(list(fg_rr, fu_rr), legendtext = plot.legend)
	ppcomp(list(fg_rr, fu_rr), legendtext = plot.legend)

}

# Function to truncate for negative values #Data: data frame, only for normal distribution
my_trunc = function (data){									
	x_output = data.frame (matrix(0,n,ncol(data)))
	for (i in 1:ncol(data)){
		if (data[1,i]==0 & data[2,i]==0){
			next 
		}
		set.seed(20)
		share = rnorm(n, data[1,i], data[2,i])
		share = sort(share)
		posi_share = which(share > 0)
		posi_share = min(posi_share)
		final_share = share [ posi_share:(length(share) - posi_share)]
		set.seed(20)
		final_share= sample(final_share, n, replace = TRUE)
		x_output[, i] =final_share
	}
	return(x_output)
}

# Function to truncate (high values)
my_trunc_high = function (data){									
	x_output = data.frame (matrix(0,n,ncol(data)))					# Only used for coefficient with values between 0 and 1
		share = data[,1]
		share = sort(share)
		posi_share = which(share < 1)
		posi_share = max(posi_share)
		final_share = share [ posi_share:(length(share) - posi_share)]
		final_share= sample(final_share, n, replace = TRUE)
		x_output[, 1] =final_share
		return(x_output)
}


######################
## MODEL START HERE ##
######################

#################
## Sample size ##
#################

n= 5000

#########P39

# ANIMAL#
#########
#1.V1       New born dairy calves                   number
#2.V2       Deceased dairy calves                   number
#3.V260     dairy cows                              number
#4.V261     Young stock (0-1 yr)                    number
#5.V262     Bulls (1-2 yr)                          number
#6.V263 	Young stock (1-2 yr)                    number
#7.V264 	Bulls (2yr+)                            number

farm_animal = data.frame (animal$V1, animal$V2, animal$V260, animal$V261, animal$V262, animal$V263, animal$V264)

animal_cate = colSums(farm_animal)

animal_cate = as.vector(as.matrix(animal_cate))

##
x_animal = data.frame(matrix(animal_cate, nrow=n, ncol=length(animal_cate), byrow=TRUE))

# Mortality 
af_dead = data.frame(matrix(988, n, 1))
ca_dead = data.frame(matrix(2829, n, 1))
mf_dead = data.frame(matrix(895, n, 1))
rf_dead = data.frame(matrix(109, n, 1))

 

# Replacement rate females
repl_rate = animal$V30
fit_function (repl_rate)
fln_rr = fitdist(repl_rate, "norm")
repl_rate_s = data.frame(matrix(c(fln_rr$estimate["mean"], fln_rr$estimate["sd"]), nrow=2, ncol=1))
x_repl_rate = my_trunc(data=repl_rate_s)
names(x_repl_rate)= c("P5")

af_exit = x_animal[,3] * x_repl_rate/100

#c. ANIMAL PARAMETERS
##########################
# DUTCH DATA from GLEAM ##
##########################

#Age at first calving in years
mean_afc = 2
sd_afc = 2*0.1
set.seed(20)
x_afc = data.frame(rnorm(n, mean_afc, sd_afc))
names(x_afc) = "P9"

#Live-weight adult female (kg)
mean_af_lw = 600
sd_af_lw = 600*0.1
set.seed(20)
x_af_lw = data.frame(rnorm(n, mean_af_lw, sd_af_lw))
names(x_af_lw) = "P10"

#Live_weight adult male (kg)
mean_am_lw = 722
sd_am_lw = 722*0.1
set.seed(20)
x_am_lw = data.frame (rnorm(n,mean_am_lw, sd_am_lw))
names(x_am_lw) = "P11"

#live_weight replacement females (kg)
mean_rf_lw = 450
sd_rf_lw = 450*0.1
set.seed(20)
x_rf_lw = data.frame(rnorm(n, mean_rf_lw, sd_rf_lw))
names(x_rf_lw) = "P12"

#live_weight replacement males (kg)
mean_rm_lw = 450  #Values in data too low, need to check
sd_rm_lw = 450*0.1
set.seed(20)
x_rm_lw = data.frame(rnorm(n, mean_rm_lw, sd_rm_lw))
names(x_rm_lw) = "P13"

#live_weight young females (1-2 years)
mean_mf_lw = 350
sd_mf_lw = 350*0.1
set.seed(20)
x_mf_lw = data.frame(rnorm(n, mean_mf_lw, sd_mf_lw))
names(x_mf_lw) = "P14"

#live_weight young males (1-2 years)
mean_mm_lw = 350
sd_mm_lw = 350*0.1
set.seed(20)
x_mm_lw = data.frame(rnorm(n, mean_mm_lw, sd_mm_lw))
names(x_mm_lw) = "P15"

#Calves weight (kg)
mean_ca_lw = 42
sd_ca_lw = 42*0.1
set.seed(20)
x_ca_lw = data.frame(rnorm(n, mean_ca_lw, sd_ca_lw))
names(x_ca_lw) = "P16"

#Growth rate
my_growth = function (w_adult,w_calf, afc){                 # w_adult: weight adult, w_calf: weight calf, afc: age at first calving
	growth = data.frame(matrix(0,n,1))
	growth[,1] = (w_adult - w_calf)/ (afc*365)
	return (growth)
}

female_growth = my_growth (w_adult= x_af_lw, w_calf= x_ca_lw, afc= x_afc)
male_growth   = my_growth (w_adult= x_am_lw, w_calf= x_ca_lw, afc= x_afc)


#Energy requirement
#1. Maintenance
my_ne_maint = function (weight, cf, const){  				# cf is a constant see IPCC guidelines Eq. 10.3, const= a correction factor
	maintenace = data.frame(matrix(0,n,1))
	maintenace[,1]= weight^0.75*cf*const
	return(maintenace)
}

af_ne_maint = my_ne_maint (weight= x_af_lw, cf=0.386, const=1)
rf_ne_maint = my_ne_maint (weight= x_rf_lw, cf=0.322, const=0.974)
mf_ne_maint = my_ne_maint (weight= x_mf_lw, cf=0.322, const=1)
am_ne_maint = my_ne_maint (weight= x_am_lw, cf=0.370, const=1)
rm_ne_maint = my_ne_maint (weight= x_rm_lw, cf=0.370, const=0.974)
mm_ne_maint = my_ne_maint (weight= x_mm_lw, cf=0.370, const=1)

#2. Activity
#Stall
#cact = 0       #cact refers to the percentage of manure deposited during grazing
#grazing
#cact=0.17*100/100
#ranging
#cact = 0.36*grazing/100

my_ne_activity = function (ne_maint,cact){
	activity = data.frame(matrix(0,n,1))
	activity [,1]= ne_maint*cact
	return(activity )
}

af_ne_activity = my_ne_activity (ne_maint= af_ne_maint, cact= 0.17)
rf_ne_activity = my_ne_activity (ne_maint= rf_ne_maint, cact= 0.17)
mf_ne_activity = my_ne_activity (ne_maint= mf_ne_maint, cact= 0.17)
am_ne_activity = my_ne_activity (ne_maint= am_ne_maint, cact= 0.17)
rm_ne_activity = my_ne_activity (ne_maint= rm_ne_maint, cact= 0.17)
mm_ne_activity = my_ne_activity (ne_maint= mm_ne_maint, cact= 0.17)

################################
## FERTILITY ##
ferti_mean = 0.8
set.seed(20)
x_ferti1 = data.frame(rnorm(n,ferti_mean, ferti_mean*0.1))
names(x_ferti1) = "P17"
x_ferti = my_trunc_high(x_ferti1)

###########################################
#3. Pregnancy
my_ne_pregnancy_af = function (ne_maint, fr){               # fr = fertility rate of females
	pregnancy = data.frame(matrix(0,n,1))
	pregnancy [,1]= ne_maint*0.10*fr/100
	return (pregnancy)
}
af_ne_pregnancy = my_ne_pregnancy_af (ne_maint= af_ne_maint , fr= x_ferti)

my_ne_pregnancy_rf = function (ne_maint, afc){
	pregnancy_rf = data.frame(matrix(0,n,1))
	pregnancy_rf [,1]= ne_maint*0.10/(0.5*afc)
	return (pregnancy_rf)
}

rf_ne_pregnancy = my_ne_pregnancy_rf (ne_maint= rf_ne_maint, afc= x_afc)

#4. Growth
my_ne_growth = function (weight, coe, mature_weight, weight_gain){
	ne_growth = data.frame(matrix(0,n,1))
	ne_growth[,1]= 22.02 * (weight/(coe*mature_weight))^0.75 * weight_gain^1.097    # Check for metabolic weight, not used in GLEAM file
	return (ne_growth)
}

rf_ne_growth = my_ne_growth (weight= x_rf_lw, coe= 0.8, mature_weight= x_af_lw, weight_gain = female_growth) # weight=replacement female, mature_weight=adult female
mf_ne_growth = my_ne_growth (weight= x_mf_lw, coe= 0.8, mature_weight= x_af_lw, weight_gain = female_growth) # weight=young female, mature_weight=adult female
rm_ne_growth = my_ne_growth (weight= x_rm_lw, coe= 1.02, mature_weight= x_am_lw, weight_gain = male_growth)  # weight=young female, mature_weight=adult female
mm_ne_growth = my_ne_growth (weight= x_mm_lw, coe= 1.02, mature_weight= x_am_lw, weight_gain = male_growth)  # weight=young female, mature_weight=adult female

## Milk production kg per year
milk_produ = animal$V17
## Fit PDF
fit_function (milk_produ)
# Selection of normal distribution
set.seed(20)
milk_year = data.frame (rnorm(n, mean(milk_produ), sd(milk_produ)))
names(milk_year) = "P18"


#Fat composition 
milkfat = animal$V31
fit_function (milkfat)
# Selection of normal distribution
fn_rm = fitdist(milkfat, "norm")
set.seed(20)
x_milkfat= data.frame(rnorm(n, fn_rm$estimate["mean"], fn_rm$estimate["sd"]))
names(x_milkfat)= c("P19")


# NE lactation
my_ne_lactation = function (milk,  milk_fat){                               # Milk density 1027–1033 kg/m3
	lactation = data.frame(matrix(0,n,1))
	lactation[,1]= (milk/365) * (1.47 + 0.40 * milk_fat)                    
	return(lactation)
}

af_ne_lactation = my_ne_lactation (milk= milk_year,  milk_fat= x_milkfat) #According to Niyireba,

# Total GE and Feed intake
af_ne_total = af_ne_maint + af_ne_activity + af_ne_pregnancy + af_ne_lactation
rf_ne_total = rf_ne_maint + rf_ne_activity + rf_ne_pregnancy
mf_ne_total = mf_ne_maint + mf_ne_activity
am_ne_total = am_ne_maint + am_ne_activity
rm_ne_total = rm_ne_maint + rm_ne_activity
mm_ne_total = mm_ne_maint + mm_ne_activity

# Digestibility
digest_mean = 0.7271
digest_sd = 0.0745
set.seed(20)
x_digest1 = data.frame(rnorm(n,digest_mean, digest_sd))
names(x_digest1) = "P20"
x_digest = my_trunc_high(x_digest1)*100


# Ratio of net energy available in diet for maintenance
my_rem = function(digest){
	rem = data.frame(matrix(0,n,1))
	rem[,1] = 1.123 - (4.092*10^-3*digest)+ (1.126*10^-5 * digest^2 - (25.4/digest))
	return (rem)
}
rem = my_rem (digest= x_digest)

#Ratio of net energy available for growth in a diet to digestible energy consumed
my_reg = function (digest){
	reg = data.frame (matrix(0,n,1))
	reg[,1] = 1.164 - (5.160*10^-3*digest) + ((1.308 *10^-5*digest^2)-(37.4/digest))
	return (reg)
}
reg = my_reg (digest= x_digest)

#Gross energy
# Make a MC for digestibility with a uniform distribution min=55 max=75 for pasture fed animals
af_ge = data.frame (matrix(0,n,1))
af_ge[,1] = (af_ne_total/rem)/(x_digest/100)

rf_ge = data.frame (matrix(0,n,1))
rf_ge[,1] = ((rf_ne_total/rem) + (rf_ne_growth/reg)) / (x_digest/100)

mf_ge = data.frame (matrix(0,n,1))
mf_ge[,1] = ((mf_ne_total/rem) + (mf_ne_growth/reg)) / (x_digest/100)

am_ge = data.frame (matrix(0,n,1))
am_ge[,1] = (am_ne_total/rem) / (x_digest/100)

rm_ge = data.frame (matrix(0,n,1))
rm_ge[,1] = ((rm_ne_total/rem) + (rm_ne_growth/reg)) / (x_digest/100)

mm_ge = data.frame (matrix(0,n,1))
mm_ge[,1] = ((mm_ne_total/rem) + (mm_ne_growth/reg)) / (x_digest/100)

# Energy content of ration
# White clover : GE Mean=18.3   sd=0.8  Min=17.5    Max=19.3
# Kikuyu        : Mean = 18.3    sd=

#Feed intake by day per animal
my_intake_day = function (ge){
	intake_day = data.frame (matrix(0,n,1))
	intake_day[,1] = ge / 18.45                     # Please estimate the energy content of the feed based on feed characteristics
	return(intake_day)                              # Replace 18.45 by a correct value (18.3)
}                                                   # Intake_day is the total intake in kg.DM/animal/day

af_intake_day = my_intake_day (ge = af_ge)
rf_intake_day = my_intake_day (ge = rf_ge)
mf_intake_day = my_intake_day (ge = mf_ge)
am_intake_day = my_intake_day (ge = am_ge)
rm_intake_day = my_intake_day (ge = rm_ge)
mm_intake_day = my_intake_day (ge = mm_ge)

#write.table(af_intake_day, "af_intake_day.csv", sep=",", row.names=FALSE)

#Annual intake
my_total_intake = function (number, day_intake){
	total_intake = data.frame (matrix(0,n,1))
	total_intake[,1] = number * day_intake * 365
	return (total_intake)
}

af_intake = my_total_intake (number= (x_animal[,3]+af_dead), day_intake= af_intake_day)
rf_intake = my_total_intake (number= (x_animal[,6]+rf_dead), day_intake= rf_intake_day)
mf_intake = my_total_intake (number= (x_animal[,4]+mf_dead), day_intake= mf_intake_day)
am_intake = my_total_intake (number= x_animal[,7], day_intake= am_intake_day)
rm_intake = my_total_intake (number= x_animal[,5], day_intake= rm_intake_day) 

######################
## Feed composition ##
######################

############
#   DAIRY ##
############
# Feed categories identified 
#---------------------------
# Concentrates uptake dairy cattle
# Wet by products uptake dairy cattle == DGGS from Barley
# Other feed uptake dairy cattle
# Maize product uptake dairy cattle
# Grass silage uptake dairy cattl
# Grass (fresh) uptake dairy cattle
# Milk products uptake dairy cattle

feed_uptake_energy = data.frame(feed$F56, feed$F57, feed$F59, feed$F60, feed$F62, feed$F63, feed$F58)

#share_feed = feed_uptake_cate / rowSums(feed_uptake_cate)
#Barley 1.27 + 9.51 + 0.68 +2.92 
#Soybean expeller 12.85 
#Maize 1.18 + 17.42
#Rapeseed expeller 0.32 + 4.52
#Palm kernel expeller 19.33 
#Molasses (Cane) 2.10 
#Beet pulp 6.34 
# Breaking down of concentrates section, based on Corina thesis page 45

concentrate_share = c(0.18332483,0.16381948,0.23712392,0.06170321,0.24643039,0.10759816)

concentrate = data.frame(matrix(0,nrow(feed_uptake_energy), ncol=6))
names(concentrate) = c("grain", "soybean", "corn", "rapeseed", "palm", "beet")
concentrate[,1] = feed_uptake_energy [,1] * concentrate_share[1]
concentrate[,2] = feed_uptake_energy [,1] * concentrate_share[2]
concentrate[,3] = feed_uptake_energy [,1] * concentrate_share[3]
concentrate[,4] = feed_uptake_energy [,1] * concentrate_share[4]
concentrate[,5] = feed_uptake_energy [,1] * concentrate_share[5]
concentrate[,6] = feed_uptake_energy [,1] * concentrate_share[6]



feed_category_new = data.frame(concentrate, feed_uptake_energy [,-1])


feed_share = feed_category_new / rowSums(feed_category_new)

# List of new 10 feed categories 
	# Grain (wheat etc) 
	# Soybean 
	# Corn
	# Rapeseed 
	# palm
	# beet
	# Wet by products uptake dairy cattle
	# Other feed uptake dairy cattle ==> Luzerne
	# Maize product uptake dairy cattle
	# Grass silage uptake dairy cattl
	# Grass (fresh) uptake dairy cattle
	# Milk products uptake dairy cattle

## Share of feed material
mean_fs = as.matrix(t(apply(feed_share,2,mean)))
sd_fs = as.matrix(t(apply(feed_share,2,sd)))

result_feed = as.data.frame(rbind(mean_fs, sd_fs))

x_share_feed = my_trunc (data=result_feed)
	names(x_share_feed) = c("P23","P24","P25","P26","P27","P28","P29","P30","P31","P32","P33","P34")

x_share_feed_new= x_share_feed/rowSums(x_share_feed)

#Use default protein content values 3.5% in the milk product: EUROLAC Schils company 18%
#x_n_milk_replacer = data.frame(rnorm(n, 180/6.25, 180/6.25*0.1)) #gN /kg milk replacer
set.seed(20)
x_n_milk_replacer = data.frame(rnorm(n, 180/6.25, 180/6.25*0.1)) #gN /kg milk replacer
	names(x_n_milk_replacer) = "P442"

milk_protein = animal$V4
fit_function (milk_protein)
fn_rp = fitdist(milk_protein, "norm")
set.seed(20)
x_proteinmilk= data.frame(rnorm(n, fn_rp$estimate["mean"], fn_rp$estimate["sd"]))
names(x_proteinmilk)= c("P36")


# Feed composition FADN data
feed_compo_fadn = data.frame(feed$F18, feed$F22, feed$F3, feed$F7, feed$F10)
feed_compo_fadn = feed_compo_fadn /1000 


#Result_co refers to mean and sd for compositions
result_co = data.frame(matrix(0, nrow=2, ncol(feed_compo_fadn)))
for (i in 1:ncol(feed_compo_fadn)){
	condition = is.na(feed_compo_fadn[,i]) == FALSE
	sss= fitdist(feed_compo_fadn[condition,i], "norm")
	result_co[,i] = summary(sss)[1]
}

x_feed_co = my_trunc (data=result_co)
names(x_feed_co) = c("P37","P38","P39","P40","P411")



# Feed composition concentrates from GLEAM (Paper 1) and   Feedpedia
	#     "grain", "soybean", "corn", "rapeseed", ::::     "palm", "molasse", "beet"

mean_ncont_conce = c(21.4,79,16,63,26.72,46.2)
sd_ncont_conce = mean_ncont_conce * 0.1

ncont_conce_t = data.frame (matrix(0,2,6))
	ncont_conce_t[1,] = mean_ncont_conce
	ncont_conce_t[2,] = sd_ncont_conce

ncont_conce = my_trunc(data=ncont_conce_t )
names(ncont_conce) = c("P41","P42","P43","P44","P45","P46")

# Feed composition N
feed_compo_n = data.frame(ncont_conce/1000, x_feed_co, x_n_milk_replacer/1000)

n_cont_feed_dairy = data.frame(matrix(rowSums(x_share_feed_new * feed_compo_n), nrow=n, ncol=1))


#1. Total N intake kg.P.cohort-1.day-1

my_n_intake = function (cohort, cont_feed){
	n_intake = data.frame (matrix(0, n, 1))
	n_intake [,1] = cohort *cont_feed 
	return (n_intake)
}

af_n_intake = my_n_intake (cohort = af_intake, cont_feed = n_cont_feed_dairy)
rf_n_intake = my_n_intake (cohort = rf_intake, cont_feed = n_cont_feed_dairy)
mf_n_intake = my_n_intake (cohort = mf_intake, cont_feed = n_cont_feed_dairy)
am_n_intake = my_n_intake (cohort = am_intake, cont_feed = n_cont_feed_dairy)
rm_n_intake = my_n_intake (cohort = rm_intake, cont_feed = n_cont_feed_dairy)

# Live animals
#3. Calculation N retention kgN.cohort-1.day-1 Following GLEAM Animal track algorithms!!

#calves
my_cn_retention = function (ckg, ne_replac, growf){
	cn_retention = data.frame (matrix(0, n, 1))
	cn_retention[,1]= ((ckg/365) * (368 -(7.03*ne_replac/growf))/1000)/6.25
	return (cn_retention)
}

c_n_ret = my_cn_retention (ckg= x_ca_lw, ne_replac= rf_ne_growth, growf= female_growth)

#Adult female
# Protein content of milk is delivered from Niyireba et al 2011 and GLEAM


my_afn_retention = function (milkday, protein_milk, number){
	afn_retention = data.frame (matrix(0, n, 1))
	afn_retention[,1]= ((milkday * protein_milk/100)/6.25 + c_n_ret) * number
	return (afn_retention)
}

af_n_ret = my_afn_retention (milkday = milk_year/365, protein_milk= x_proteinmilk, number=(x_animal[,3]+af_dead))

#Adult male, retention is 0
am_n_ret = data.frame (matrix(0, n, 1))
am_n_ret[,1] = 0

#Replacement females
my_rfn_retention = function (growf, ne_replac,  afc, number){
	rfn_retention = data.frame (matrix(0, n, 1))
	rfn_retention[,1]= (((growf * (268 - (7.03* ne_replac/growf)))/1000/6.25) + 1/afc *c_n_ret) * number 
	return (rfn_retention)
}

rf_n_ret = my_rfn_retention (growf = female_growth, ne_replac = rf_ne_growth, afc= x_afc, number= (x_animal[,6]+rf_dead))

# Replacement males
my_other_n_retention = function (grow, ne_replac, number){
	other_n_retention = data.frame (matrix(0, n, 1))
	other_n_retention[,1]= ((grow * (268-(7.03*ne_replac/grow)))/1000/6.25) * number 
	return (other_n_retention)
}

rm_n_ret = my_other_n_retention (grow= male_growth, ne_replac = rm_ne_growth, number= x_animal[,5])
mf_n_ret = my_other_n_retention (grow= female_growth, ne_replac = rf_ne_growth, number= (x_animal[,4]+mf_dead))

# 7. Development of formulas of N losses during manure management
# loss dead animals
my_loss_dead = function (grow, ne_replac, dead){
	loss_dead = data.frame (matrix(0, n, 1))
	loss_dead [,1] = ((grow * (268-(7.03*ne_replac/grow)))/1000/6.25) * dead
}

af_n_dead = my_afn_retention (milkday = milk_year/365, protein_milk= x_proteinmilk, number= af_dead)
am_n_dead = my_loss_dead (grow= male_growth, ne_replac = rm_ne_growth, dead= 0)
rf_n_dead = my_loss_dead (grow = female_growth, ne_replac = rf_ne_growth, dead= rf_dead)
rm_n_dead = my_loss_dead (grow= male_growth, ne_replac = rm_ne_growth, dead= 0)
mf_n_dead = my_loss_dead (grow= female_growth, ne_replac = mf_ne_growth, dead= mf_dead)
ca_n_dead = c_n_ret * ca_dead


#EXCRETION #######
my_excre = function (intake, retent, dead){
	excre_an = data.frame (matrix(0, n, 1))
	excre_an[,1] = intake - retent*365 - dead*365
	return (excre_an)
}
# Manure excreted
af_n_ex = my_excre (intake=af_n_intake, retent=af_n_ret, dead=(af_dead + ca_dead))
am_n_ex = my_excre (intake=am_n_intake, retent=am_n_ret, dead=0)
rf_n_ex = my_excre (intake=rf_n_intake, retent=rf_n_ret, dead=rf_dead)
rm_n_ex = my_excre (intake=rm_n_intake, retent=rm_n_ret, dead=0)
mf_n_ex = my_excre (intake=mf_n_intake, retent=mf_n_ret, dead=mf_dead)

#############################
## MANURE MANAGEMENT SYSTEM #
#############################
# Grazing time
# 1. calves(mm, mf)
# 2. Adult stock 
# 3. Young stock (rf, rm)
#% Grazing time: calves (young stock)
#% Grazing time july-august: dairy cows
#% Grazing time september and later: dairy cows
#% Grazing time: before july: dairy cows
#% Grazing time: youngstock 1-2 yr = replacement
grazing = data.frame(animal$V21, animal$V22,animal$V23,animal$V24,animal$V25)
grazing_time = data.frame(matrix(0, nrow=nrow(grazing), ncol=3))
	grazing_time[,1] = grazing[,1]                                              #Calves
	grazing_time[,2] = (grazing[,2] * 2 + grazing[,3]*4 + grazing[,4]*6 )/12    #Dairy cows
	grazing_time[,3] = grazing[,5]                                              #Young stock

result_grazing = data.frame(matrix(0, nrow=2, ncol(grazing_time)))
for (i in 1:ncol(grazing_time)){
	ssss= fitdist(grazing_time[,i], "norm")
	result_grazing[,i] = summary(ssss)[1]
}

x_grazing_time = my_trunc(data=result_grazing)
names(x_grazing_time) = c("P48", "P49", "P50")

stall_time = 1- x_grazing_time                          # Values for manure management are from GLEAM 
#share_manure_solid = stall_time *0.03                   # Solid
#share_manure_liquid = stall_time *0.97                  # Liquid

# Manure loss on Non agricultural land
# Emission Factors: DIRECT : Use Dutch country specific value from Feedprint
di_vol_min = c(0.0027, 0.001)    # Solid # Liquid
di_vol_max = c(0.010, 0.003)
di_vol_mean = c(0.005, 0.002)

x_di_vol = data.frame(matrix(0,n,length(di_vol_min)))

for (i in 1:length(di_vol_max)){
	set.seed(20)
	di_vol = log10(rltriangle(n, 10^di_vol_min[i], 10^di_vol_max[i], 10^di_vol_mean[i]))
	x_di_vol[,i] = di_vol
	names(x_di_vol) = c("P51", "P52")
}

# Emission factors : INDIRECT #DEFAUT DUTCH VALUES
# AF barn, AF grazing, MF, MM

# Solid manure
indi_vol_solid_mean  = c(0.105, 0.332, 0.117, 0.117)    # Solid manure
indi_vol_solid_sd = indi_vol_solid_mean *0.1
indi_vol_solid_t = data.frame(matrix(0,2,4))
	indi_vol_solid_t[1,] = indi_vol_solid_mean
	indi_vol_solid_t[2,] = indi_vol_solid_sd

x_indi_solid = my_trunc (data=indi_vol_solid_t)
names(x_indi_solid) = c("P53", "P54", "P200", "P201")



## Liquid manure
indi_vol_liquid_mean = c(0.102, 0.124, 0.112, 0.117)     # Liquid manure
indi_vol_liquid_sd = indi_vol_liquid_mean *0.1
indi_vol_liquid_t = data.frame(matrix(0,2,4))
	indi_vol_liquid_t[1,] = indi_vol_liquid_mean
	indi_vol_liquid_t[2,] = indi_vol_liquid_sd

x_indi_liquid = my_trunc (data=indi_vol_liquid_t)
names(x_indi_liquid) = c("P202", "P203", "P204", "P205")
                                                                  

## Leaching of manure 
leaching_ma = data.frame(matrix(c(0.020,0.002, 0.040,0.004), nrow=2, ncol=2))
x_leaching_ma = my_trunc (data=leaching_ma)
names(x_leaching_ma) = c("P55", "P56")

####
my_manure_loss = function (excre_cate, share_gr,share_man_gr,efgr, share_sta,share_man_sta,efsta){
	manure_loss = data.frame (matrix(0, n, 1))
	manure_loss[,1] = excre_cate *((share_gr*share_man_gr*efgr) +(share_sta*share_man_sta*efsta))
	return (manure_loss)
}

# Manure loss liquid : direct volatilization
af_manure_direct_liquid  = my_manure_loss (excre_cate = af_n_ex,share_gr= x_grazing_time[,2],share_man_gr=0.97, efgr=x_di_vol[,2], share_sta=stall_time[,2],share_man_sta=0.97,efsta=x_di_vol[,2])
am_manure_direct_liquid  = my_manure_loss (excre_cate = am_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.97, efgr=x_di_vol[,2], share_sta=stall_time[,1],share_man_sta=0.97,efsta=x_di_vol[,2])
rf_manure_direct_liquid  = my_manure_loss (excre_cate = rf_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.97, efgr=x_di_vol[,2], share_sta=stall_time[,3],share_man_sta=0.97,efsta=x_di_vol[,2])
rm_manure_direct_liquid  = my_manure_loss (excre_cate = rm_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.97, efgr=x_di_vol[,2], share_sta=stall_time[,3],share_man_sta=0.97,efsta=x_di_vol[,2])
mf_manure_direct_liquid  = my_manure_loss (excre_cate = mf_n_ex,share_gr= x_grazing_time[,1],share_man_gr=0.97, efgr=x_di_vol[,2], share_sta=stall_time[,1],share_man_sta=0.97,efsta=x_di_vol[,2])

# Manure loss liquid : indirect volatilization
af_manure_indirect_liquid  = my_manure_loss (excre_cate = af_n_ex,share_gr= x_grazing_time[,2],share_man_gr=0.97, efgr=x_indi_liquid[,2], share_sta=stall_time[,2],share_man_sta=0.97,efsta=x_indi_liquid[,2])
am_manure_indirect_liquid  = my_manure_loss (excre_cate = am_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.97, efgr=x_indi_liquid[,2], share_sta=stall_time[,1],share_man_sta=0.97,efsta=x_indi_liquid[,2])
rf_manure_indirect_liquid  = my_manure_loss (excre_cate = rf_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.97, efgr=x_indi_liquid[,3], share_sta=stall_time[,3],share_man_sta=0.97,efsta=x_indi_liquid[,3])
rm_manure_indirect_liquid  = my_manure_loss (excre_cate = rm_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.97, efgr=x_indi_liquid[,4], share_sta=stall_time[,3],share_man_sta=0.97,efsta=x_indi_liquid[,4])
mf_manure_indirect_liquid  = my_manure_loss (excre_cate = mf_n_ex,share_gr= x_grazing_time[,1],share_man_gr=0.97, efgr=x_indi_liquid[,3], share_sta=stall_time[,1],share_man_sta=0.97,efsta=x_indi_liquid[,3])

# Manure loss solid : direct volatilization
af_manure_direct_solid  = my_manure_loss (excre_cate = af_n_ex,share_gr= x_grazing_time[,2],share_man_gr=0.03, efgr=x_di_vol[,1], share_sta=stall_time[,2],share_man_sta=0.03,efsta=x_di_vol[,1])
am_manure_direct_solid  = my_manure_loss (excre_cate = am_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.03, efgr=x_di_vol[,1], share_sta=stall_time[,1],share_man_sta=0.03,efsta=x_di_vol[,1])
rf_manure_direct_solid  = my_manure_loss (excre_cate = rf_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.03, efgr=x_di_vol[,1], share_sta=stall_time[,3],share_man_sta=0.03,efsta=x_di_vol[,1])
rm_manure_direct_solid  = my_manure_loss (excre_cate = rm_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.03, efgr=x_di_vol[,1], share_sta=stall_time[,3],share_man_sta=0.03,efsta=x_di_vol[,1])
mf_manure_direct_solid  = my_manure_loss (excre_cate = mf_n_ex,share_gr= x_grazing_time[,1],share_man_gr=0.03, efgr=x_di_vol[,1], share_sta=stall_time[,1],share_man_sta=0.03,efsta=x_di_vol[,1])

# Manure loss solid : indirect volatilization
af_manure_indirect_solid  = my_manure_loss (excre_cate = af_n_ex,share_gr= x_grazing_time[,2],share_man_gr=0.03, efgr=x_indi_solid[,2], share_sta=stall_time[,2],share_man_sta=0.03,efsta=x_indi_solid[,2])
am_manure_indirect_solid  = my_manure_loss (excre_cate = am_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.03, efgr=x_indi_solid[,2], share_sta=stall_time[,1],share_man_sta=0.03,efsta=x_indi_solid[,2])
rf_manure_indirect_solid  = my_manure_loss (excre_cate = rf_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.03, efgr=x_indi_solid[,3], share_sta=stall_time[,3],share_man_sta=0.03,efsta=x_indi_solid[,3])
rm_manure_indirect_solid  = my_manure_loss (excre_cate = rm_n_ex,share_gr= x_grazing_time[,3],share_man_gr=0.03, efgr=x_indi_solid[,4], share_sta=stall_time[,3],share_man_sta=0.03,efsta=x_indi_solid[,4])
mf_manure_indirect_solid  = my_manure_loss (excre_cate = mf_n_ex,share_gr= x_grazing_time[,1],share_man_gr=0.03, efgr=x_indi_solid[,3], share_sta=stall_time[,1],share_man_sta=0.03,efsta=x_indi_solid[,3])

######
# Leaching
my_manure_leach = function (excre_cate, efleach, share){
	manure_leach = data.frame (matrix(0, n, 1))
	manure_leach[,1] = excre_cate *efleach *share
	return (manure_leach)
}

# Liquid
af_manure_leaching_liquid  = my_manure_leach (excre_cate = af_n_ex, share=0.97, efleach=x_leaching_ma[,2])
am_manure_leaching_liquid  = my_manure_leach (excre_cate = am_n_ex, share=0.97, efleach=x_leaching_ma[,2])
rf_manure_leaching_liquid  = my_manure_leach (excre_cate = rf_n_ex, share=0.97, efleach=x_leaching_ma[,2])
rm_manure_leaching_liquid  = my_manure_leach (excre_cate = rm_n_ex, share=0.97, efleach=x_leaching_ma[,2])
mf_manure_leaching_liquid  = my_manure_leach (excre_cate = mf_n_ex, share=0.97, efleach=x_leaching_ma[,2])

#Solid
af_manure_leaching_solid  = my_manure_leach (excre_cate = af_n_ex, share=0.03, efleach=x_leaching_ma[,1])
am_manure_leaching_solid  = my_manure_leach (excre_cate = am_n_ex, share=0.03, efleach=x_leaching_ma[,1])
rf_manure_leaching_solid  = my_manure_leach (excre_cate = rf_n_ex, share=0.03, efleach=x_leaching_ma[,1])
rm_manure_leaching_solid  = my_manure_leach (excre_cate = rm_n_ex, share=0.03, efleach=x_leaching_ma[,1])
mf_manure_leaching_solid  = my_manure_leach (excre_cate = mf_n_ex, share=0.03, efleach=x_leaching_ma[,1])

########
## Total 
af_manure_loss  = af_manure_direct_liquid + af_manure_indirect_liquid + af_manure_direct_solid + af_manure_indirect_solid + af_manure_leaching_liquid + af_manure_leaching_solid
am_manure_loss  = am_manure_direct_liquid + am_manure_indirect_liquid + am_manure_direct_solid + am_manure_indirect_solid + am_manure_leaching_liquid + am_manure_leaching_solid
rf_manure_loss  = rf_manure_direct_liquid + rf_manure_indirect_liquid + rf_manure_direct_solid + rf_manure_indirect_solid + rf_manure_leaching_liquid + rf_manure_leaching_solid
rm_manure_loss	= rm_manure_direct_liquid + rm_manure_indirect_liquid + rm_manure_direct_solid + rm_manure_indirect_solid + rm_manure_leaching_liquid + rm_manure_leaching_solid
mf_manure_loss  = mf_manure_direct_liquid + mf_manure_indirect_liquid + mf_manure_direct_solid + mf_manure_indirect_solid + mf_manure_leaching_liquid + mf_manure_leaching_solid


###########
## Total ##
########### 

my_total = function (af,am,rf,rm,mf,mm, ca, days){
	total_an = data.frame (matrix(0, n, 1))
	total_an[,1] = (af + am + rf + rm + mf + mm + ca) * days/1000
	return (total_an)
}


# Total dead animals 
herd_n_dead = my_total (af=af_n_dead, am=am_n_dead, rf=rf_n_dead, rm=rm_n_dead, 
				   mf=mf_n_dead, mm=0, ca= ca_n_dead, days=365) 

# 4. Total N intake, N retention, N excretion per Herd.year-1, already in per year
herd_n_intake =  my_total (af=af_n_intake, am=am_n_intake, rf=rf_n_intake, rm=rm_n_intake, 
							   mf=mf_n_intake, mm=0, ca=0, days=1) 

herd_n_ret = my_total (af=af_n_ret, am=am_n_ret, rf=rf_n_ret, rm=rm_n_ret, 
				   mf=mf_n_ret, mm=0, ca=0, days=365)

herd_n_ex = my_total(af=af_n_ex, am=am_n_ex, rf=rf_n_ex, rm=rm_n_ex, 
							   mf=mf_n_ex, mm=0, ca=0, days=1)

herd_n_loss = my_total(af=af_manure_loss, am=am_manure_loss, rf=rf_manure_loss, rm=rm_manure_loss, 
							   mf=mf_manure_loss, mm=0, ca=0, days=1)

animal_input = herd_n_intake

#####
#5. Fraction return N
frac_n_return = data.frame (matrix(0, n, 1))
frac_n_return[,1] = herd_n_ex / herd_n_intake

#Total annual parameters MILK and TISSUES
milk_ret_year = (af_n_ret - c_n_ret) * 365/1000

#Edible and non-edible tissues
tissu_ret  = data.frame (matrix(0, n, 1))
tissu_ret[,1] = herd_n_ret - milk_ret_year

# Nutrient losses at animal level
animal_loss = herd_n_loss + herd_n_dead

# Recycled in crop production
manure_recycle = herd_n_ex - herd_n_loss

## Output animals
animal_output = herd_n_intake - animal_loss

##########################################
####  III. Animal products processing stage
#milk_processing

milk_input = milk_ret_year

my_proc_output = function (prod, loss){
	proc_out = data.frame (matrix(0, n, 1))
	proc_out[,1] = prod  - prod * loss
	return (proc_out)
}

milk_output = my_proc_output (prod = milk_input, loss = 0.0015)
#output meat See FAO document: http://www.fao.org/wairdocs/lead/x6114e/x6114e00.htm #Contents
#0.0015 amount of P loss per kg of P in milk

#Meat processing : Culling cows, meat female and meat male
#For af, mf, mm, we consider 
af_n_ret_slaughter = (c_n_ret * af_exit) * 365/1000
mf_n_ret_slaughter = my_other_n_retention (grow= female_growth, ne_replac = rf_ne_growth, number= af_exit)*365/1000


slaughter_input = data.frame (matrix(0, n, 1))
slaughter_input[,1] = af_n_ret_slaughter + mf_n_ret_slaughter      #input meat

# See if we can add also mm and mf animals
# Dressing carcass
carcassd = data.frame(matrix(c(0.55, 0.05), nrow=2,ncol=1))              #Value in GLEAM
x_carcassd = my_trunc(data=carcassd)
names(x_carcassd)= "P57"
#max_carcassd = 0.60                 #Value from GLEAM


# output meat See FAO document: http://www.fao.org/wairdocs/lead/x6114e/x6114e00.htm #Contents
#                               http://www.fao.org/wairdocs/lead/x6114e/x6114e04.htm
# Carcass and edible products   62-64
# Edible Offal                  3-4%
# Blood                         3-4%
# Inedible raw material         8-10%
# Hide                          7%
# manure                        8%
# Shrinkage                     2-10%

use_inedible = data.frame(matrix(c(0.16, 0.016), nrow=2,ncol=1)) # this range refers to the % of non-edible products used in other industries
x_use_inedible = my_trunc(data=use_inedible)
names(x_use_inedible) = "P58"

organic_waste = data.frame(matrix(c(0.07, 0.007), nrow=2,ncol=1))       # Organic waste refers to manure from slaughterhouse use for fertilizer
x_organic_waste = my_trunc(data=organic_waste)
names(x_organic_waste) = "P59"

my_sproc_output = function (prod, dress_car, use_ined, organic_waste){
	sproc_out = data.frame (matrix(0, n, 1))
	sproc_out[,1] = prod * (dress_car + use_ined + organic_waste )
	return (sproc_out)
}

slaughter_output = my_sproc_output (prod = slaughter_input, dress_car= x_carcassd, 
									use_ined= x_use_inedible, organic_waste= x_organic_waste)

#Final products
process_input = data.frame (matrix(0, n, 1))
process_input[,1] = milk_input + slaughter_input

process_output = data.frame (matrix(0, n, 1))
process_output[,1] = milk_output + slaughter_output

process_loss = data.frame (matrix(0, n, 1))
process_loss[,1] = process_input  - process_output


############################################
## CROP PRODUCTION : Grass + White clover ##
############################################

#total_area = sum(animal[,32] - animal[,34])
# At farm level
# Feed categories identified 
	#1"grain"
	#2"soybean"
	#3"corn"
	#"4rapeseed"
	#"5palm"
	#6 "molasse" #"beet
	#7 Wet by products uptake dairy cattle
	#8 Other feed uptake dairy cattle (Luzerne)
	#9 Maize product uptake dairy cattle
	#10 Grass silage uptake dairy cattl
	#11 Grass (fresh) uptake dairy cattle
	#12 Milk products uptake dairy cattle

# Fertiliser use
fertiliser_fadn = data.frame("maize"=animal$V159, "grass"=animal$V164)
result_fe = data.frame(matrix(0, nrow=2, ncol(fertiliser_fadn)))
for (i in 1:ncol(fertiliser_fadn)){
	condition = fertiliser_fadn[,i]>0
	sss_fe= fitdist(fertiliser_fadn[condition,i], "norm")
	result_fe[,i] = summary(sss_fe)[1]
}

fertiliser_n = data.frame(matrix(0,nrow=2,ncol=11)) # Value from GLEAM SD and van Middelaar - mean*0.1
#GRAINS CORN    MLSOY   MLRAPE
#97.2   4.9 51.6    180 4.9 108
	grain_fe = c(97.2, 97.2*0.1)
	soy_fe   = c(51.6, 51.6*0.1)
	corn_fe  = c(150, 150*0.1)
	rape_fe  = c(180, 180*0.1)
	palm_fe  = c(104, 104*0.1)
	beet_fe  = c(108, 108*0.1)
	ddgs_fe  = c(97.2, 97.2*0.1)

	fertiliser_n[,1] = grain_fe
	fertiliser_n[,2] = soy_fe
	fertiliser_n[,3] = corn_fe
	fertiliser_n[,4] = rape_fe
	fertiliser_n[,5] = palm_fe
	fertiliser_n[,6] = beet_fe
	fertiliser_n[,7] = ddgs_fe
	fertiliser_n[,8] = result_fe[,1]
	fertiliser_n[,9] = result_fe[,1]
	fertiliser_n[,10] = result_fe[,2]
	fertiliser_n[,11] = result_fe[,2]

x_fertiliser = my_trunc(data=fertiliser_n)
names(x_fertiliser) = c("P60","P61","P62","P63","P64","P65","P66","P67","P68","P69","P70")

##############
# Manure use
manure_fadn = data.frame("maize"=animal$V162, "grass"=animal$V163)
result_ma = data.frame(matrix(0, nrow=2, ncol(manure_fadn)))
for (i in 1:ncol(manure_fadn)){
	condition = manure_fadn[,i]>0
	sss_ma= fitdist(manure_fadn[condition,i], "norm")
	result_ma[,i] = summary(sss_ma)[1]
}


manure_n = data.frame(matrix(0,nrow=2,ncol=11)) # Use GLEAM /Van Middelar C.data for imported feed
#GRAINS CORN    MLSOY   MLRAPE
#52.72  52.72   9   52.72
	grain_ma = c(52.72, 52.72*0.1)
	soy_ma   = c(9, 9*0.1)
	corn_ma  = c(50, 50*0.1)
	rape_ma  = c(28, 28*0.1)
	palm_ma  = c(0,0)
	beet_ma  = c(116, 116*0.1)
	ddgs_ma  = c(52.72, 52.72*0.1)

	manure_n[,1] = grain_ma
	manure_n[,2] = soy_ma
	manure_n[,3] = corn_ma
	manure_n[,4] = rape_ma
	manure_n[,5] = palm_ma
	manure_n[,6] = beet_ma
	manure_n[,7] = ddgs_ma
	manure_n[,8] = result_ma[,1]
	manure_n[,9] = result_ma[,1]
	manure_n[,10] = result_ma[,2]
	manure_n[,11] = result_ma[,2]

#manure_new = manure_n[,-5]
x_manure_co = my_trunc(data=manure_n)
names(x_manure_co) = c("P71","P72","P73","P74","P175","P76","P77","P78","P79","P80","P81")



# Organic fertliser 
organic_fert = animal$V99       # Less important but we can add it
condition = organic_fert>0
fn_ro = fitdist (organic_fert[condition], "norm")
summary(fn_ro)
orga_co = rnorm(n, fn_ro$estimate["mean"], fn_ro$estimate["sd"])
	orga_co = sort(orga_co)
	posi_orga_co = which(orga_co > 0)
	posi_orga_co = min(posi_orga_co)
	final_orga_co = orga_co [ posi_orga_co:(length(orga_co) - posi_orga_co)]
x_n_organic = data.frame(matrix(sample(final_orga_co, n, replace = TRUE), nrow=n, ncol=1))
names(x_n_organic) = "P82"
# Fill table
organic_fertiliser = data.frame(matrix(0,n,11))
organic_fertiliser[,8] = x_n_organic
organic_fertiliser[,9] = x_n_organic

#MANURE TOTAL
x_manure_f = x_manure_co + organic_fertiliser
###########

## N fixation 
fixation_fadn =  animal$V117
condition = fixation_fadn>0
sss_fi= fitdist(fixation_fadn[condition], "norm")
result_fix = summary(sss_fi)[1]


fixation_n = data.frame(matrix(0,nrow=2,ncol=10))
	grain_fix = c(5, 5*0.1)
	soy_fix   = c(115, 115*0.1)
	corn_fix  = c(5, 5*0.1)
	rape_fix  = c(110, 110*0.1)
	palm_fix  = c(0,0)
	beet_fix  = c(116, 116*0.1)
	ddgs_fix  = c(5, 5*0.1)
	alfa_fix  = c(160, 160*0.1) # See Dolman email
	fixation_n[,1] = grain_fix
	fixation_n[,2] = soy_fix
	fixation_n[,3] = corn_fix
	fixation_n[,4] = rape_fix
	fixation_n[,5] = palm_fix
	fixation_n[,6] = beet_fix
	fixation_n[,7] = ddgs_fix
	fixation_n[,8] = result_fix[1]
	fixation_n[,9] = alfa_fix
	fixation_n[,10] =result_fix[1]
	fixation_n[,11] = result_fix[1]

x_fixation_f = my_trunc(data=fixation_n)
names(x_fixation_f) = c("P83","P84","P85","P86","palm", "P87","P88","P89","P90","P91","P92")


#################################
# N deposition
deposition_fadn = animal$V118
condition = deposition_fadn>0
fn_d = fitdist(deposition_fadn[condition], "norm")
result_dep= summary(fn_d)[1]

deposition_n = data.frame(matrix(0,nrow=2,ncol=11))
	grain_dep = c(15.87, 15.87*0.1)
	soy_dep   = c(4.32, 4.32*0.1)
	corn_dep  = c(11.56, 11.56*0.1)
	rape_dep  = c(12.52, 12.52*0.1)
	palm_dep  = c(4.91,4.91*0.1)
	beet_dep = c(25.87, 25.87*0.1)
	ddgs_dep = c(15.87, 15.87*0.1)

	deposition_n[,1] = grain_dep
	deposition_n[,2] = soy_dep
	deposition_n[,3] = corn_dep
	deposition_n[,4] = rape_dep
	deposition_n[,5] = palm_dep
	deposition_n[,6] = beet_dep
	deposition_n[,7] = ddgs_dep
	deposition_n[,8] = result_dep[1]
	deposition_n[,9] = result_dep[1]
	deposition_n[,10] = result_dep[1]
	deposition_n[,11] = result_dep[1]

x_at_deposi = my_trunc(data=deposition_n)
names(x_at_deposi) = c("P93","P94","P95","P96","P97","P98","P99","P100","P101","P102", "P210")

################################
# Yield

yield_fadn = data.frame (animal$V177,animal$V180)
result_yie = data.frame(matrix(0, nrow=2, ncol(yield_fadn)))
for (i in 1:ncol(yield_fadn)){
	condition = yield_fadn[,i]>0
	sss_yie= fitdist(yield_fadn[condition,i], "norm")
	result_yie[,i] = summary(sss_yie)[1]
}

#3087.54    5799.1  2434.25 1440.26 5799.1  3113.6
	grain_yield = c(6783, 6783*0.1)
	soy_yield   = c(2434.25, 2434.25*0.1)
	corn_yield  = c(7652.59, 7652.59*0.1)
	rape_yield  = c(1440.26, 1440.26*0.1)
	palm_yield  = c(21210, 21210*0.1)
	beet_yield  = c(3113.6, 3113.6*0.1)
	ddgs_yield  = c(6783, 6783*0.1)

	yield = data.frame(matrix(0,nrow=2,ncol=11))
	yield[,1] = grain_yield
	yield[,2] = soy_yield
	yield[,3] = corn_yield
	yield[,4] = rape_yield
	yield[,5] = palm_yield
	yield[,6] = beet_yield
	yield[,7] = ddgs_yield
	yield[,8] = result_yie[,1]  # get new yield for other
	yield[,9] = result_yie[,1]
	yield[,10] = result_yie[,2]
	yield[,11] = result_yie[,2]

x_yield = my_trunc(data=yield)
names(x_yield) = c("P103","P104","P105","P106","P107","P108","P109","P110","P111","P112","P113")

x_crop_resid = data.frame (matrix(0, nrow=n, ncol=11))
	#"grain"
	#"soybean"
	#"corn"
	#"rapeseed"
	#"palm"
	#"molasse" #"beet
	# Wet by products uptake dairy cattle
	# Other feed uptake dairy cattle
	# Maize product uptake dairy cattle
	# Grass silage uptake dairy cattl
	# Grass (fresh) uptake dairy cattle
	# Milk products uptake dairy cattle
	x_crop_resid[,1] = ((520 +1.51*yield[,1])*(0.24*0.009+(1-0.9)*0.006))
	x_crop_resid[,2] = ((1350+0.93*yield[,2])*(0.19*0.008+0.008))
	x_crop_resid[,3] = ((610+1.03*yield[,3])*(0.19*0.007+0.006))
	x_crop_resid[,4] = ((1350+0.93*yield[,4])*(0.19*0.008+0.008))
	x_crop_resid[,5] = data.frame(matrix(40, nrow=n, ncol=1))
	x_crop_resid[,6] = data.frame(matrix(20, nrow=n, ncol=1))
	x_crop_resid[,7] = ((520 +1.51*yield[,7])*(0.24*0.009+(1-0.9)*0.006))
	x_crop_resid[,8] = (0.29*yield[,8] *0.027)
	x_crop_resid[,9] = ((610+1.03*yield[,9])*(0.19*0.007+0.006))
	x_crop_resid[,10] = 0.3*yield[,10]*0.015
	x_crop_resid[,11] = 0.3*yield[,11]*0.015

# Monte carlo simulation for crop residues, see Vellinga et al 2013 Feedpedia page 36

#x_crop_resid = data.frame(matrix(0,n,ncol(crop_resid)))
#for (i in 1:ncol(crop_resid)){
#	resid = rnorm (n, crop_resid[,i], crop_resid[,i]*0.1)
#	x_crop_resid[,i] = resid
#	names(x_crop_resid) = c("P114","P115","P116","P117","P118","P119","P120","P121","P122","P123","P125")
#}

# For the sensitivity analysis, use 1:11 remove(7)
##############################
### COMPUTATION START HERE  ##
##############################

# Emission factors, range from IPCC guidelines chapter 11
set.seed(20)
ef_dir_fr = data.frame(matrix(rep(log10(rltriangle(n, 10^0.003, 10^0.03, 10^0.01)),11), nrow=n, ncol=11)) #0.01 direct volatilization for grass (fert, crop residues)
#names (ef_dir_fr)= "P126"
set.seed(20)
ef_dir_m = data.frame(matrix(rep(log10(rltriangle(n, 10^0.007, 10^0.06, 10^0.02)),11), nrow=n, ncol=11))  #0.02 direct volatilization for grass (manure)
#names (ef_dir_m)= "P127"
set.seed(20)
ef_ind_m = data.frame(matrix(rep(log10(rltriangle(n, 10^0.05, 10^0.5, 10^0.2)),11), nrow=n, ncol=11))    #0.2 indirect volatilization for all crops (manure)
#names (ef_ind_m)=  "P128"
set.seed(20)
ef_ind_fr = data.frame(matrix(rep(log10(rltriangle(n, 10^0.03, 10^0.3, 10^0.1)),11), nrow=n, ncol=11))   # 0.1 indirect volatilization for all crops (fertilizer, crop residues)
#names (ef_ind_fr)=  "P129"


# Calculation inputs per ha
my_input_n_crop = function (fert, crop_resid, manur, at_dep, biof){
	input_n_crop = data.frame(matrix (0, n, 11))
	input_n_crop =  fert + crop_resid + manur + at_dep + biof
	return (input_n_crop)
}

total_input_ha = my_input_n_crop (fert=x_fertiliser, crop_resid= x_crop_resid, manur= x_manure_f, 
								at_dep= x_at_deposi, biof= x_fixation_f)

#Calculation of outputs per ha
my_output_n_crop = function (yield, ncont) {
	output_n_crop = data.frame(matrix (0, n, 11))
	output_n_crop = yield * ncont
	return (output_n_crop)
}

total_output_ha = my_output_n_crop (yield= x_yield, ncont= feed_compo_n[1:11])     #Crop yield data from GLEAM Model

#########   #"grain"
	#"soybean"
	#"corn"
	#"rapeseed"
	#"palm"
	#"molasse"
	#"beet
	# Wet by products uptake dairy cattle
	# Other feed uptake dairy cattle
	# Maize product uptake dairy cattle
	# Grass silage uptake dairy cattl
	# Grass (fresh) uptake dairy cattle
	# Milk products uptake dairy cattle
#0.03   0.03    0.05    0.03    0.03    0.03
runoff_mean = c(0.03,0.05,0.01,0.05,0.03,0.03,0.03,0.03,0.03,0.03,0.03)
runoff_sd = runoff_mean *0.1
runoff = data.frame(matrix(0,2,11))
	runoff[1,] = runoff_mean
	runoff[2,] = runoff_sd 

x_runoff = my_trunc(data=runoff)
names(x_runoff) = c("P130","P131","P132","P133","P134","P135","P136","P137","P138","P139","P140")

## LEACHING
leaching_mean = c(0.25,0.25,0.09,0.25,0.33,0.33,0.33,0.33,0.33,0.33,0.33)
leaching_sd = leaching_mean *0.1

leaching = data.frame(matrix(0,2,11))
	leaching[1,]= leaching_mean
	leaching[2,]= leaching_sd

x_leaching = my_trunc(data=leaching)
names(x_leaching) = c("P141","P142","P143","P144","P145","P146","P147","P148","P149","P150","P151")


# 1st phase of N losses in crop
my_loss_crop = function (fert, crop_resid, manur, runoff, effr, efm, efindfr, efindm){
	soil_loss = data.frame(matrix (0, n, 11))
	dir_loss = data.frame(matrix (0, n, 11))
	dir_loss[,10:11] = (fert[,10:11] + crop_resid[,10:11]) * effr[,10:11] + manur[,10:11] * efm[,10:11] # Grassf
	dir_loss[,1:9] = (fert[,1:9] + crop_resid[,1:9]+ manur[,1:9] ) * effr[,1:9]  # Crop
	ind_loss = (fert + crop_resid) * efindfr + manur * efindm
	run_off = (fert + manur) * runoff #Values estimated from Velthof et al., 2009 MITERRA Model
	soil_loss  = dir_loss + ind_loss + run_off 
	return (soil_loss)
}

loss_crop = my_loss_crop (fert = x_fertiliser, crop_resid= x_crop_resid, manur= x_manure_f, runoff= x_runoff, 
						  effr=ef_dir_fr, efm= ef_dir_m, efindfr=ef_ind_fr , efindm= ef_ind_m)

# Mineralization fraction for organic N
miner_grass = data.frame(matrix(rep(rnorm (n, 0.10, 0.01),11), nrow=n, ncol=11))
names(miner_grass) = c("P152")
miner_crop = data.frame(matrix(rep(rnorm(n, 0.3,0.03),11), nrow=n, ncol=11))
names(miner_crop) = c("P153")

# Organic N directly stored based on (Dollé and Smati, 2005; Velthof et al., 2009
my_soil_stock_organic = function( manur, runoff, crop_resid, efm, efindm, miner_g, miner_c){
  stock_manure = data.frame(matrix (0, n, 11))
  stock_manure[,10:11] = (manur[,10:11] - (manur[,10:11] * runoff[,10:11] + manur[,10:11]*(efm[,10:11]+
						  efindm[,10:11])))*miner_g[,10:11]
  stock_manure[,c(1:9)] = (manur[,c(1:9)] - (manur[,c(1:9)] * runoff[,c(1:9)] + 
							  manur[,c(1:9)]*(efm[,c(1:9)] + efindm[,c(1:9)])))*miner_c[,c(1:9)]
  stock_resid = data.frame(matrix (0, n, 11))
  stock_resid[,10:11] =  (crop_resid[,10:11] - (crop_resid[,10:11] * runoff[,10:11] + 
						 crop_resid[,10:11]*(efm[,10:11] + efindm[,10:11])))*miner_g[,10:11]
  stock_resid[,c(1:9)] = (crop_resid[,c(1:9)] - (crop_resid[,c(1:9)] * runoff[,c(1:9)] + 
						 crop_resid[,c(1:9)]*(efm[,2] + efindm[,2])))*miner_c[,c(1:9)]
  stock_organic  = stock_manure  + stock_resid
  return (stock_organic)
}

organic_stock  = my_soil_stock_organic (manur=x_manure_f, runoff=x_runoff, crop_resid= x_crop_resid, 
										efm= ef_dir_m, efindm= ef_ind_m, miner_g= miner_grass, miner_c= miner_crop)

# Soil N surplus
my_surplus = function (){
  surplus = total_input_ha - loss_crop - organic_stock - total_output_ha
  return (surplus)
}
soil_surplus = my_surplus()

# Leaching
my_leaching = function (surpl, leaching_rate){  #Leaching parameter does not contain the 1st colomn [,-1]
	soil_leaching = data.frame(matrix (0, nrow(surpl), ncol(surpl)))
	for (i in seq_len(nrow(surpl))){
		for (j in seq_len(ncol(surpl))){
			if (surpl [i,j] > 0) {
			soil_leac = surpl [i,j]*leaching_rate[i,j] +(surpl [i,j] *(1-leaching_rate[i,j]))*0.7  #70% of the surplus is lost via leaching Velthof
			if (surpl [i,j]<=0){
			soil_leac = 0
		}
		soil_leaching [i,j] = soil_leac
	  }
	}
  }
  return (soil_leaching)
}

leaching_soil = my_leaching ( surpl = soil_surplus, leaching_rate=x_leaching)

######
# Total N losses
loss_crop_ha = leaching_soil  + loss_crop


#stock N change
my_stock_change =  function (surpl, leaching_rate){
  posi_change = data.frame(matrix (0, nrow(surpl), ncol(surpl)))
  surpl = soil_surplus
  for (i in seq_len(nrow(surpl))){
	for (j in seq_len(ncol(surpl))){
	  if (surpl[i,j] > 0){
		posi_cha = (surpl [i,j] *(1-leaching_rate[i,j]))*0.3
	  }
	  if (surpl [i,j] <=0){
		posi_cha = 0
	  }
	  posi_change[i,j] = posi_cha
	}
  }
  neg_change = data.frame(matrix (0, nrow(surpl), ncol(surpl)))
  for (i in seq_len(nrow(surpl))){
	for (j in seq_len(ncol(surpl))){
	  if (surpl[i,j] <0){
		neg_cha = surpl[i,j]
	  }
	  if (surpl[i,j]>=0){
		neg_cha = 0
	  }
	  neg_change[i,j] = neg_cha
	}
  }
  stock_change= organic_stock  + posi_change + neg_change 
  return (stock_change)
}

stock_change_ha = my_stock_change (surpl = soil_surplus, leaching_rate=x_leaching)

#############
# Total herd intake by feed material
n_feed_mat_dairy = x_share_feed_new * feed_compo_n

my_n_intake_fe = function (cohort, cont_feed){
	n_intake_fe = data.frame (matrix(0, n, 12))
	cohort = as.vector(as.matrix(cohort))
	cohort = data.frame(matrix(cohort, nrow=n, ncol=12, byrow=FALSE))
	n_intake_fe  = cohort *cont_feed
	return (n_intake_fe)
}

herd_n_intake_total = data.frame(matrix(0, n, 12))
herd_n_intake_total = my_n_intake_fe (cohort = af_intake, cont_feed = n_feed_mat_dairy)+
					  my_n_intake_fe (cohort = rf_intake, cont_feed = n_feed_mat_dairy)+
					  my_n_intake_fe (cohort = mf_intake, cont_feed = n_feed_mat_dairy)+
					  my_n_intake_fe (cohort = am_intake, cont_feed = n_feed_mat_dairy)+
					  my_n_intake_fe (cohort = rm_intake, cont_feed = n_feed_mat_dairy)

feed_loss = feed$F97
fit_function (feed_loss)
fln_feed = fitdist(feed_loss, "lnorm")
set.seed(20)
x_feed_loss = data.frame(rlnorm (n,fln_feed$estimate["meanlog"], fln_feed$estimate["sdlog"] ))
names(x_feed_loss) = "P154"

herd_n_loss_total = data.frame(matrix(0, n, 12))
herd_n_loss_total =  my_n_intake_fe (cohort = x_feed_loss, cont_feed = x_share_feed_new)
					
herd_n_feed_use = herd_n_intake_total + herd_n_loss_total					

###############
	#"soybean"
	#"corn"
	#"rapeseed"
	#"palm"
	#"molasse"
	#"beet
	# Wet by products uptake dairy cattle
	# Other feed uptake dairy cattle
	# Maize product uptake dairy cattle
	# Grass silage uptake dairy cattl
	# Grass (fresh) uptake dairy cattle
	# Milk products uptake dairy cattle
# Feed use and production
fue_max = c(0.85, 0.54)
fue_min = c(0.75, 0.45)


x_fue = data.frame(matrix(0,n,length(fue_min)))
for (i in 1:length(fue_min)){
	set.seed(20)
	fue = runif(n, fue_min[i], fue_max[i]) # IPCC uncertainty level
	x_fue[,i]=fue
	names(x_fue) = c("P998","P999")
} 
fue_other = data.frame(matrix(1,n,1))

x_fue_final = data.frame (fue_other,fue_other,fue_other,fue_other,fue_other,fue_other, fue_other,fue_other,x_fue, x_fue[,2])

# Proportion of feed to residues
prop_res_max = c(0.72,0.20,0.46)
prop_res_min = c(0.60,0.18,0.38)
x_prop = data.frame(matrix(0,n,length(prop_res_max)))

for (i in 1:length(prop_res_max)){
	set.seed(20)
	prop = runif(n, prop_res_min[i], prop_res_max[i]) # IPCC uncertainty level
	x_prop[,i]=prop
	names(x_prop) = c("P995","P996", "P997")
} 
x_prop_final = data.frame (fue_other,fue_other,x_prop[,1],fue_other,fue_other,x_prop[,2], x_prop[,3],fue_other,fue_other,fue_other,fue_other)


my_total_crop = function (harvest, fue, propr,item) {
	area = data.frame(matrix (0, n, 11))
	area = harvest / (total_output_ha*fue*propr)
	total_c = data.frame(matrix (0, n, 1))
	total_c[,1] = rowSums(item * area /1000, na.rm=TRUE)				
	return (total_c)
}
# P input
crop_total_input= my_total_crop( harvest=herd_n_feed_use[,-12], item=total_input_ha, fue=x_fue_final, propr=x_prop_final)
# P output
crop_total_output = my_total_crop (harvest=herd_n_feed_use[,-12], item= total_output_ha, fue=x_fue_final, propr=x_prop_final)
# soil_stock change
crop_total_stock_change = my_total_crop (harvest=herd_n_feed_use[,-12], item= stock_change_ha, fue=x_fue_final, propr=x_prop_final)
# P losses
crop_total_loss = my_total_crop (harvest=herd_n_feed_use[,-12], item= loss_crop_ha, fue=x_fue_final, propr=x_prop_final)

#########################
#  MATRIX CONSTRUCTION ##
#########################

new_table = data.frame(matrix(0, n, 1))

# 1. Crop residues
crop_res_recycle =  my_total_crop(harvest = herd_n_feed_use[,-12], item= x_crop_resid, fue=x_fue_final, propr=x_prop_final) * 0.3   #add crop residues to resource

# 2. Manure
manu_recycle = my_total_crop (harvest = herd_n_feed_use[,-12], item= x_manure_f, fue=x_fue_final, propr=x_prop_final) * 0.25 #Based on manure use

# 3. organic waste
org_waste_recycle = new_table

# 4. feed intake
crop_feed_intake = herd_n_intake - data.frame(herd_n_feed_use[,12]/1000) 

# 5. feed from animal origin
anim_prod_feed  = new_table

# 6. feed from processing of Animal products
processed_prod_feed  = new_table

#7. crop_processing
crop_to_process  = new_table

# 8. Primary animal processing (live animal+milk)
anim_prod_process  = process_input

# 9. Animal product processing (0)
anim_pro_process = new_table

#Imported products

# 10. Imported crop
crop_import = new_table

# 11. feed to crop (0)
ani_to_crop = new_table

# 12. animal products to crop
processed_to_crop = new_table

# 13. crop to breeding
crop_to_breed = new_table

# 14. feed to breeding
anipro_to_breed = new_table

# 15. animal prod to breeding
ani_to_breed = new_table

# 16. fertilizer to processing =0
fert_to_process = new_table

# 17. Feed to processing = 0
feed_to_process = new_table

#18. animal prod to processing
ani_to_process = new_table

######################
#PRODUCT
#19.  Crop output (feed)
crop_output = crop_total_output 

# 20. animal to crop
crop_breeding = new_table

# 21. crop_process
crop_processing = new_table

#22. animal to cropping
ani_to_cropping = new_table

#23. Animal production output
anim_output = animal_output

#24. ani to animal process
ani_to_processing = new_table

#25. animal product to cropping
process_ani_crop = new_table

#26. Animal product to breeding
process_ani_breed = new_table

#27. Animal product to processing
process_output = process_output

####################
#Resources extraction

#28. Resource to cropping
resource_crop  = crop_total_input - crop_res_recycle  - manu_recycle

#29. resource to breeding
resource_anim = data.frame(herd_n_feed_use[,12]/1000)

#30. resource to processing
resource_process = new_table

##
#Stock change : Accumulation (negative); Remove (positive) - s matrix
#31. stock change crop
delta_crop  = crop_total_stock_change *-1

#32. stock change animal
delta_breed = new_table

#33. Stock change processing
delta_process = new_table

####
#Wastes
# 34 cropping
waste_crop = crop_total_loss # in Gg

# 35. breeding
waste_anim =  animal_loss

# 36. Processing
waste_process = process_loss

################################
##    INDICATORS ESTIMATION  ###
################################

###################
# Final model calculation
# Checking the errors in the model

my_verf_stage = function () {
	output_1 = data.frame(matrix(0, n, 3))

	for(i in 1:n) {
		u = matrix(c(crop_res_recycle [i,], manu_recycle [i,], org_waste_recycle [i,], crop_feed_intake [i,], anim_prod_feed [i,],
			   processed_prod_feed [i,], crop_to_process [i,], anim_prod_process [i,], anim_pro_process [i,]),
			 nrow=3, ncol=3, byrow=F)

		  #Resource
		r = matrix (c(resource_crop [i,], resource_anim [i,], resource_process [i,]), nrow=3, ncol=1, byrow=F)
		  #change in stock
		s = matrix (c(delta_crop [i,]*-1, delta_breed[i,]*-1, delta_process [i,]*-1), nrow=3, ncol=1, byrow=F) #-1 because in the table are -s'
		  #waste
		w = matrix (c(waste_crop [i,], waste_anim [i,], waste_process [i,]), nrow=3, ncol=1) #-1 because in the table are -w'
		  #import
		m = matrix (c(crop_import [i,], ani_to_crop [i,], processed_to_crop [i,] , crop_to_breed [i,],
						anipro_to_breed [i,], ani_to_breed [i,], fert_to_process [i,],
						feed_to_process [i,], ani_to_process[i,]), nrow=3, ncol=3, byrow=F)
		  #Product
		v = matrix (c(crop_output[i,], crop_breeding [i,], crop_processing [i,] , ani_to_cropping [i,],
				anim_output [i,], ani_to_processing [i,], process_ani_crop [i,], process_ani_breed [i,],
				process_output [i,] ), nrow=3, ncol=3, byrow=F)

		#Check validity
		mt=t(m)
		ut=t(u)
		verf_1 = rowSums(ut, na.rm=F)+rowSums(mt, na.rm=F)+t(r)
		verf_2 = rowSums(v, na.rm=F)+t(s)+t(w)

		output_1[i,] = verf_1 - verf_2
		}
	return (output_1)
}

verf_stage = my_verf_stage ()
verf_stage[1:10,]

###########

# NUE at each stage
my_nue_stage = function (){
	output_2 = data.frame(matrix(0, n, 3))
	names(output_2) = c("Nuec", "nuea", "nuep" )
	for(i in 1:n) {
		u = matrix(c(crop_res_recycle [i,], manu_recycle [i,], org_waste_recycle [i,], crop_feed_intake [i,], anim_prod_feed [i,],
			   processed_prod_feed [i,], crop_to_process [i,], anim_prod_process [i,], anim_pro_process [i,]),
			 nrow=3, ncol=3, byrow=F)

		  #Resource
		r = matrix (c(resource_crop [i,], resource_anim [i,], resource_process [i,]), nrow=3, ncol=1, byrow=F)
		  #change in stock
		s = matrix (c(delta_crop [i,]*-1, delta_breed[i,]*-1, delta_process [i,]*-1), nrow=3, ncol=1, byrow=F) #-1 because in the table are -s'
		  #waste
		w = matrix (c(waste_crop [i,], waste_anim [i,], waste_process [i,]), nrow=3, ncol=1) #-1 because in the table are -w'
		  #import
		m = matrix (c(crop_import [i,], ani_to_crop [i,], processed_to_crop [i,] , crop_to_breed [i,],
						anipro_to_breed [i,], ani_to_breed [i,], fert_to_process [i,],
						feed_to_process [i,], ani_to_process[i,]), nrow=3, ncol=3, byrow=F)
		  #Product
		v = matrix (c(crop_output[i,], crop_breeding [i,], crop_processing [i,] , ani_to_cropping [i,],
				anim_output [i,], ani_to_processing [i,], process_ani_crop [i,], process_ani_breed [i,],
				process_output [i,] ), nrow=3, ncol=3, byrow=F)

		#Check validity
		mt=t(m)
		ut=t(u)
		nue = (rowSums(v, na.rm=F)+t(s))/(rowSums(ut, na.rm=F)+rowSums(mt, na.rm=F)+t(r))
		output_2[i,] = nue*100
	}
	return (output_2)
}

nue_stage = my_nue_stage ()

windows (width=4, height=4)
boxplot(nue_stage)

#####################
#  Life-cycle NUE   #
#####################

my_life_cycle_nue = function (){
	output_3 = data.frame(matrix(0, n, 3))
	names(output_3) = c("r1", "r2", "r3" )

	for(i in 1:n) {
		u = matrix(c(crop_res_recycle [i,], manu_recycle [i,], org_waste_recycle [i,], crop_feed_intake [i,], anim_prod_feed [i,],
			   processed_prod_feed [i,], crop_to_process [i,], anim_prod_process [i,], anim_pro_process [i,]),
			 nrow=3, ncol=3, byrow=F)

		  #Resource
		r = matrix (c(resource_crop [i,], resource_anim [i,], resource_process [i,]), nrow=3, ncol=1, byrow=F)
		  #change in stock
		s = matrix (c(delta_crop [i,]*-1, delta_breed[i,]*-1, delta_process [i,]*-1), nrow=3, ncol=1, byrow=F) #-1 because in the table are -s'
		  #waste
		w = matrix (c(waste_crop [i,], waste_anim [i,], waste_process [i,]), nrow=3, ncol=1) #-1 because in the table are -w'
		  #import
		m = matrix (c(crop_import [i,], ani_to_crop [i,], processed_to_crop [i,] , crop_to_breed [i,],
						anipro_to_breed [i,], ani_to_breed [i,], fert_to_process [i,],
						feed_to_process [i,], ani_to_process[i,]), nrow=3, ncol=3, byrow=F)
		  #Product
		v = matrix (c(crop_output[i,], crop_breeding [i,], crop_processing [i,] , ani_to_cropping [i,],
				anim_output [i,], ani_to_processing [i,], process_ani_crop [i,], process_ani_breed [i,],
				process_output [i,] ), nrow=3, ncol=3, byrow=F)

		#Check validity
		mt=t(m)
		ut=t(u)
		gs = matrix (c(delta_crop [i,]*-1, 1:8*0),nrow=3, ncol=3, byrow=T)

		rf = t(r) %*% ginv (t(v) - u - m + gs ) #Calculation of RES* #Always use mass package to get ginv function
		output_3[i,] = rf
	}
	life_cycle = new_table
	life_cycle[,1] =(1/output_3[,3])*100

	return (life_cycle)
}

life_cycle_nue = my_life_cycle_nue ()


windows (width=4, height=4)
plot(density(life_cycle_nue[,1]))
#savePlot(paste('4- Graph/density_lcnue.emf', sep = ''), type = 'emf')

windows (width=4, height=4)
boxplot(life_cycle_nue)

windows (width=4, height=4)
hist(life_cycle_nue[,1], probability=TRUE, xlim=c(0,70))


############################
# Net Nutrient loss per ha #
############################

my_net_loss = function (harvest, feed_eff, prop_cr) {
	area = data.frame(matrix (0, n, 1))
	area[,1] = rowSums(harvest / (total_output_ha * feed_eff* prop_cr), na.rm=TRUE)
	net_loss = data.frame(matrix (0, n, 1))
	af_crop = (crop_res_recycle + crop_feed_intake)/crop_output
	af_animal = (manu_recycle + anim_prod_process)/anim_output
	af_crop = my_trunc_high(data=af_crop)
	af_animal = my_trunc_high(data=af_animal)
	net_loss[,1]= (waste_crop* af_crop + waste_anim*af_animal  + waste_process)*1000/area	
	return(net_loss)
}

net_nutrient_loss = my_net_loss (harvest= herd_n_feed_use[,-12], feed_eff= x_fue_final, prop_cr= x_prop_final)

###############
##    NHI  ####
###############

my_final_loss = function (){
	#Calculation of allocation factors based on N content
	af_crop = (crop_res_recycle + crop_feed_intake)/crop_output
	af_animal = (manu_recycle + anim_prod_process)/anim_output
	af_crop = my_trunc_high(data=af_crop)
	af_animal = my_trunc_high(data=af_animal)
	final_loss = data.frame(matrix(0, n, 3))
	names (final_loss) = c("loss_crop", "loss_animal", "loss_process")
	final_loss [,1] = waste_crop* af_crop
	final_loss [,2] = waste_anim* af_animal
	final_loss [,3] = waste_process
	return (final_loss)
}

final_loss_table = my_final_loss()

my_nhi = function (){
	mean_loss = data.frame(matrix(0, n, 1))
	mean_loss[,1] = rowMeans(final_loss_table, na.rm=TRUE)
	sd_loss = new_table
	sd_loss[,1] = apply (final_loss_table, 1, sd, na.rm=TRUE)
	nhi= data.frame(matrix(0, n, 1))
	names (nhi)[names(nhi)== "loss_crop"] = "NHI"
	nhi = (sd_loss/mean_loss)*100
	return (nhi)
}

nhi = my_nhi()

#########################
## OUTPUT FINAL RESULTS #
#########################

final_indicator = data.frame(life_cycle_nue, net_nutrient_loss, nhi)
write.table(final_indicator, "3- Results/Paper_3_results/result_UA_N_FADN_NL.csv", sep=",", row.names=FALSE)



##########################################
## sensitivity Analysis - Life-cycle NUE #
##########################################

h_nl = data.frame (x_repl_rate, x_afc, x_af_lw, x_am_lw,
		x_rf_lw, x_rm_lw, x_mf_lw, x_mm_lw, x_ca_lw, 
		x_ferti1, milk_year, x_milkfat, x_digest1, x_share_feed, x_proteinmilk, x_n_milk_replacer,
		feed_compo_n[,-c(7,8,12)], x_grazing_time,
		x_di_vol, x_indi_solid, x_indi_liquid, x_leaching_ma,
		x_carcassd, x_use_inedible, x_organic_waste,
		x_fertiliser[,-c(7,9,11)], x_manure[,-c(7,9,11)], x_n_organic, x_fixation[,-c(3,6,9,10)], 
		x_at_deposi[,c(1:6,8)], x_yield[,-c(7,9,11)], "P213"=ef_dir_fr[,1],
		"P214"=ef_dir_m [,1],"P215"=ef_ind_m [,1],"P216"=ef_ind_fr[,1], x_runoff[,c(1:4)],
		x_leaching[,c(1,3,5)], x_feed_loss, x_fue, x_prop)

corre_20nl = cor(h_nl)
write.table(corre_20nl, "3- Results/Paper_3_results/NL/corrfadn.csv", sep=",", row.names=FALSE)

y = with (h_nl, (my_life_cycle_nue ()))

x = src (X=h_nl, y=y,  nboot = 5000) 
print (x)

windows (width=4, height=4)
plot(x)
######

sensitivity_table = data.frame(as.matrix (x$SRC))
sensitivity_table$variable = row.names(sensitivity_table)
row.names(sensitivity_table) = NULL
new_data = data.frame(sensitivity_table$variable, sensitivity_table$original^2)  #Squared Regression coefficienct
names(new_data) = c("fact","coer_nue")
new_data = new_data[order(new_data[,2], decreasing = TRUE),]

windows (width=4, height=4)
barplot(new_data[c(1:20),2], xlab= "", ylab = "SRC", ylim = ,names = new_data[c(1:20),1],
		main= "Sensitivity analysis for Life-cycle NUE-N - Intensive")
#plot (new)
sum(new_data[,2])


##########################################
## sensitivity Analysis - Life-cycle NNB #
##########################################


h1_loss = data.frame (x_repl_rate, x_afc, x_af_lw, x_am_lw,
		x_rf_lw, x_rm_lw, x_mf_lw, x_mm_lw, x_ca_lw, 
		x_ferti1, milk_year, x_milkfat, x_digest1, x_share_feed, x_proteinmilk, x_n_milk_replacer,
		feed_compo_n[,-c(7,8,12)], x_grazing_time,
		x_di_vol, x_indi_solid, x_indi_liquid, x_leaching_ma,
		x_carcassd, x_use_inedible, x_organic_waste,
		x_fertiliser[,-c(7,9,11)], x_manure[,-c(7,9,11)], x_n_organic, x_fixation[,-c(3,6,9,10)], 
		x_at_deposi[,c(1:6,8)], x_yield[,-c(7,9,11)], "P213"=ef_dir_fr[,1],
		"P214"=ef_dir_m [,1],"P215"=ef_ind_m [,1],"P216"=ef_ind_fr[,1], x_runoff[,c(1:4)],
		x_leaching[,c(1,3,5)], x_feed_loss, x_fue, x_prop)


y_loss = with (h1_loss, (my_net_loss (harvest= herd_n_feed_use[,-12],  feed_eff= x_fue_final, prop_cr= x_prop_final)))

x_loss = src (X=h1_loss, y=y_loss,  nboot = 5000)
print (x_loss)

windows (width=4, height=4)
plot(x_loss)

sensitivity_loss = data.frame(as.matrix (x_loss$SRC))
sensitivity_loss$variable = row.names(sensitivity_loss)
row.names(sensitivity_loss) = NULL
new_data_loss = data.frame(sensitivity_loss$variable, sensitivity_loss$original^2)  #Squared Regression coefficienct
names(new_data_loss) = c("fact","coer_loss")
new_data_loss = new_data_loss[order(new_data_loss[,2], decreasing = TRUE),]

windows (width=4, height=4)
barplot(new_data_loss[c(1:20),2], xlab= "", ylab = "SRC", ylim = ,names = new_data_loss[c(1:20),1],
		main= "Sensitivity analysis for NNB-N - Intensive")
#plot (new)
sum(new_data_loss[,2])


####################################
## sensitivity Analysis  - NHI    ##
####################################

h1_nhi = data.frame (x_repl_rate, x_afc, x_af_lw, x_am_lw,
		x_rf_lw, x_rm_lw, x_mf_lw, x_mm_lw, x_ca_lw, 
		x_ferti1, milk_year, x_milkfat, x_digest1, x_share_feed, x_proteinmilk, x_n_milk_replacer,
		feed_compo_n[,-c(7,8,12)], x_grazing_time,
		x_di_vol, x_indi_solid, x_indi_liquid, x_leaching_ma,
		x_carcassd, x_use_inedible, x_organic_waste,
		x_fertiliser[,-c(7,9,11)], x_manure[,-c(7,9,11)], x_n_organic, x_fixation[,-c(3,6,9,10)], 
		x_at_deposi[,c(1:6,8)], x_yield[,-c(7,9,11)], "P213"=ef_dir_fr[,1],
		"P214"=ef_dir_m [,1],"P215"=ef_ind_m [,1],"P216"=ef_ind_fr[,1], x_runoff[,c(1:4)],
		x_leaching[,c(1,3,5)], x_feed_loss, x_fue, x_prop)


y_nhi = with (h1_nhi, (my_nhi()))

x_nhi = src (X=h1_nhi, y=y_nhi,  nboot = 5000)
print (x_nhi)

windows (width=4, height=4)
plot(x_nhi)

sensitivity_nhi = data.frame(as.matrix (x_nhi$SRC))
sensitivity_nhi$variable = row.names(sensitivity_nhi)
row.names(sensitivity_nhi) = NULL
new_data_nhi = data.frame(sensitivity_nhi$variable, sensitivity_nhi$original^2)  #Squared Regression coefficienct
names(new_data_nhi) = c("fact","coer_nhi")
new_data_nhi = new_data_nhi[order(new_data_nhi[,2], decreasing = TRUE),]

windows (width=4, height=4)
barplot(new_data_nhi[c(1:20),2], xlab= "", ylab = "SRC", ylim = ,names = new_data_nhi[c(1:20),1],
		main= "Sensitivity analysis for NHI-N_FADN")
#plot (new)
sum(new_data_nhi[,2])

# Final results table
result = data.frame( new_data ,new_data_loss,new_data_nhi)
write.table (result, "3- Results/Paper_3_results/result_all_parameter_N_FADN_NL.csv", sep=",", row.names=FALSE)

############################################
## GRAPHIC REPRESENTATIONS FINAL RESULTS  ##
############################################

int_table    = merge(new_data, new_data_loss, by.x = "fact", by.y = "fact", all.x=TRUE)

final_table  = merge(int_table, new_data_nhi, by.x = "fact", by.y = "fact", all.x=TRUE)
final_table  = final_table [order(final_table [,2], decreasing = TRUE),]
#final_table = final_table[1:15,]

f_table = as.matrix(final_table[,-1])
condition = f_table >= 0.01
sum_co = rowSums(condition)
final_resu = final_table [sum_co >0,]
#final_resul = rbind(final_resu[1:5,],  final_resu[9,], final_resu[c(6:8,10:11),])

#names_result = data.frame( c("P1", "P2", "P3", "P4", "P5", "P6",
							#"P7", "P8", "P9", "P10", "P11"))
#names(names_result) = c("Categ")

final_result = final_resu #data.frame(names_result, final_resul[,-1])
write.table(final_result, "3- Results/Paper_3_results/result_sensitivity_N_FADN_NL.csv", sep=",", row.names=FALSE)


windows (width=16, height=6)
par (mar=c(4,4,4,4)+0.1)
plot(final_result[,2] , xlab="Parameters", ylab="Squared standardized regression coefficients", ylim=c(0,0.70), col="red", cex=1.5, pch=13, xaxt = "n")
par(new=T)
plot(final_result[,3], xlab="", ylab="", ylim=c(0,0.70), col="blue", cex=1.5, pch=17, xaxt = "n")
par(new=T)
plot(final_result[,4], xlab="", ylab="", ylim=c(0,0.70), col="green", cex=1.5, pch=19, xaxt = "n")
axis(1, at = 1:nrow(final_result), labels = final_result$fact, cex = 1.5, las=3)
abline(v=c(1:nrow(final_result)), col="grey", lty=1)
abline(h=0.1, col="red", lty=1)
legend ("topleft", c(" Life-cycle NUE", " Life-cycle net nitrogen balance", " Nitrogen hotspot index"), cex = 1,
		pt.cex = 1.5, pch=c(13,17,19), bg = "white", col = c("red", "blue", "green"), border="black" , bty="o")
savePlot("3- Results/Paper_3_results/Nitrogen_SRC_FADN_NL_2010.emf",type="emf")



r_square = data.frame(matrix(c(sum(new_data[,2]), sum(new_data_loss[,2]), sum(new_data_nhi[,2]))), nrow=1, ncol=3)
names(r_square) = c("lcnue", "nnb", "nhi")
write.table(r_square, "3- Results/Paper_3_results/result_r_square_FADN.csv", sep=",", row.names=FALSE)




#######################################
#http://stackoverflow.com/questions/25070547/ggplot-side-by-side-geom-bar
library("reshape2")

 
fadn    = read.table("3- Results/Paper_3_results/NL/FADN/result_sensitivity_N_FADN_NL.csv", header = T, sep = ",", row.names = NULL) 


names(fadn) = c("para","NUE", "NNB", "NHI")


fadnm = melt(fadn)
windows (width=16, height=6)
ggplot(fadnm, aes(x = para, y= value, fill = variable), xlab="Input Parameters") +
	geom_bar(stat="identity", width=.7, position = "dodge")+
 	theme(axis.text = element_text(size=15))+
 	ylim(c(0,0.5))+
  	theme_bw(base_size = 12, base_family = "")+
 	theme(plot.margin = unit(c(1,1,1,1),"cm"))+
  	theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  savePlot("3- Results/Paper_3_results/ALL_Nitrogen_SRC_fadn_NL_2010.emf",type="emf")


# iNDIVIDUAL plot
####################
## Life-cycle-NUE ##
####################


fadn_nue = fadn[,1:2]

fadn_nue$nue2 <- order(fadn_nue$NUE, decreasing = TRUE)

windows (width=8, height=6)

ggplot(fadn_nue, aes( nue2, NUE)) +
	geom_bar(stat = "identity", fill="yellowgreen", colour="black")+
	scale_x_discrete("Para", breaks=sort(fadn_nue$nue2), labels=fadn_nue$para[order(fadn_nue$nue2)])+
	theme(axis.text = element_text(size=15))+
 	ylim(c(0,0.5))+
  	theme_bw(base_size = 12, base_family = "")+
 	theme(plot.margin = unit(c(1,1,1,1),"cm"))+
  	theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())
savePlot("3- Results/Paper_3_results/life_cycle_nue_fadn_NL.emf",type="emf")


####################
## Life-cycle-NNB ##
####################

fadn_nnb = fadn[,c(1,3)]
fadn_nnb = fadn_nnb[order(fadn_nnb$NNB, decreasing=TRUE),]


fadn_nnb$nnb2 <- order(fadn_nnb$NNB, decreasing = TRUE)

windows (width=8, height=6)

ggplot(fadn_nnb, aes( nnb2, NNB)) +
	geom_bar(stat = "identity", fill="tomato3", colour="black")+
	scale_x_discrete("Para", breaks=sort(fadn_nnb$nnb2), labels=fadn_nnb$para[order(fadn_nnb$nnb2)])+
	theme(axis.text = element_text(size=15))+
 	ylim(c(0,0.5))+
  	theme_bw(base_size = 12, base_family = "")+
 	theme(plot.margin = unit(c(1,1,1,1),"cm"))+
  	theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())
savePlot("3- Results/Paper_3_results/life_cycle_nnb_fadn_NL.emf",type="emf")


####################
## NHI ##
####################

fadnnhi = fadn[,c(1,4)]
fadnnhi = fadnnhi[order(fadnnhi$NHI, decreasing=TRUE),]


fadnnhi$nhi2 <- order(fadnnhi$NHI, decreasing = TRUE)

windows (width=8, height=6)

ggplot(fadnnhi, aes( nhi2, NHI)) +
	geom_bar(stat = "identity", fill="darkcyan", colour="black")+
	scale_x_discrete("Para", breaks=sort(fadnnhi$nhi2), labels=fadnnhi$para[order(fadnnhi$nhi2)])+
	theme(axis.text = element_text(size=15))+
 	ylim(c(0,0.5))+
  	theme_bw(base_size = 12, base_family = "")+
 	theme(plot.margin = unit(c(1,1,1,1),"cm"))+
  	theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())
savePlot("3- Results/Paper_3_results/nhi_fadn_NL.emf",type="emf")

