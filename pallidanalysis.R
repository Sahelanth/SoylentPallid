#Is eating Soylent 1.5 for every meal more dangerous than eating normal foods,
#in terms of lead and cadmium exposure?

#Background:
#You've probably heard of Soylent. A few years back, I was eating a lot of Soylent 1.5.
#It came out that 1.5 had lead and cadmium levels above California thresholds:
#https://www.eater.com/2015/8/16/9162301/soylent-as-you-sow-cadmium-lead-california
#Now, I knew California has some overly cautious laws re cancer, so I initially wasn't *that* concerned.

#Then I read Rosa Labs (Soylent's manufacturer's) response
#https://faq.soylent.com/hc/en-us/articles/204197379-California-Proposition-65
#It seemed obvious to me that the meals they cited as a comparison, heavy on fish, were chosen to make
#1.5 look better and unlikely to be comparable to the metal levels in what people normally eat. So I
#decided to investigate.

##Plan
#My original analysis plan was to compare a 3-soylent-meal day to 3 "normal" meals constructed from other 
#foods. Use hierarchical clustering to see if 3 Soylent meals 1.5 naturally "falls in the same category" 
#as some other diet plan, or whether it's an outlier.

#My final plan was to compare the mean and median exposures like so:
#Single meal, lead
#Single meal, cad
#3 meals mean, lead
#3 meals mean, cad
#3 meals median, lead
#3 meals median, cad
#and do Bonferroni to test for multiple comparisons.

#Key question is qualitative - whether soylent 3x per day has comparable content to other foods.
#Whether it's in the same cluster as a bunch of others, or in a cluster of its own.


##I abandoned that analysis plan because it was obvious from the one-meal results that the rest
#weren't necessary.

##############################Extract FDA food rows

#Get data
download.file("http://www.fda.gov/downloads/Food/FoodScienceResearch/TotalDietStudy/UCM455198.zip", destfile="./yearstandards.zip")

stds <- read.table(unz("yearstandards.zip", "Elements 2011.txt"), header=TRUE, sep="\t")

#Focus on only the rows containing lead or cadmium
leadandcad <- subset(stds, Element=="CADMIUM" | Element=="LEAD")

#Check whether any entries are in micrograms rather than milligrams
grep("ug/kg", leadandcad$Unit)
#None are, so don't need to adjust concentrations

#Variables I care about: Food.No, Food.Name, Element, Conc

averaged <- aggregate(Conc ~ Food.Name + Element, leadandcad, mean)

###############################add a row for Soylent 1.5.
#Manually add in lead/ cad levels from here: 
#https://faq.soylent.com/hc/en-us/article_attachments/203065059/prop-65-meals.v2.rev_D.jpg

#Make data frame for Soylent using Rosa Labs concentrations from https://docs.google.com/spreadsheets/d/1xS8bAQKZoksJfONrEb1N_kFD5Sp-0W-TOJBOJQyLJk8/edit#gid=0

Food.Name <- c("Soylent1.5", "Soylent1.5", "Soylent2.0", "Soylent2.0")
Element <- c("CADMIUM", "LEAD", "CADMIUM", "LEAD")
Conc <- c(0.186, 0.043, 0.010, 0.005)
soyframe <- data.frame(Food.Name, Element, Conc)

averaged <- rbind(averaged, soyframe)

#Okay, now I should split into lead and cadmium frames, presumably using aggregate().

leadset <- averaged[averaged$Element=="LEAD",]
cadset <- averaged[averaged$Element=="CADMIUM",]

#Simplify these sets into just food name and concentration

simpleleadset <- leadset[,3]
names(simpleleadset) <- leadset$Food.Name

simplecadset <- cadset[,3]
names(simplecadset) <- cadset$Food.Name

#I should reshape averaged into having lead conc in one column and cadmium conc in another, and use them as x and y
library(tidyverse)
data_for_hclust <- averaged %>% spread(key=Element, value=Conc)
row.names(data_for_hclust) <- make.names(data_for_hclust[,1])
data_for_hclust <- data_for_hclust[,-1]

#Now take pairwise distances
pairwise_distances <- dist(data_for_hclust)
hclustering <- hclust(pairwise_distances)
lead_cad_pairwise_distance_clusters <- plot(hclustering, cex=0.5)
#Holy shit is there an outlier. Has rows 248 and 272.
#That's Soylent 1.5 and sunflower seads, both because of cadmium.

pairwise_distances_lead_only <- dist(simpleleadset)
hclustering_lead_only <- hclust(pairwise_distances_lead_only)
plot(hclustering_lead_only)
#WELP

pairwise_distances_cad_only <- dist(simplecadset)
hclustering_cad_only <- hclust(pairwise_distances_cad_only)
plot(hclustering_cad_only, cex=0.5)


###So, average case is not looking good for soylent. What if aggregate the max, rather than average,
###levels of lead and cadmium found in foods?
max <- aggregate(Conc ~ Food.Name + Element, leadandcad, FUN=max)
max <- rbind(max, soyframe)

#Okay, now I should split into lead and cadmium frames, presumably using aggregate().

maxleadset <- max[max$Element=="LEAD",]
maxcadset <- max[max$Element=="CADMIUM",]

#Simplify these sets into just food name and concentration

simplemaxleadset <- maxleadset[,3]
names(simplemaxleadset) <- maxleadset$Food.Name

simplemaxcadset <- maxcadset[,3]
names(simplemaxcadset) <- maxcadset$Food.Name

max_data_for_hclust <- max %>% spread(key=Element, value=Conc)
row.names(max_data_for_hclust) <- make.names(max_data_for_hclust[,1])
max_data_for_hclust <- max_data_for_hclust[,-1]

max_pairwise_distances <- dist(max_data_for_hclust)
max_hclustering <- hclust(max_pairwise_distances)
max_lead_cad_pairwise_distance_clusters <- plot(max_hclustering, cex=0.5)
#Here it's no longer in the highest cluster with sunflower seeds... but it's still an outgroup to
#everything except sunflower seeds

pairwise_distances_max_lead_only <- dist(simplemaxleadset)
hclustering_max_lead_only <- hclust(pairwise_distances_max_lead_only)
plot(hclustering_max_lead_only)
#Still outlier, though closer to things like raw avocado and canned fruit - can discuss that

pairwise_distances_max_cad_only <- dist(simplemaxcadset)
hclustering_max_cad_only <- hclust(pairwise_distances_max_cad_only)
plot(hclustering_max_cad_only, cex=0.5)
#Potato chips and boiled spinach also relatively bad on cadmium, but still far from normal


###################################Will it hurt you?
#Now, three-meal analysis is beside the point. Nothing but eating all-sunflower-seed meals will match it.
#But key question: will it hurt you? 
#Find out mass of 3 soylent 1.5 powder meals per day, multiply by concentration, see if it takes you over the limit.


#files.soylent.com/pdf/soylent-nutrition-facts.pdf
#115 grams per serving, 4 servings per pouch.
#460 grams per day -> 0.46 kg

#Cadmium intake per day 
0.46*max_data_for_hclust[272,1]
#That's 0.08556 mg

#Lead:
0.46*max_data_for_hclust[272,2]
#That's 0.01978 mg


##FDA limits:
#Limits for solid food - lead 6 ppm, liquid food 1 ppm
#Cadmium - 0.1 ppm for "cereals and vegetables"
#The Prop 65 safe harbor level for lead is 0.5 micrograms per day. 60x lower than avg exposure EPA permits in drinking water -> amount they found in Soylent was 1/5 to 1/2 lower than what EPA permits in drinking water

#mg/kg converts 1:1 to ppm, as you can see through dimensional analysis.
#So, lead concentration is 1/2 what FDA allows in liquid food, and is well over the (more cautious) Prop 65
#safe harbor level.
#The cadmium level is higher than the FDA cereals-and-vegetables limit, and if liquid makes it more
#bioavailable that is a real problem.


