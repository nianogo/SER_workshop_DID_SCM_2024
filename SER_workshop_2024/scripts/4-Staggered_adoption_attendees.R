#------------------------------------------------------------------------------#
#-----An overview of Difference-in-Difference and Synthetic Control Methods----- 
#-------------------------Classical and Novel Approaches-----------------------#
#-----------Society for Epidemiologic Research (SER)---------------------------# 
#---------------------------------Workshop: Part 4/4---------------------------#
#-----------------------------Date:04/22/24------------------------------------#
#Roch Nianogo (niaroch@ucla.edu) & Tarik Benmarhnia (tbenmarhnia@health.ucsd.edu)
#------------------------------------------------------------------------------#

#------------------------------------Checking the directory---------------------
getwd()


#-------------------------------------Installing packages-----------------------


if (!require("pacman")){
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
} # a nice package to load several packages simultaneously


p_load("tidyverse","magrittr","broom",        #manipulate data
       "cleaR",                               #clear workspace 
       "here",                                #directory managment
       "Synth", "gsynth",                     #synthetic control
       "panelView", "lme4", "estimatr",       #multi-level modeling
       "did"                                  #for did w/ staggered adoption
)                           

#------------------Mise en place (clear everything)------------------------------
clear() #same as remove(list=ls())




#------------------------------------Load the data------------------------------

sim_data_hte_staggered <- read_csv(here("data", "sim_data_hte_staggered.csv"))



#-------------------------Some info about the data------------------------------

# 50 states
# 15 are treated, 35 are controls
# staggered adoption

  # state.name[1:5] ("Alabama", "Alaska", "Arizona", "Arkansas", "California") 
    #enacted policy in 2000

  # state.name[6:10] ("Colorado", "Connecticut", "Delaware", "Florida", "Georgia") 
    #enacted policy in 2003

  # state.name[6:10] ("Hawaii", "Idaho", "Illinois", "Indiana", "Iowa") 
    #enacted policy in 2006

#-----------------------------------Visualize the data--------------------------
mydata = sim_data_hte_staggered

p_load("panelView")
panelview(y ~ treatedpost, data = mydata,
          index = c("state","year"), 
          pre.post = TRUE) 


panelview(y ~ treatedpost, data = mydata, 
          index = c("state","year"), 
          type = "outcome",  
          by.group = TRUE)


#For all the data
mydata %>% 
  ggplot(aes(x=year, y=y, group=state)) + 
  annotate("rect", fill = "gray", alpha = 0.5,
           xmin = 2000, xmax = Inf,
           ymin = -Inf, ymax = Inf) +
  labs(title = paste("Outcome by year"),
       x = "Year", 
       y = "Outcome",
       color = "Treatment") +
  geom_line(aes(color=factor(treated)), size=0.5) +
  scale_color_discrete(labels=c("Control", "Treated")) +
  geom_vline(xintercept = 2000, lty=2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 



#---------------Analyze the data: Callaway and Sant’Anna------------------------
#Tutorial here: https://cran.r-project.org/web/packages/did/vignettes/did-basics.html
#Paper here: https://www.sciencedirect.com/science/article/pii/S0304407620303948?via%3Dihub

#load package
p_load(did)

#Create new variables
mydata1 <- mydata %>% 
  mutate(first_treated = case_when( (state %in% state.name[1:5])   ~ 2000,
                                    (state %in% state.name[6:10])  ~ 2003,
                                    (state %in% state.name[11:15]) ~ 2006,
                                    TRUE~0))



group_time_effects <- att_gt( yname  = "y",
                              tname  = "year",
                              idname = "state_num",
                              gname  = "first_treated",
                              xformla = ~ xit,
                              data = mydata1)


#Group-time effects
summary(group_time_effects)
ggdid(group_time_effects)

#Simple Aggregation
agg.simple <- aggte(group_time_effects, type = "simple")
summary(agg.simple)

#Dynamic Effects and (Event Studies): Effect by length of exposure
agg.es <- aggte(group_time_effects, type = "dynamic")
summary(agg.es)
ggdid(agg.es)

#Group-Specific Effects : : Effect by group
agg.gs <- aggte(group_time_effects, type = "group")
summary(agg.gs)
ggdid(agg.gs)

#Calendar Time effects
agg.ct <- aggte(group_time_effects, type = "calendar")
summary(agg.ct)
ggdid(agg.ct)

# Now that we have seen different effects that could be obtained using the method 
# above, let us see what other models give us.

#-------------------------Analyze the data: Simple DID--------------------------


p_load("estimatr")
dta <- lm_robust(y ~ treatedpost + factor(year) + xit, 
                 data = mydata,
                 fixed_effects=state,
                 clusters = state, 
                 se_type = "stata")
dta

did <- round(data.frame(ATT     = dta$coefficients["treatedpost"], 
                        se      = dta$std.error["treatedpost"],
                        low_ci  = dta$conf.low["treatedpost"],
                        high_ci = dta$conf.hig["treatedpost"]),2)
did


p_load("lmerTest", "lme4")
fit1 <- lmerTest::lmer(y ~ treated:post + 
                         xt + xi + xit + 
                         (1| state) + 
                         (1 | year), 
                       data = mydata, 
                       REML = T)
broom.mixed::tidy(fit1) %>% mutate_if(is.numeric, round, 2) %>% 
  flextable::flextable()
# This quantity above is not giving us the group-specific overall ATT and could also
# be biased.  


#----------------------Analyze the data: Generalized SCM------------------------


y <- gsynth(y ~ treatedpost + xit, 
            data = mydata,  
            EM = F, 
            index = c("state","year"), 
            inference = "parametric", 
            se = TRUE,
            nboots = 100,  #so that it can run faster, default is 200
            r = c(0, 5), 
            CV = TRUE, 
            seed = 123,
            force = "two-way", 
            parallel = FALSE)

y1 <- round(data.frame(y$est.avg),2)


#Period-specific ATT
y$est.att

#average ATT
y$est.avg

#Weights from the generalized synthetic controls
y[["alpha.co"]]

plot(y, type = "counterfactual", raw = "none", main="")

plot(y, type = "counterfactual", raw = "band", main="")
plot(y, type = "counterfactual", raw = "all")
plot(y, type = "ct", raw = "none", main = "", shade.post = FALSE)



#-----Estimate individual state (or time) effects and then Pooled the results----
#Using DID methods
#Using the generalized SCM
#Using the augmented SCM


