### Excercise 11.1 Baby response to his name
rm(list=ls())
options(digits=3, width=60)
library(dplyr)
library(rstan)
library(brms)
library(ggplot2)
library(ape)
library(caper)
library(matrixStats)
library(sjPlot)
library(rstan)
rstan_options(auto_write = TRUE)
library(performance)





# Inspired by Rasmus Baath's Swedish Fish problem. 
# A mother wants to know whether her baby is responding to his name at month 7 (a development milestone).

rnorm(1)
rpois(2, 1)
# build a bayesian model 
# part 1: data
n_trail <- 10 
n_success <- 6

# part 2: generative model, the story about how my data was generated, or which distribution should my response variable follow?
# Hint: similar problem as coin toss, what's the odds of me getting a head/tail?







 
generative_model <- function(response_rate) {
    return(rbinom(1,n_trail, response_rate))
}

# part 3: set the prior for the parameters we are estimating
# what is it?
# what do we know about it?









# set prior 
prior_size <- 20000
prior <- runif(prior,0,1) # Here you sample prior_size draws from the prior  
par(mfrow=c(2,1))
hist(prior) # It's always good to eyeball the prior to make sure it looks ok.


# All ingredients ready! Let's fit our model (infer our parameter) in a bayesian fashion
# create a dataframe to store the simulated number of responses
sim_data <- rep(NA, prior_size)
# Here we put each prior into our generative model to simulate one case of the response (y)
for(i in 1:prior_size) {
    sim_data[i] <- generative_model(prior[i])
}
hist(sim_data)
# We only want to keep the priors that generated our actual data. Here you filter off all draws that do not match the data.
n_success = 6
#
sim_data == n_success
posterior <- prior[sim_data == n_success] 
length(posterior)
hist(posterior) # Eyeball the posterior
length(posterior) # See if we got enough draws left after the filtering.
# There are no rules here, but you probably want to aim
# for >1000 draws. 

# Now you can summarize the posterior, where a common summary is to take the mean
# or the median posterior, and perhaps a 95% quantile interval.
median(posterior)
mean(posterior)

# Question: Is there a better chance of the baby responding to his name than not to?          

# Way 1: evidence ratio (Evid.Ratio). For a one-sided hypothesis, this is just the posterior probability (Post.Prob) under the hypothesis against its alternative. That is, when the hypothesis is of the form a > b, the evidence ratio is the ratio of the posterior probability of a > b and the posterior probability of a < b.

higher <- length(posterior[posterior > 0.5])/length(posterior)    

# Way 2: Credible interval

# Is the baby's response significantly different from 0.5?
quantile(posterior, c(0.025, 0.975))






                                                                                                                                   

### Question III) If I call the baby 100 times, what would be number of responses?

#**Hint:** The answer is again not a single number but a distribution over probable number of responses.











# rbinom is vectorized, we can give it a vector of probabilities:
responses <- rbinom(n = length(posterior), size = 100, prob = posterior)

hist(responses, xlim = c(0, 100))
quantile(responses, c(0.025, 0.975))
# So a decent guess is that is would be between 30 and 85 responses.

# Back to slides

## Excercise 11.2 Expert advice from grandma

# set prior 
prior_size <- 10000
prior <- runif(prior_size,0.5,1) 
hist(prior) 


# All ingredients ready! Let's fit our model (infer our parameter) in a bayesian fashion
# create a dataframe to store the simulated number of responses
sim_data <- rep(NA, prior_size)
# Here we put each prior into our generative model to simulate one case of the response (y)
for(i in 1:prior_size) {
    sim_data[i] <- generative_model(prior[i])
}

# We only want to keep the priors Here you filter off all draws that do not match the data.
n_success = 6
posterior <- prior[sim_data == n_success] 

hist(posterior) # Eyeball the posterior
length(posterior) # See that we got enought draws left after the filtering.
# There are no rules here, but you probably want to aim
# for >1000 draws. 
median(posterior)
mean(posterior)
quantile(posterior, c(0.025, 0.975))

# advice from american pediatric association
# set prior 
prior_size <- 10000
prior <- rbeta(prior_size,7,3) # play with rbeta to see what shape1 shape2 does to the distribution, we want prior centers around 70-80%. 
hist(prior) # It's always good to eyeball the prior to make sure it looks ok.


# All ingredients ready! Let's fit our model (infer our parameter) in a bayesian fashion
# create a dataframe to store the simulated number of responses
sim_data <- rep(NA, prior_size)
# Here we put each prior into our generative model to simulate one case of the response (y)
for(i in 1:prior_size) {
    sim_data[i] <- generative_model(prior[i])
}

# We only want to keep the priors Here you filter off all draws that do not match the data.
n_success = 6
posterior <- prior[sim_data == n_success] 

hist(posterior) # Eyeball the posterior
length(posterior) # See that we got enought draws left after the filtering.
# There are no rules here, but you probably want to aim
# for >1000 draws. 
median(posterior)
mean(posterior)
quantile(posterior, c(0.025, 0.975))

### back to slides

### Excercise 11.2 Baby vs. Boarder Collie 
### Question: Is my boarder collie more responsive to his name than my baby?

# part 1: data
n_trail <- 10 
n_baby <- 6
n_dog <- 9

# part 2: generative model, the story about how my data was generated, or which distribution should my response variable follow?
# Hint: similar problem as coin toss, what's the odds of me getting a head/tail?

generative_model <- function(response_rate) {
    return(rbinom(1,n_trail, response_rate))
}

# part 3: set the prior for the parameters we are estimating
# what is it?
# what do we know about it?

# set prior 
prior_size <- 20000
prior_baby <- runif(prior_size,0,1)  
prior_dog <- runif(prior_size,0,1) 

# All ingredients ready! Let's fit our model (infer our parameter) in a bayesian fashion
# create a dataframe to store the simulated number of responses
sim_data <- data.frame('baby'=rep(NA, prior_size),'dog'=rep(NA, prior_size),
                       'rate_baby' = prior_baby, 'rate_dog' = prior_dog)
# Here we put each prior into our generative model to simulate one case of the response (y)
for(i in 1:prior_size) {
    sim_data$baby[i] <- generative_model(prior_baby[i])
    sim_data$dog[i] <- generative_model(prior_dog[i])
}

# We only want to keep the priors Here you filter off all draws that do not match the data.

posterior <- sim_data %>% filter((sim_data$baby == n_baby) & (sim_data$dog == n_dog)) 
# plot the posterior
hist(posterior$rate_baby) # Eyeball the posterior
hist(posterior$rate_dog) # Eyeball the posterior
# Is Bingo more responsive than Xiaobai?
posterior$diff <- posterior$rate_baby - posterior$rate_dog
# Now you can summarize the posterior, where a common summary is to take the mean
# or the median posterior, and perhaps a 95% quantile interval.
median(posterior$diff)
mean(posterior$diff)
quantile(posterior$diff, c(0.025, 0.975))
mean(posterior$diff<0)
nrow(posterior) # See that we got enought draws left after the filtering.
# There are no rules here, but you probably want to aim
# for >1000 draws. 

# Very inefficient!
plot(sim_data$rate_baby,sim_data$rate_dog)
plot(posterior$rate_baby,posterior$rate_dog)

# We need more effective ways to search through the parameter space!

##########back to slides#############

## 

## R code for in-class test example in Ravenzwaaij 2018.
# ?dExample: inâ€“class test



# Learning stan!
# redo the baby response problem with rstan
library(rstan)

# The Stan model as a string.
model_string <- "
// Here we define the data we are going to pass into the model
data {
int n; // Number of trials
int s;  // Number of successes
}

// Here we define what 'unknowns' aka parameters we have.
parameters {
real<lower=0, upper=1> rate;
}

// The generative model
model {
rate ~ uniform(0, 1);
s ~ binomial(n, rate);
}

// In the generated quantiles block you can calculate 'derivatives' of
// the parameters. Here is a silly example calculating the square of the 
// rate. Variables have to be defined before they are assigned to.
generated quantities {
real rate_sqrt;
rate_sqrt = rate^0.5;
}
"

data_list <- list(n = 10, s = 6)

# Compiling and producing posterior samples from the model.
stan_samples <- stan(model_code = model_string, data = data_list)

# Plotting and summarizing the posterior distribution

stan_samples
# lp: "log density up to a constantâ€? in the manual. Here, log density should be log likelihood of the observations conditioned on the posterior parameters: p(y | p_post). In the MCMC sampling paradigm, this could be obtained by plugging the parameter samples in each iteration into the likelihood and compute the probability of each observation. Therefore, â€œlpâ€? is actually a vector in length equals to the number of observations. â€œlp__â€? in Stanâ€™s output would be the sum of the â€œlpâ€? vector, which essentially quantifies how well the model match the data, hereby useful for model comparison purposes.

rstan::traceplot(stan_samples)
plot(stan_samples)

# Export the samples to a data.frame for easier handling.
posterior <- as.data.frame(stan_samples)

# Now we could, for example, calculate the probability that the rate is higher
# than, say, 20%
sum(posterior$rate > 0.2) / length(posterior$rate )

# More excercise
## Bayesian A/B testing for baby vs. dog with Stan

# The Stan model as a string.
model_string <- "
// Here we define the data we are going to pass into the model
data {
    int n; // Number of trials
    int sbaby;  // Number of successes for baby
    int sdog; // Number of successes for dog
}

// Here we define what 'unknowns' aka parameters we have.
parameters {
    real<lower=0, upper=1> rateB;
    real<lower=0, upper=1> rateD;
}

// The generative model
model {
    rateB ~ uniform(0, 1);
    rateD ~ uniform(0, 1);
    sbaby ~ binomial(n, rateB);
    sdog ~ binomial(n, rateD);
}

// In the generated quantiles block you can calculate 'derivatives' of
// the parameters. Here is a silly example calculating the square of the 
// rate. Variables have to be defined before they are assigned to.
generated quantities {
    real rate_diff;
    rate_diff = rateD - rateB;
}
"

data_list <- list(n = 10, sbaby = 6, sdog = 9)

# Compiling and producing posterior samples from the model.
stan_samples <- stan(model_code = model_string, data = data_list)

# Plotting and summarizing the posterior distribution
stan_samples
rstan::traceplot(stan_samples)
plot(stan_samples)

# Export the samples to a data.frame for easier handling.
posterior <- as.data.frame(stan_samples)

# Now we could, for example, calculate the probability that the rate is higher
# than, say, 20%
sum(posterior$rate_diff > 0.2) / length(posterior$rate_diff )



### Question II: Change the model so that it uses a more informative prior. 

### code data in another way

library(rstan)

# The Stan model as a string.
model_string <- "
data {
# Number of data points
int n1;
int n2;
# Number of successes
int y1[n1];
int y2[n2];
}

parameters {
real<lower=0, upper=1> theta1;
real<lower=0, upper=1> theta2;
}

model {  
theta1 ~ beta(1, 1);
theta2 ~ beta(1, 1);
y1 ~ bernoulli(theta1);
y2 ~ bernoulli(theta2); 
}

generated quantities {
}
"

y1 <- c(1, 1, 1, 1, 1, 0, 1, 0, 0, 0)
y2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
data_list <- list(y1 = y1, y2 = y2, n1 = length(y1), n2 = length(y2))
# Compiling and producing posterior samples from the model.
stan_samples <- stan(model_code = model_string, data = data_list)

# Plotting and summarizing the posterior distribution
stan_samples
# log likelihood of the observations conditioned on the posterior parameters: p(y | p_post). 
# https://www.jax.org/news-and-insights/jax-blog/2015/october/lp-in-stan-output#
plot(stan_samples)

#using brms 
# Setting up data
y1 <- c(1, 1, 1, 1, 1, 0, 1, 0, 0, 0)
y2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
data_list <- list(y1 = y1, y2 = y2, n1 = length(y1), n2 = length(y2))
data <- data.frame('response' = c(y1,y2),
                   'name'=c(rep('baby',length(y1)),rep('dog',length(y2))))

# visualize the data
ggplot(data, aes(x = name, y = response)) + 
    geom_point() +
    geom_jitter() 

## generalized linear model
glm1 <- glm(response ~ name, data = data, family = "binomial")
summary(glm1)

## brms with default setting
model_brms <- brm(response~name,
                  data = data, 
                  family = bernoulli(),
                  chains = 4, # default settings, 
                  iter = 2000, # default settings, 
                  warmup = 1000, # default settings, 
                  thin = 1, # default setting, no thinning. Must be a positive integer. Set thin > 1 to save memory and computation time if iter is large.
                  ) 

summary(model_brms)

# Rhat and Effective Sample Size (ESS) are convergence and efficiency diagnostics for Markov Chains
#Rhat < 1.01, the estimation is convergent
# Bulk_ESS: bulk of the distribution (related e.g. to efficiency of mean and median estimates)
# tail ESS: computing the minimum of effective sample sizes for 5% and 95% quantiles
# Both bulk-ESS and tail-ESS should be at least 100 (approximately) per Markov Chain in order to be reliable and indicate that estimates of respective posterior quantiles are reliable.

## plot the model for diagnostics.
plot(model_brms)
# left column: distribution of posterior for model parameters. Q: is it normal?
# right colum: it's value on each chain. Q: why only 1000 not 2000?
## summarizing the posterior: box plot
mcmc_plot(model_brms)
# What does this mean?

## Hypothesis testing

hypothesis(model_brms, "namedog > 0",
           class = 'b')

# Prediction: plot the posterior fixed effect estimates with credible intervals over the data
plot(conditional_effects(model_brms), points = TRUE)
# check the usage, more useful in a mixed model, recall conditional and marginal r-squared.
?conditional_effects

#R-squared: variance explained by model
bayes_R2(model_brms)

library(sjPlot)
tab_model(model_brms)

## More informative priors
# get the default priors
prior_summary(model_brms)
# student_t, df degrees of freedom, mu and sigma(the later two are non-centrality optional parameters).
# vecterization: When putting the same prior on all population-level effects (fixed effects for categorical variable) at once, those priors can be vectorized (sampled at the same time) and the sampling process is faster.

# Usually we focus on priors for the regression coefficients and not on the error and variance terms, since we are most likely to actually have information on the size and direction of a certain effect and less (but not completely) unlikely to have prior knowledge on the unexplained variances. 
model_brms_1 <- brm(response~name,
                  data = data, 
                  family = bernoulli(),
                  chains = 4, # default settings, 
                  iter = 2000, # default settings, 
                  warmup = 1000, # default settings, 
                  thin = 1, # default setting
                  prior = prior(normal(1,1), class="b",coef="namedog")
) 
summary(model_brms_1)
summary(model_brms)
## Continuous data, milk and cow
# Earlier this year farmer John ran an experiment where he gave 10 cows a special diet that he had heard could make them produce more milk. He recorded the number of liters of milk from these "diet" cows and from 15 "normal" cows during one month.
library(rstan)

diet_milk <- c(651, 679, 374, 601, 401, 609, 767, 709, 704, 679)
normal_milk <- c(798, 1139, 529, 609, 553, 743, 151, 544, 488, 555, 257, 692, 678, 675, 538)

# The Stan model as a string.
model_string <- "
data {
# Number of data points
int n1;
int n2;
# the number of liters of milk
vector[n1] v1;
vector[n2] v2;
}

parameters {
real<lower=0, upper=2000> mu1;
real<lower=0, upper=2000> mu2;
real<lower=0, upper=2000> sigma1;
real<lower=0, upper=2000> sigma2;
}

model {  
mu1 ~ normal(mean(v1), mean(v1));
mu2 ~ normal(mean(v2), mean(v2));
sigma1 ~ normal(sd(v1), sd(v1));
sigma2 ~ normal(sd(v2), sd(v2));
v1 ~ normal(mu1,sigma1);
v2 ~ normal(mu2,sigma2);
}

generated quantities {
}
"

data_list <- list(v1 = diet_milk, v2 = normal_milk, n1 = length(diet_milk), n2 = length(normal_milk))

# Compiling and producing posterior samples from the model.
stan_samples <- stan(model_code = model_string, data = data_list)

# Plotting and summarizing the posterior distribution
stan_samples
# log likelihood of the observations conditioned on the posterior parameters: p(y | p_post). 
# https://www.jax.org/news-and-insights/jax-blog/2015/october/lp-in-stan-output#
plot(stan_samples)

# milk and cow in different format
d <- data.frame(
    milk = c(651, 679, 374, 601, 401, 609, 767, 709, 704, 679, 798, 1139,
             529, 609, 553, 743, 151, 544, 488, 555, 257, 692, 678, 675, 538),
    group = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 
              2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2))

data_list <- list(y = d$milk, x = d$group, n = length(d$milk), 
                  n_groups = max(d$group))

model_string <- "
data {
# Number of data points
int n;
int n_groups;
# the number of liters of milk
vector[n1] y;
vector[n2] v2;
}

parameters {
real<lower=0, upper=2000> mu1;
real<lower=0, upper=2000> mu2;
real<lower=0, upper=2000> sigma1;
real<lower=0, upper=2000> sigma2;
}

model {  
mu1 ~ normal(mean(v1), mean(v1));
mu2 ~ normal(mean(v2), mean(v2));
sigma1 ~ normal(sd(v1), sd(v1));
sigma2 ~ normal(sd(v2), sd(v2));
v1 ~ normal(mu1,sigma1);
v2 ~ normal(mu2,sigma2);
}

generated quantities {
}
"""

# Chicken and eggs
model_string <- "
data {
# Number of data points
int n1;
int n2;
# the number of eggs
int y1[n1];
int y2[n2];
}

parameters {
real<lower=0, upper=14> lambda1;
real<lower=0, upper=14> lambda2;
}

model {  
y1 ~ poisson(lambda1);
y2 ~ poisson(lambda2);
}

generated quantities {
}
"

diet_eggs <- c(6, 4, 2, 3, 4, 3, 0, 4, 0, 6, 3)
normal_eggs <- c(4, 2, 1, 1, 2, 1, 2, 1, 3, 2, 1)

data_list <- list(y1 = diet_eggs, y2 = normal_eggs, n1 = length(diet_eggs), n2 = length(normal_eggs))

# Compiling and producing posterior samples from the model.
stan_samples <- stan(model_code = model_string, data = data_list)

# Plotting and summarizing the posterior distribution
stan_samples
# log likelihood of the observations conditioned on the posterior parameters: p(y | p_post). 
# https://www.jax.org/news-and-insights/jax-blog/2015/october/lp-in-stan-output#
plot(stan_samples)
```
#**Hint 1:** If you have a vector of samples representing a probability distribution, which you should have from the last question, calculating the amount of probability above a certain value is done by simply *counting* the number of samples above that value and dividing by the total number of samples.
#**Hint 2:** when the response_rate is 0.5, there is equal chance of response vs. no response.







