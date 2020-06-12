

load("~/Desktop/panel_change/data/gsslong.Rdata")
gss_vars <- names(gss.long)[6:length(names(gss.long))]
gss_vars <- gss_vars[gss_vars %!in% c("colmslm", "libmslm", "spkmslm",
                                      "mslmscale")]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
cleanGSS <- function(var, type = "full") {
  data <- gss.long %>%
    select(idnum, ds, minage, wtpannr123, wave, var) %>%
    mutate(wave = paste("y", wave, sep = "")) %>%
    spread(wave, var) %>%
    na.omit() %>%
    mutate(var = var) %>%
    mutate(change = ifelse(y1 != y2, 1, 0)) 
  return(data)
}



library(lavaan)

sdm <- '
y1 ~ 0
y2 ~ 0
y3 ~ 0
Mu ~ 1

Mu =~ 1*y1 + 1*y2 + 1*y3
'

aum <- '
y1 ~ 0
y2 ~ 0
y3 ~ 0
Mu ~ 1

Mu =~ 1*y1 + 1*y2 + 1*y3

y1 ~~ cov*y2
y2 ~~ cov*y3
'

sem.results <- vector(mode = "list", length = length(gss_vars))
for (i in 1:length(gss_vars)) {
  var <- gss_vars[i]
  
  df <- cleanGSS(var)
  
  sdm.fit <- sem(sdm, df, orthogonal = TRUE, missing = "ML",
                 warn = FALSE, sampling.weights = "wtpannr123")
  sdm.bic <- BIC(sdm.fit)
  aum.fit <- sem(aum, df, orthogonal = TRUE, missing = "ML",
                 warn = FALSE, sampling.weights = "wtpannr123")
  aum.bic <- BIC(aum.fit)
  
  sem.results[[i]] <- data.frame(var = var, sdm.bic = sdm.bic, aum.bic = aum.bic)
  print(paste(i, var))

}


bind_rows(sem.results) %>%
  mutate(bic.diff = aum.bic - sdm.bic,
         or.aum = exp(-bic.diff/2),
         prob.aum = or.aum/(or.aum + 1)) %>%
  arrange(desc(prob.aum))


both.results <- bind_rows(sem.results) %>%
  mutate(bic.diff = aum.bic - sdm.bic,
         or.aum = exp(-bic.diff/2),
         prob_aum_sem = or.aum/(or.aum + 1)) %>%
  arrange(desc(prob_aum_sem)) %>%
  left_join(res.df)

both.results %>%
  ggplot(aes(x = prob_aum_sem)) + 
  geom_histogram(color = "black", fill = "gray", bins = 20) + 
  theme_classic() + 
  labs(x = "Probability AUM", y = "Count",
       title = "Probability of Active Updating for GSS Questions",
       subtitle = "Probabilities generated using a structural equation models")

disagrees <- both.results %>%
  mutate(disagree = ifelse(prob_aum_nls > .5 & prob_aum_sem < .5, 1, 0)) %>%
  filter(disagree == 1)



both.results %>%
  filter(prob_aum_nls > .5, prob_aum_sem < .5) %>%
  ggplot(aes(x = prob_aum_nls, y = prob_aum_sem)) + 
  geom_point(shape = 21, fill = "gray") +
  geom_hline(yintercept = .5, linetype = 2, color = "gray") + 
  geom_vline(xintercept = .5, linetype = 2, color = "gray") + 
  geom_text_repel(data = disagrees, aes(label = var), size = 3) + 
  theme_bw() + 
  labs(x = "Probability AUM under regression",
       y = "Probability AUM under SEM",
       title = "Comparison of likelihood of AUM under different models")



both.results %>%
  filter(prob_aum_nls < .50, prob.aum > .50)

both.results %>%
  filter(prob_aum_nls > .50, prob.aum < .50)
  
  ggplot(aes(x = prob.aum)) + 
  geom_histogram()



cp <- cleanGSS("grass") %>%
  mutate(pattern = paste(w1, w2, w3, sep = "-"))

cp %>%
  group_by(pattern) %>%
  summarise(n = n())

reversions: 48 + 63
change: 86 + 63


k <- both.results %>%
  mutate(agree = ifelse(prob_aum_sem > .5 & prob_aum_nls > .5, 0,
                        ifelse(prob_aum_sem < .5 & prob_aum_nls < .5, 0, 1))) %>%
  filter(agree == 1)

