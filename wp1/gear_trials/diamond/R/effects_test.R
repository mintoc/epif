library(effects)
data(Arrests)

Arrests$year <- as.factor(Arrests$year)
arrests.mod <- glm(released ~ employed + citizen + checks +
                   colour*year + colour*age,
                   family=binomial, data=Arrests)

summary(arrests.mod)

