rice <- read.csv("C:/University/Year 3 Sem 2/Time Series/rice_yield_D.csv", stringsAsFactors = FALSE)

# Convert date column to Date class
rice$date <- as.Date(rice$date, format = "%d/%m/%Y")

start_year  <- as.numeric(format(min(rice$date), "%Y"))
start_month <- as.numeric(format(min(rice$date), "%m"))
start_year
start_month

rice.ts <- ts(rice$x, start = c(start_year, start_month), frequency = 12)

#1. Time Plot & Description
plot(rice.ts, main = "Monthly Rice Yield Index", ylab = "Rice Yield Index", xlab = "Time")

boxplot(rice.ts ~ cycle(rice.ts), xlab = "Month", ylab = "Rice Yield Index",
        main = "Seasonal Boxplot of Rice Yield Index")

#2. Time Series Decomposition
decomp_add <- decompose(rice.ts, type = "additive")
plot(decomp_add)

decomp_mul <- decompose(rice.ts, type = "multiplicative")
plot(decomp_mul)

acf(decomp_add$random[7:138])
acf(decomp_mul$random[7:138])



#3. Correlation Structure
#acf(rice.ts, main = "ACF of Rice Yield Index")
acf(as.numeric(rice.ts), lag.max = 36, main = "ACF of Rice Yield Index")

#4. Time Series Regression Model
#An additive regression model with a linear trend and monthly seasonal indicators 
#(OLS Method)
Time   <- time(rice.ts)
Season <- cycle(rice.ts)

rice.lm <- lm(rice.ts ~ 0 + Time + factor(Season))

acf(resid(rice.lm))
pacf(resid(rice.lm))

#(GLS Method)
library(nlme)
rice.gls <- gls(rice.ts ~ 0 + Time + factor(Season), cor = corARMA(p = 4))

#Harmonic Seasonal Model
SIN <- COS <- matrix(nr = length(rice.ts), nc = 6)
for (i in 1:6){
  COS[, i] <- cos(2 * pi * i * time(rice.ts))
  SIN[, i] <- sin(2 * pi * i * time(rice.ts))
}
meantime <- mean(time(rice.ts))
sdtime <- sd(time(rice.ts))
TIME <- (time(rice.ts) - meantime) / sdtime

rice.harmonic.gls <- gls(rice.ts ~ TIME + COS[, 1] + SIN[, 1]
                   + COS[, 2] + SIN[, 2] + COS[, 3] + SIN[, 3]
                   + COS[, 4] + SIN[, 4] + COS[, 5] + SIN[, 5]
                   + COS[, 6] + SIN[, 6], cor = corARMA(p = 4), method = "ML")
coef(rice.harmonic.gls)/sqrt(diag(vcov(rice.harmonic.gls)))

t_ratio <- coef(rice.harmonic.gls)/sqrt(diag(vcov(rice.harmonic.gls)))
t_ratio[t_ratio > 2 | t_ratio < -2]


#Harmonic Seasonal Model with only retain the significant sine and / or cosine terms
rice.harmonic.gls.reduced <- gls(rice.ts ~ TIME + COS[, 1] + SIN[, 1] + COS[, 2] + SIN[, 2] + COS[, 5], 
                   cor = corARMA(p = 4), method = "ML")
summary(rice.harmonic.gls.reduced)


#Compare 4 Regression Model
AIC(rice.lm)
AIC(rice.gls)
AIC(rice.harmonic.gls)
AIC(rice.harmonic.gls.reduced)


#Detail Analysis of The best model
summary(rice.gls)

z.gls <- resid(rice.gls)
acf(z.gls)
pacf(z.gls)



z.gls.ar <- ar(z.gls, method = "mle")
z.gls.ar
z.gls.ar$ar
mean(z.gls)

w.gls <- z.gls.ar$resid
acf(w.gls[-(1:4)])
qqnorm(w.gls[-(1:4)])
qqline(w.gls[-(1:4)])


# Fitted values from GLS model
fitted.gls <- fitted(rice.gls)

fitted.gls.ts <- ts(fitted.gls,
                    start = start(rice.ts),
                    frequency = frequency(rice.ts))
plot(rice.ts,
     main = "Observed and Fitted Rice Yield Index (GLS Model)",
     ylab = "Rice Yield Index",
     xlab = "Time")

lines(fitted.gls.ts, col = "red", lwd = 2)

legend("topleft",
       legend = c("Observed", "Fitted (GLS)"),
       col = c("black", "red"),
       lty = c(1,1),
       lwd = c(1,2))

#AR Model
rice.arima <- arima(rice.ts, order = c(1,1,1), method = "CSS")
rice


#SARIMA 
get.best.sarima <- function(rice.ts, maxord = c(1,1,1,1,1,1))
























