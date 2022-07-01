#Slope test functions

slopetest <- function(v1, v2, v3, v4) {
  slope1 <- as.numeric(lm(v1 ~ v2)$coefficients[2][1])
  slope2 <- as.numeric(lm(v3 ~ v4)$coefficients[2][1])
  ym1 <- mean(v1)
  ym2 <- mean(v3)
  xm1 <- mean(v2)
  xm2 <- mean(v4)
  ssy1 <- sum((v1 - ym1)^2)/(length(v1)-2)
  ssy2 <- sum((v3 - ym2)^2)/(length(v3)-2)
  c1 <- var(v2)*(length(v2)-1)
  c2 <- var(v4)*(length(v4)-1)
  z <- (slope1 - slope2)/(((ssy1*ssy2)/(c1*c2))^(0.5))
  return(z)
}

slopetest2 <- function(v1, v2, v3, v4) {
  slope1 <- as.numeric(lm(v1 ~ v2)$coefficients[2][1])
  slope2 <- as.numeric(lm(v3 ~ v4)$coefficients[2][1])
  j1 <- summary(lm(v1 ~ v2))
  j2 <- summary(lm(v3 ~ v4))
  SE1 <- as.numeric(j1$coefficients[2,2])
  SE2 <- as.numeric(j2$coefficients[2,2])
  z <- (slope1 - slope2)/((SE1^2+SE2^2)^(0.5))
  return(z)
}


