
library(aplpack)

C1 <- read.csv("/Users/Lorenzo/Desktop/Tesi/C1.csv")

C2 <- read.csv("/Users/Lorenzo/Desktop/Tesi/C2.csv")

C3 <- read.csv("/Users/Lorenzo/Desktop/Tesi/C3.csv")

C4 <- read.csv("/Users/Lorenzo/Desktop/Tesi/C4.csv")

C5 <- read.csv("/Users/Lorenzo/Desktop/Tesi/C5.csv")

C6 <- read.csv("/Users/Lorenzo/Desktop/Tesi/C6.csv")

clients = 

# C1

C1=C1[,-1]

C1$Gender=as.factor(C1$Gender)
C1$Job=as.factor(C1$Job)
C1$Area=as.factor(C1$Area)
C1$CitySize=as.factor(C1$CitySize)
C1$Investments=as.factor(C1$Investments)

hist(C1$LifeStyle)

