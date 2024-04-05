
#Week_1

#H_W1_1

v1 <- c(1, 5, 12, 10)
v2 <- rep(c(0,1),10)
v3 <- seq(from = 0, to = 20, by = 0.5)
v4 <- 10:30
v4_1 <- v4*rev(v2)
v4_2 <- v4_1+v2
v4_2 <- v4_1^rev(v2)

#H_W1_2

v1%%2
v4%%2
v3%%2
which(v1%%2==1)

#H_W1_3
v5 <- c(v1,sum(v1))
v5

#H_W1_4

vect_4 <-1:10000000%%2&1:10000000%%3&1:10000000%%5&1:10000000%%7
sum(vect_4)

10000000-sum(which(1:10000000%%210==0)^0)

length(which(1:10000000%%2 !=0&
               1:10000000%%3 !=0&
               1:10000000%%5 !=0&
               1:10000000%%7  !=0))

length(which(1:10000000%%210 !=0
              ))

#H_W1_5

matr_1 <- matrix(1:100,nrow =100,ncol = 100, byrow = T )
matr_2 <- t(matr_1)
matr <- matr_2/matr_1
matr

vv1 <- rep(1:100,times=100)
vv2 <- rep(1:100,each=100)
matrix(vv1/vv2,100,100)


#H_W1_6

v6 <- v3*v4
v6
sum(v6<100 | v6>300)
sum(v6>=200 & v6<=400)
sum(v6[-(11:(length(v6)-10))])


#H_W1_7

v3
floor(length(v3)*0.1)
floor(length(v3)*0.9)
v7 <- v3[-(floor(length(v3)*0.1):floor(length(v3)*0.9))]

#H_W1_8

sort(v1)[3]

#H_W1_9

round(pi, digit=15)

v8 <- c(1:1000)
sqrt(sum(1/v8^2)*6)

v9 <- pi

floor(pi*10^15)%%10


#Week_2

#H_W2_1

Vec <- c(1:8,200,300)


Vec <- c(Vec, round(rnorm(n =10,mean = mean(Vec),
               sd = var(Vec)),digits = 2) )
Vec


#H_W2_2

rpois(n = 1,lambda = 1000)

Vec1 <- runif(n = rpois(n = 1,lambda = 1000),
      min = 1,max =100 )

length(Vec1)

max(Vec1) - min(Vec1)
  
sd(Vec1)/ mean(Vec1)
  

3*(mean(Vec1)-median(Vec1))/sd(mean(Vec1))


(quantile(Vec1, probs = 0.1)+ quantile(Vec1, probs = 0.9) - 2* median(Vec1))/
  (quantile(Vec1, probs = 0.9) - quantile(Vec1, probs = 0.1))

sum(quantile(Vec1, probs = c(0.1,0.9)))


#H_W2_3

letters
LETTERS


sample(x = letters,size = 100,replace = T)
letters[sample(x =1:26,size =100,replace = T)]

sum(sort(unique(sample(x = letters,size = 100,replace = T)))==letters)
sum(sort(unique(letters[sample(x =1:26,size =100,replace = T)]))==letters)

all(letters %in% sample(x = letters,size = 100,replace = T))
length(unique(sample(x = letters,size = 100,replace = T)))==length(letters)



paste0(sample(x = letters,size = floor(runif(n = 1,min =1,max = 10 )) ,
              replace = T), collapse = "")


paste(
  
paste0(sample(x = LETTERS,size = 2 , replace = T), collapse = ""),
paste0(sample(x = 0:9,size = 3 , replace = T), collapse = ""),
paste0(sample(x = LETTERS,size = 2 , replace = T), collapse = "")

)


#H_W2_4

(rnorm(n = 10000,mean = 100,sd = 20) - 100)/20

mean((rnorm(n = 10000,mean = 100,sd = 20) - 100)/20)
sd((rnorm(n = 10000,mean = 100,sd = 20) - 100)/20)


#H_W2_5

rnorm(n = 2, mean = c(1, 100), sd = c(0.01, 30))

Vec2 <- rpois(n = 10000,lambda = runif(n = 10000,min = 0,max = 10))

mean(Vec2)
var(Vec2)

mean(
  sample(Vec2, size=1000, replace =T )
)

var(
  sample(Vec2, size=1000, replace =T)
)

Vec_lambda <- runif(n = 10000,min = 0,max = 10)
mean(Vec_lambda)
mean(Vec2)


#H_W2_6

options(digits = 9)

Vec3 <- rpois(n = 9,lambda = 1:9)

(Vec3[1]*10^8+Vec3[2]*10^7+Vec3[3]*10^6+Vec3[4]*10^5+Vec3[5]*10^4+
  Vec3[6]*10^3+Vec3[7]*10^2+Vec3[8]*10+Vec3[9])/10^9

as.numeric(paste0(rpois(n = 9,lambda = 1:9), collapse = ""))/10^9

#H_W2_7

p1 <- 2^1*exp(-2)/factorial(1)
p2 <- 2^2*exp(-2)/factorial(2)
p3 <- 2^3*exp(-2)/factorial(3)


#H_W2_8

Vec4 <- rnorm(n = 10^6, mean = 150,sd = 20)
sum(abs(Vec4-150)<60)/10^6*100

sum(abs(Vec4-150)<3*sd(Vec4))/length(Vec4)*100


#H_W2_9

aa)

V_pois <- rpois(n = 10^6,lambda = 2)


length(which(V_pois>=4))/10^6
mean(V_pois)/4

a <- 4

sum(V_pois>=a)/length(V_pois)
mean(V_pois)/a



length(which(rpois(n = 10^6,lambda = 2)>=3))/10^6
mean(rpois(n = 10^6,lambda = 2))/3

length(which(rpois(n = 10^6,lambda = 2)>=2.1))/10^6
mean(rpois(n = 10^6,lambda = 2))/2.1


V_unif <- runif(n = 10^6,min = 1,max = 100)
length(which(V_unif>=60))/10^6
mean(V_unif)/60

length(which(runif(n = 10^6,min = 1,max = 100)>=90))/10^6
mean(runif(n = 10^6,min = 1,max = 100))/90



V_norm <- rnorm(n = 10^6,mean = 10,sd = 4)
length(which(V_norm>=11))/10^6
mean(V_norm)/11

length(which(rnorm(n = 10^6,mean = 10,sd = 4)>=100))/10^6
mean(rnorm(n = 10^6,mean = 10,sd = 4))/100

length(which(rnorm(n = 10^6,mean = 10,sd = 4)>=10.1))/10^6
mean(rnorm(n = 10^6,mean = 10,sd = 4))/10.1


bb)

length(which(V_pois-2>=3))/10^6
var(V_pois)/3^2

length(which(V_pois-2>=1.5))/10^6
var(V_pois)/1.5^2

length(which(V_pois-2>=1.5))/10^6<= var(V_pois)/1.5^2



X <- rnorm(n = 10^1,mean = 10,sd = 10)
h_X <- X^2

mean(X)^2
mean(h_X)

mean(X)^2<=mean(h_X)



# Week_3

#H_W3_1

H3V1 <- rnorm(n = 10000,mean = 10,sd = 10)

H3V2 <- cummax(H3V1)-cummin(H3V1)

min(which(H3V2==H3V2[10000]))


#H_W3_2

H3vec <- rnorm(n = 10000,mean = 100,sd = 20)
H3vec_w <- rpois(n = 10000,lambda = 100)
weighted.mean(H3vec, H3vec_w)

mean(rep(x = H3vec, times= H3vec_w))

#/rep - ფუნქციაში, როდესაც არგუმენტებად მოცემულია ვექტორები, მაშ პირველიინ ვექტორის 
#პირველ ელემეტს აიღებსმეორე  ვექტორის  პირველ ელემენტჯერ, მერეს მეორე ვექტორ
#ის მეორე ელემენჯერ და ა.შ. მათი საშუალო არის იგივე შეწონილი საშუალო.

weighted.mean(H3vec,sqrt(H3vec_w))

H3vec <- c(rnorm(n = 10000,mean = 100,sd = 20),10^6)
H3vec_w <- c(rpois(n = 10000,lambda = 100),10^6)
weighted.mean(H3vec, H3vec_w)
weighted.mean(H3vec,sqrt(H3vec_w))


#H_W3_3

rexp(n = 10^4,rate = 0.1)
hist(rexp(n = 10^4,rate = 0.1))

EX=1/0.1=10
DX=1/0.1^2=100
mean(rexp(n = 10^4,rate = 0.1))
var(rexp(n = 10^4,rate = 0.1))
  

rexp(n = 10^4,rate = 10)
hist(rexp(n = 10^4,rate = 10))

EX=1/10=0.1
DX=1/10^2=0.01
mean(rexp(n = 10^4,rate = 10))
var(rexp(n = 10^4,rate = 10))


round(quantile(x =rexp(n = 10^4,rate = 10),probs = seq(0,1,0.1)),4)

round(-log(1-0.2)/10,4)
round(-log(1-0.5)/10,4)
round(-log(1-0.9)/10,4)

min(rexp(n = 10^4,rate = 10))
max(rexp(n = 10^4,rate = 10))

#/uniroot(f = func, interval = c(0,5))$root

H3V3 <- rexp(n = 10^6,rate = 1/10)
mean(H3V3)
var(H3V3)
round(quantile(H3V3, probs = seq(0,1,0.1)), 3)

-log(1-0.8)*10

func <- function(x){-log(1-x)*10-max(H3V3)}
uniroot(f = func,interval = c(0,1))$root

min(H3V3)
max(H3V3)





#H_W3_4

H3V6 <- rweibull(n = 10^6,shape = 2,scale = 30)
EX <- 30*gamma(1+1/2)
medX <- 30*log(2)^(1/2)
mean(H3V6)
median(H3V6)

func1 <- function(x){gamma(1+1/x)/log(2)^(1/x)-mean(H3V6)/median(H3V6)}
k <- uniroot(f = func1,interval = c(0,100))$root

func2 <- function(x){x*log(2)^(1/k)-median(H3V6)}
lambda <- uniroot(f = func2,interval = c(0,100))$root




#H_W3_5

cut(x = c(1,2,5,10),
    breaks = c(0, 4, Inf),
    right = F)

cut(x = c(1,2,5,10),
    breaks = 3,
    right = F)

library(epiDisplay)

H3V4 <- rnorm(n = 10000,mean = 100,sd = 1000)

H3V5 <- cut(x = H3V4,breaks = 10,right = T)
tab1(H3V5)

H3V6 <- cut(x = H3V4,breaks = 10,right = T, labels = 1:10)

tab1(H3V6)





# Week_4

library(epiDisplay)
library(psych)


#H_W4_1

iris_tr <- iris
colnames(iris_tr)
sum(iris_tr$Petal.Length>=median(iris_tr$Petal.Length))

#H_W4_2

iris_tr$Prop_Petal.Width <- prop.table(iris_tr$Petal.Width)*100
iris_tr
iris_tr[iris_tr$Prop_Petal.Width>1,]

which(iris_tr$Prop_Petal.Width>1)


#H_W4_3

iris_tr$SP.Area <- (iris_tr$Sepal.Length*iris_tr$Sepal.Width+
                   iris_tr$Petal.Length*iris_tr$Petal.Width)*3

describeBy(x = iris_tr$SP.Area,group = iris_tr$Species)

sum(iris_tr[grepl("setosa",iris_tr$Species) ,"SP.Area"])
sum(iris_tr[iris_tr$Species=="setosa" ,"SP.Area"])


sum(iris_tr[grepl("versicolor",iris_tr$Species) , "SP.Area"])
sum(iris_tr[grepl("virginica",iris_tr$Species) ,"SP.Area"])

sum(iris_tr[grepl("setosa",iris_tr$Species) ,"SP.Area"])/
  min(iris_tr[grepl("setosa",iris_tr$Species) ,"SP.Area"])
sum(iris_tr[grepl("versicolor",iris_tr$Species) , "SP.Area"])/
  min(iris_tr[grepl("versicolor",iris_tr$Species) , "SP.Area"])
sum(iris_tr[grepl("virginica",iris_tr$Species) ,"SP.Area"])/
  min(iris_tr[grepl("virginica",iris_tr$Species) ,"SP.Area"])

sum(iris_tr$SP.Area)/min(iris_tr$SP.Area)

quantile(iris_tr$SP.Area,probs = 0.9)
quantile(iris_tr$SP.Area,probs = 0.1)

iris_tr_Quant <- iris_tr[(iris_tr$SP.Area>quantile(iris_tr$SP.Area,probs = 0.1))&
          (iris_tr$SP.Area<quantile(iris_tr$SP.Area,probs = 0.9)),]

sum(iris_tr_Quant$SP.Area)/min(iris_tr_Quant$SP.Area)

mean(iris_tr_Quant$SP.Area)

mean(iris_tr$SP.Area, trim = 0.1)



#H_W4_4

iris_tr2 <- iris
which((iris_tr2$Sepal.Length*10)%%10==0|(iris_tr2$Sepal.Width*10)%%10==0
         ,(iris_tr2$Petal.Length*10)%%10==0 |(iris_tr2$Petal.Width*10)%%10==0)

ceiling(15.1)
floor(14.1)

#H_W4_5

levels(iris_tr2$Species) <- c(levels(iris_tr2$Species), "variegate")

iris_tr2$Species[141:150] <- "variegate"
iris_tr2$Species[(length(iris_tr2$Species)-9):length(iris_tr2$Species)] <- "variegate"

#H_W4_6

iris_tr$length_av <- (iris_tr$Sepal.Length+iris_tr$Petal.Length)/2

1:2+3:4
sum(1:2,3:4)
median(1:2,3:4)
median(c(1:2,3:4))

#H_W4_7

describe(iris_tr$Sepal.Length)

Sample <- sample(x = iris_tr$Sepal.Length,size = 100, replace = T)
sd(Sample)/sqrt(100)

#H_W4_8

iris_tr$characterization <- ifelse(test = iris_tr$Sepal.Length>5,1,
                                   ifelse(test = iris_tr$Sepal.Width<3,2,
                                          ifelse(iris_tr$Petal.Length>4&iris_tr$Petal.Length<5,3,
                                          ifelse(iris_tr$Petal.Width%%0.2==0,4,5))))
  
table(iris_tr$characterization) 


# Week_5  

library(palmerpenguins)
library(lubridate)
  
#H_W5_1

penguins_tr <- penguins

str(penguins_tr)
summary(penguins_tr)

#H_W5_2

penguins_tr[is.na(penguins_tr$bill_length_mm) | is.na(penguins_tr$bill_depth_mm) | 
              is.na(penguins_tr$flipper_length_mm) |
              is.na(penguins_tr$body_mass_g) | is.na(penguins_tr$sex) , ]


penguins_tr$NN <- apply(X = is.na(penguins_tr)[,c(1:8)], MARGIN = 1, FUN = sum)
penguins_tr[penguins_tr$NN>0, ]


penguins_tr[rowSums( is.na(penguins_tr)[,c(1:8)])>0 ,]


#H_W5_3

aggregate(x = body_mass_g~year, data = penguins_tr,FUN = mean)

aggregate(x = body_mass_g~year, data = penguins_tr,FUN = mean,na.action = na.pass)

aggregate(x = body_mass_g~sex, data = penguins_tr,FUN = mean)

aggregate(x = body_mass_g~sex, data = penguins_tr,FUN = mean,na.action = na.pass)


penguins_tr$sex_with_NA <- penguins_tr$sex

levels(penguins_tr$sex_with_NA) <- c(levels(penguins_tr$sex_with_NA), "NA")

penguins_tr$sex_with_NA[is.na(penguins_tr$sex_with_NA)] <- "NA"

sum(is.na(penguins_tr$sex_with_NA))

aggregate(x = body_mass_g~penguins_tr$sex_with_NA, data = penguins_tr,FUN = mean)



#H_W5_4

f <- function(x){c(cv=sd(x)/mean(x), range=max(x)-min(x), 
                   IQT = quantile(x, 0.75, name=F) - quantile(x, 0.25, name=F))}

aggregate(x = cbind(penguins_tr$bill_length_mm, penguins_tr$flipper_length_mm)~species , 
          data = penguins_tr, FUN = f)


#H_W5_5

penguins_tr$rtbill <- penguins_tr$bill_length_mm/penguins_tr$bill_depth_mm

mean(penguins_tr$rtbill, na.rm = T)
weighted.mean(penguins_tr$rtbill,w = penguins_tr$bill_depth_mm, na.rm=T)
sum(penguins_tr$bill_length_mm, na.rm = T)/sum(penguins_tr$bill_depth_mm, na.rm=T)

length_per_depth <- sum(penguins_tr$bill_length_mm, na.rm = T)/sum(penguins_tr$bill_depth_mm, na.rm=T)



#H_W5_6

penguins_tr$bill_length_mm_adj <- penguins_tr$bill_depth_mm*length_per_depth
penguins_tr$bill_length_mm_max <- max(penguins_tr$bill_length_mm, penguins_tr$bill_length_mm_adj  , na.rm=T)

length(which(penguins_tr$bill_length_mm_max==penguins_tr$bill_length_mm))
sum(penguins_tr$bill_length_mm_max==penguins_tr$bill_length_mm, na.rm = T)

#H_W5_7

c(split(penguins_tr, penguins_tr$sex),
  split(penguins_tr, penguins_tr$island))



list(split(penguins_tr, penguins_tr$sex),
     split(penguins_tr, penguins_tr$island))


#H_W5_8

lst <- list(3, c(rep(1, 100), rep(2,50)) , matrix(data = 1:100, nrow = 10, ncol = 10),
            penguins_tr[penguins_tr$island=="Torgersen", c(2,3,4)], rep(NA, 10) )

numbers <- as.data.frame(lst[3])

mean(numbers)

mean(matrix(data = 1:100, nrow = 10, ncol = 10))



#H_W5_9

penguins_raw_tr <- penguins_raw

penguins_tr[, c(5,6,8)]

penguins_tr$ID <- penguins_tr$flipper_length_mm*penguins_tr$body_mass_g*penguins_tr$year
penguins_tr$ID[which(is.na(penguins_tr)[,"ID"]==T)] <- 0

length(penguins_tr$ID)
length(unique(penguins_tr$ID))


penguins_raw_tr$year <- year(penguins_raw_tr$`Date Egg`)

penguins_raw_tr$ID <- penguins_raw_tr$`Flipper Length (mm)`*penguins_raw_tr$`Body Mass (g)`*penguins_raw_tr$year
penguins_raw_tr$ID[which(is.na(penguins_raw_tr)[,"ID"]==T)] <- 0

length(penguins_raw_tr$ID)
length(unique(penguins_raw_tr$ID))


merge(x = penguins_tr ,y = penguins_raw_tr, by = "ID")

penguins_tr1 <- penguins_tr[penguins_tr$ID>0, ]
penguins_raw_tr1 <- penguins_raw_tr[penguins_raw_tr$ID>0, ]


length(penguins_tr1$ID)
length(unique(penguins_tr$ID))

length(penguins_raw_tr1$ID)
length(unique(penguins_raw_tr1$ID))

merge(x = penguins_tr1 ,y = penguins_raw_tr1, by = "ID")


# Week_6  

library(dplyr)
library(tidyr)

iris_tr <- iris

#H_W6_1

iris_tr %>% mutate(SL_bn=case_when(Sepal.Length<5~1, Sepal.Length>6~3, TRUE~2))

#H_W6_2

f <- function(x){c(mn=mean(x), md=median(x), sd=sd(x), cv=mean(x)/sd(x),d9=quantile(x,0.9, name=F))}
iris_tr%>% group_by(Species) %>% summarise_if(.predicate = is.numeric,.funs = f)

#H_W6_3

table(iris_tr %>% mutate(PL=case_when(Petal.Length>Sepal.Length~1, TRUE~0)) %>% select(PL))
iris_tr %>% mutate(PL=case_when(Petal.Length>Sepal.Length~1, TRUE~0)) %>% select(PL) %>% table()

table(iris_tr %>% rowwise() %>% mutate(var3=mean(c(Sepal.Length,Petal.Length))) %>% 
  mutate(SPmn=case_when(var3<4~1, TRUE~0)) %>% select(SPmn))

iris_tr %>% rowwise() %>% mutate(var3=mean(c(Sepal.Length,Petal.Length))) %>% 
  mutate(SPmn=case_when(var3<4~1, TRUE~0)) %>% select(SPmn) %>% table()

iris_tr %>% arrange(Sepal.Length) %>% slice(length(Sepal.Length)-2) %>% select(Sepal.Length)
iris_tr %>% arrange(desc(Sepal.Length)) %>% slice(3) %>% select(Sepal.Length) %>% as.numeric()
iris_tr %>% arrange(desc(Sepal.Length)) %>% slice(3) %>% select(Sepal.Length) %>% c() %>% unlist()

(iris_tr %>% select(Sepal.Length) %>% arrange(Sepal.Length))[length(iris_tr$Sepal.Length)-2,]


iris_tr %>% group_by(Species) %>% 
  summarise(mn=mean(Sepal.Length), sd=sd(Sepal.Length), 
            Low=mean(Sepal.Length)-1.96*sd(Sepal.Length),
            UP=mean(Sepal.Length)+1.96*sd(Sepal.Length)) %>% 
            inner_join(y = iris_tr,by = "Species")


#H_W6_4

iris_tr %>% pivot_longer(cols = Sepal.Length:Petal.Width, 
                         names_to = "SP", values_to ="LW" ) %>% distinct(LW) %>% nrow()

iris_tr %>%distinct(Sepal.Length,Sepal.Width,Petal.Length,Petal.Width) %>% nrow()

#H_W6_5

iris_tr%>% group_by(Species) %>% summarise(SLmean=mean(Sepal.Length)) %>% 
           inner_join(y = iris_tr,by = "Species")
  


