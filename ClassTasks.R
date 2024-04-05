*Day1

2+2
5*8
2^7
8/4
numb <- 5
numb <- 6
numb <- 2*numb
numb^numb
rm(numb)

vect1 <- c(1,2,5,10)
1:100
vect2 <- 0.4:100.1
0.4:100.1
vect2
seq(0, 10, 0.3)
seq(0,1,0.1)
seq(0,1, length.out =13)

seq(from)
round(digits =0, 2.32654 )
help (seq)
?rep
rep(1,10)
rep(1,times=10)
rep(1,length.out=10)
rep(c(1,2),length.out=10)
rep(c(1,2),each=10)
vec <- c(rep(4,5), seq(0,10,2.3), 4:1, pi,3)
vec
vec <- c(rep(4,5), seq(0,10,2.3), 4:1, round(digits =2, pi ),3)
vec

vect1+vec

vect1^vec

length(vec)
rev(vec)
unique(vect1)
sort()
?sort

sort(vec, T)
sort(vec, F)
sort(vec)
vect2
vect1[c(1,4)]
subvect <- vect1[c(1,4)]
vect1[c(T,F,F,T)]
vect1>=5
vect1[vect1>=5]
vect1[3:4]
which(vect1==5)
which(vect1>=5)

vect1[vect1>5 | vect1==1]

vect1[vect1>5 & vect1==1]
which(vect1>5 | vect1==1)

5>6
6==6

c(5>6,
6==6)


matrix(data = 1:4, nrow = 2,ncol = 2,byrow = F)
matr <- matrix(data = 1:4, nrow = 2,ncol = 2,byrow = F)
matrix(1:25, 5,5)
matr1 <- matrix(1:25, 5,5)
matr1[2:4,1:2]
array(1:25,dim = c(5,5))
array(1:125,dim = c(5,5,5))
ar <- array(1:125,dim = c(5,5,5))


dim(matr1)
dim(ar)
nrow(matr1)
dim(matr1)[1]


1:5 == c(1,2,3,4,10)

all(1:5 == c(1,2,3,4,10))
all(c(T,F,F))
any(c(T,F,F))

which((1:5 == c(1,2,3,4,10))==F)

which((1:5 == c(1,2,3,40,10))==F)
sort(matr1)
matr2 <- matrix(25:1,5,5)
matr1

cbind(matr1,matr2)

rbind(matr1,matr2)

vv <- c(1:1000)
sqrt(sum(1/vv^2)*6)

*Day2

mean(c(1,2,3))
min(c(1,2,3))
median(c(1,2,3))

var(c(1,2,3,100))

sd(c(1,2,3,100))


vv1 <- 1:10
vv2 <- c(rep(5,7),1,2,3)

cov(vv1,vv2)
cor(vv1,vv2)
cor(vv1,vv2,method = "kendal")
cor(vv1,vv2,method = "spearman")
vv3 <- seq(-1,1,0.01)
cor(vv3,vv3^2)

vv4 <- seq(0,2*pi,0.01)
cor(vv4,vv4^2)

vv5 <- seq(0,2*pi,0.01)
cor(vv5,sin(vv5))

cor(vv1, 5*vv1+150)
cor(vv1, -5*vv1+150)

cov(vv1, -5*vv1+150)/(sd(vv1)*sd(-5*vv1+150))

vv6 <- c(1,5,10,12)
mean((vv6-mean(vv6))^2)
var(vv6)
sum((vv6-mean(vv6))^2)/(length(vv6)-1)

v_values <- c(1,5,20)
v_probs <- c(0.1,0.3,0.6)

sample(x=v_values, size=10, replace =T, prob =v_probs )

hist(sample(x=v_values, size=1000, replace =T, prob =v_probs ))


sample(x=1:1000000, size=10, replace =F)

EX <- sum(v_values*v_probs)
DX <- sum((v_values-EX)^2*v_probs)

mean(
  sample(v_values, size=1000000, replace =T, prob =v_probs )
)

var(
  sample(v_values, size=1000000, replace =T, prob =v_probs )
)

rpois(n = 1000,lambda = 3)

hist(rpois(n = 1000,lambda = 3))
mean(rpois(n = 10000,lambda = 4))
var(rpois(n = 10000,lambda = 4))

runif(n = 1000,min = 1,max =17 )
hist(runif(n = 1000,min = 1,max =17 ))
(1+17)/2
(17-1)^2/12
mean(runif(n = 1000,min = 1,max =17 ))
var(runif(n = 1000,min = 1,max =17 ))

rnorm(n =1000,mean = 10,sd = 5 )
hist(rnorm(n =1000,mean = 10,sd = 5 ))
mean(rnorm(n =1000,mean = 10,sd = 5 ))
var(rnorm(n =1000,mean = 10,sd = 5 ))
sd(rnorm(n =1000,mean = 10,sd = 5 ))


*Day3

vec1 <- sample(1:4,size = 1000, replace =T, prob = c(0.5,0.2,0.2,0.1) )


substr(x = "statistics", start = 2, stop = 5)

substr(x = "statistics", start = c(2,3), stop = c(5,4))

substring(text = "statistics", first =c(2,3),last = c(5,4) )

vec2 <-substring(paste0(letters[1:3][sample(3,10000,T)],collapse = ""),first = seq(1,10000,2),
                 last = seq(2,10000,2))

t1 <- table(vec1)
table(vec2)
View(t1)

class(t1)

t2 <- prop.table(t1)

t3 <- cumsum(prop.table(t1))

freq_table <- data.frame(varianti=sort(unique(vec1)), freq=as.data.frame(t1)[,2],
           prop_freq=as.data.frame(t2)[,2],
           cum_prop=as.data.frame(t3)[,1]
           )
str(freq_table)


####*Day4####

vec1 <- sample(1:4,size = 1000, replace =T, prob = c(0.5,0.2,0.2,0.1) )


vec2 <-substring(paste0(letters[1:3][sample(3,10000,T)],collapse = ""),
                 first = seq(1,10000,2),
                 last = seq(2,10000,2))

table(vec1)
prop.table(vec1)
cumsum(vec1)

prop.table(table(vec1))
cumsum(prop.table(table(vec1)))


weighted.mean(x = sort(unique(vec1)), table(vec1))
mean(vec1)

weighted.mean(x = sort(unique(vec1)), prop.table(table(vec1)))

weighted.mean(x = sort(unique(vec1)), prop.table(table(vec1))*1459)

library(janitor)
install.packages("janitor")

tabyl(vec1)
library(epiDisplay)
tab1(vec1)

View(tab1)

barplot(table(vec1))
barplot(table(vec1),horiz = F)
barplot(table(vec1),main = "diagrama",xlab = "variantebi",ylab = "sixshire",
        col = c("red","yellow",  "blueviolet" , "violetred1"  ), ylim = c(0,600))

colours()

quantile(x =vec1,probs = seq(0,1,0.1) )

quantile(x =rnorm(10000, 100,20),probs = seq(0,1,0.1))
quantile(x =rnorm(10000, 100,20),probs = seq(0,1,0.1),type = 9)

round (quantile(x =rnorm(10000, 100,0.1),probs = seq(0,1,0.1)), 2)

prop.table(matrix(1:4,2,2))
prop.table(matrix(1:4,2,2),margin = 1)
prop.table(matrix(1:4,2,2),margin = 2)

table(vec1, vec2[1:1000])
prop.table(table(vec1, vec2[1:1000]))
prop.table(table(vec1, vec2[1:1000]),margin = 2)


dataset <- data.frame(qalaqi=c("Tbil", "Rust", "Qut"), 
           mosaxleoba=c(1400,130,150),
           dedaqalaqi=c(T,F,F))
class(dataset)           
nrow(dataset)
ncol(dataset)
dim(dataset)
names(dataset)
str(dataset)

iris
summary(iris)
iris_tr <- iris


summary(dataset)

####*Day7####

db <- data.frame(v1=c(40,50),
           v2=c("a","b"))
nrow(db)
dim(db)
names(db)

iris_tr <- iris

iris_tr$Sepal.Length
iris_tr$Sepal.Length[1:10]
sum(iris_tr$Sepal.Length[1:10])

iris_tr[,1:2]
iris_tr[,-(1:2)]
iris_tr[,c("Sepal.Length","Petal.Length")]

colnames(iris_tr)
names(iris_tr)

setdiff(1:10,5:20)
?setdiff
setdiff(5:20,1:10)


iris_tr[,!(colnames(iris_tr)%in% c("Sepal.Length","Petal.Length"))]

"%!in%" <- Negate("%in%")

iris_tr[,colnames(iris_tr)%!in% c("Sepal.Length","Petal.Length")]

iris_tr[,c(F,F,F,T,T)]

iris_tr[1:5,]
iris_tr[iris_tr$Sepal.Length>5.5,]
iris_tr[iris_tr$Sepal.Length>5.5 |iris_tr$Sepal.Width>3 ,]

iris_tr[(iris_tr$Sepal.Length>5.5 |iris_tr$Sepal.Width>3)&
          iris_tr$Species=="setosa",]


iris_tr[(iris_tr$Sepal.Length>5.5 |iris_tr$Sepal.Width>3)&
          iris_tr$Species=="setosa", 1:3]

iris_tr[(iris_tr$Sepal.Length>5.5 |iris_tr$Sepal.Width>3)&
          iris_tr$Species=="setosa", 1:3][1:5,]

vec1 <- c(1:9,NA, 10:12)

vec1

as.data.frame(vec1)

View(as.data.frame(vec1))

sum(vec1, na.rm = T)

is.na(NA)
is.na(vec1)

NA>5
is.na(NA>5)
iris_tr[!is.na(iris_tr$Sepal.Length),]

df1 <- data.frame(x=c(NA,1,2), y=c(1,NA,2))
is.na(df1)


vec1[1] <- 100

vec1[14] <- 5

vec1

vec1[19] <- 20

vec1

vec1[is.na(vec1)]

vec1[is.na(vec1)] <- mean(vec1,na.rm = T)


vec1[vec1>10] <- 10

colnames(iris_tr) <- paste0("var" , 1:5)
iris_tr[1,1] <- NA
is.na(iris_tr)
iris_tr[is.na(iris_tr)]
iris_tr[iris_tr<5]

iris_tr[is.na(iris_tr)] <- 0

data.frame((a <- rnorm(1000)), a^2)
data.frame((a <- rnorm(1000)), a^2)[1:5,]

####*Day8####

iris_tr
str(iris_tr)
str(data.frame(x=c("a","b")))

iris_tr$Species <- as.character(iris_tr$Species)
str(iris_tr)


iris_tr[1,"Species"] <- "other"

iris_tr$Species <- as.factor(iris_tr$Species)

levels(iris_tr$Species) <- c(levels(iris_tr$Species), "gulyviTela")


iris_tr <- iris

summary(iris_tr)

library(psych)
describe(iris_tr)
describe(iris_tr[,1:4], quant = c(0.1,0.9), IQR = T, trim = 0.2)

ifelse(test = 4>5,yes = 1,no = 2)

ifelse(iris_tr$Sepal.Length>5,yes = 1,no = 2)

iris_tr$new_var <- ifelse(iris_tr$Sepal.Length>5,yes = 1,no = 2)

ifelse(test = 4>5,yes = c(1,2),no = c(4,5))


if(4>5) c(1,2) else c(4,5)

if(4>5) c(1,2) else matrix(1:4,2,2)


####*Day9####

db1 <- data.frame(ID = c(1,2), y = c("a", "b"))
db2 <- data.frame(ID = c(1,3), z = c("alpha", "betta"))
db3 <- data.frame(ID1 = c(1,3), z = c("alpha", "betta"))

db4 <- data.frame(ID = c(1,1,2), code=c(1,2,1), y = c("a", "b" , "c"))
db5 <- data.frame(ID1 = c(1,1,3), code=c(1,2,1), z = c("alpha", "betta", "gamma"))


merge(x = db1,y = db2)

merge(x = db1,y = db2, by = "ID")

merge(x = db1,y = db3, by.x = "ID", by.y = "ID1")

merge(x = db1,y = db3, by.x = "ID", by.y = "ID1", all = T)

merge(x = db1,y = db3, by.x = "ID", by.y = "ID1", all.x = T)

merge(x = db1,y = db3, by.x = "ID", by.y = "ID1", all.y = T)


merge(x = db4,y = db5, by.x = "ID", by.y = "ID1")

db4 <- data.frame(ID = c(1,1,2), code=c(1,2,1), y = c("a", "b" , "c"))
db5 <- data.frame(ID = c(1,1,3), code=c(1,2,1), z = c("alpha", "betta", "gamma"))

merge(x = db4,y = db5, by = "ID")

merge(x = db4,y = db5)

merge(x = db4,y = db5, by =c("ID", "code") )

db4 <- data.frame(ID = c(1,1,2), code=c(1,2,1), y = c("a", "b" , "c"))
db5 <- data.frame(ID1 = c(1,1,3), code1=c(1,2,1), z = c("alpha", "betta", "gamma"))

merge(x = db4,y = db5, by.x = c("ID", "code"), by.y = c("ID1", "code1"))

merge(x = db4,y = db5, by.x = c("ID", "code"), by.y = c("ID1", "code1"), all = T)


####*Day10####

iris_tr <- iris

?aggregate

aggregate(x = iris_tr$Sepal.Length, by=list(iris_tr$Species), FUN=mean )

aggregate(x = iris_tr$Sepal.Width, by=list(iris_tr$Species), FUN=sum )

aggregate(x = iris_tr$Sepal.Length, by=list(iris_tr$Species), FUN=quantile )

aggregate(x = iris_tr$Sepal.Length, 
          by=list(iris_tr$Species), 
          FUN=quantile, probs = 0.6 )



aggregate(x = cbind(iris_tr$Sepal.Length, iris_tr$Sepal.Width), 
          by=list(iris_tr$Species), 
          FUN=quantile, probs = 0.6 )


aggregate(x = cbind(iris_tr$Sepal.Length, iris_tr$Sepal.Width), 
          by=list(iris_tr$Species), 
          FUN=sum)


iris_tr$size <- ifelse(iris_tr$Sepal.Length>6,1,0)

aggregate(x = cbind(iris_tr$Sepal.Length, iris_tr$Sepal.Width), 
          by=list(iris_tr$Species, iris_tr$size), 
          FUN=quantile, probs = 0.6 )



dob_f <- function(x){c(sum=sum(x),mean=mean(x),sd=sd(x))}

dob_f(iris_tr$Sepal.Length)

stat <- aggregate(x = cbind(iris_tr$Sepal.Length, iris_tr$Sepal.Width), 
          by=list(iris_tr$Species, iris_tr$size), 
          FUN=dob_f )

*/შედეგის ბაზად გარდაქმნა

stat <- as.data.frame(as.matrix(stat))


stat1 <- aggregate(x = cbind(iris_tr$Sepal.Length, iris_tr$Sepal.Width), 
                  by=list(iris_tr$Species, iris_tr$size), 
                  FUN= function(x){c(sum=sum(x),mean=mean(x, trim = 0.3),sd=sd(x))} )


aggregate(x = Sepal.Length ~ Species + size, data = iris_tr, FUN = mean)

aggregate(x = Sepal.Length ~Species+size, data = iris_tr, FUN = sum)

aggregate(x = cbind(Sepal.Length, Sepal.Width) ~ Species + size, data = iris_tr, FUN = mean)



aggregate(x = . ~ Species + size, data = iris_tr, FUN = function(x){c(sum=sum(x),
                                                                      mean=mean(x),
                                                                      sd=sd(x))})

*/striqonuli operaciebi

?apply

apply(X = iris_tr[,1:4], MARGIN = 1, FUN = mean)

apply(X = iris_tr[,1:4], MARGIN = 2, FUN = mean)

apply(X = iris_tr[,1:4], MARGIN = 2, FUN = function(x){c(sum=sum(x),
                                                         mean=mean(x),
                                                         sd=sd(x))})


lst <- list(data.frame(x=1:2, y=3:4), iris_tr[1,], 1:10, "statistika",
     NA, function(x){2*x})

class(lst)

lst[[1]]$y

nchar(lst[[4]])

lst[[6]] (16)

split(x=iris_tr, f = iris_tr$Species)
split(x = iris_tr,f = iris_tr$size)

names(split(x=iris_tr, f = iris_tr$Species))
names(split(x = iris_tr,f = iris_tr$size))

ist_irit <- split(x=iris_tr, f = iris_tr$Species)

ist_irit$setosa



library(dplyr)

iris_tr %>% filter(Sepal.Length>6.5) %>% select(Sepal.Length,Sepal.Width) %>% 
  filter(Sepal.Width<3)

iris_tr %>% filter(Sepal.Length>7) %>% select(Sepal.Length,Sepal.Width)



iris_tr %>% filter(size==0) %>% select(-c(Species)) %>% 
  summarize(sum=sum(Sepal.Length),mean=mean(Petal.Length),
            corel=cor(Sepal.Length,Petal.Length))
  
iris_tr %>% group_by(size) %>% select(-c(Species)) %>% 
  summarize(sum=sum(Sepal.Length),mean=mean(Petal.Length),
            corel=cor(Sepal.Length,Petal.Length))



iris_tr %>% filter(Sepal.Length>7) %>% select(1,2) 

iris_tr %>% filter(Sepal.Length>7) %>% select(-1) 

iris_tr %>% filter(Sepal.Length>7) %>% select(-c(Sepal.Length,Sepal.Width)) 

iris_tr %>% slice(1:10) %>% select(-c(Sepal.Length,Sepal.Width)) 

iris_tr %>% slice(1:10) %>% select_if(is.numeric)

iris_tr %>% select(Sepal.Length, Petal.Length) %>% slice(1:10)


iris_tr %>% select(is.factor) %>% slice(c(1,20,50))


#/dplyr::select()


iris_tr %>% slice(1:10) %>% select_if(is.factor)

iris_tr %>% summarize(cnt=n(),mn=mean(Sepal.Length, trim = 0.1),
                      sd=sd(Sepal.Length),
                      sd2=sd(Petal.Length),
                      covv=cor(Sepal.Length,Petal.Length))

iris_tr %>% summarise_at(.vars = c("Sepal.Length", "Petal.Length"), .funs = c( mn = mean, sd = sd))

iris_tr %>% summarise_at(.vars = c("Sepal.Length", "Sepal.Width"),  .funs = c(mn = mean, sd = sd))



iris_tr %>% group_by(Species) %>% summarize(cnt=n(), mn=mean(Sepal.Length, trim = 0.1),
                      sd=sd(Sepal.Length),
                      sd2=sd(Petal.Length),
                      covv=cor(Sepal.Length,Petal.Length))




iris_tr %>% distinct(Species)



(iris_tr %>% rowwise() %>% mutate(mnL=mean(c(Sepal.Length, Petal.Length)),
                        mnW=mean(c(Sepal.Width, Petal.Width)),
                        SP=Sepal.Length/Petal.Length) %>% 
                        select(Sepal.Length,Petal.Length,Sepal.Width,Petal.Width,
                               mnL,mnW, SP ))[1:5,]

select(.data = mutate(.data = rowwise(data = iris_tr), mnL=mean(c(Sepal.Length, Petal.Length))),
                      Sepal.Length, Petal.Length, mnL)               



library(haven)


iris_tr %>% filter(Sepal.Length>5&Petal.Length<3) %>% 
  select(Sepal.Length,Petal.Length,Species)

iris_tr %>% select(-c(1,4))

iris_tr %>% slice(1:5) %>% select(-c(4:6))

iris_tr %>% select_if(is.factor) %>% slice(c(2,100,150))
  
iris_tr %>% summarise(cnt=n(), mn=mean(Sepal.Length), sd=sd(Sepal.Length))  
  
iris_tr %>% summarise_at(.vars = c("Sepal.Length","Petal.Length"),
                         .funs = c( mn=mean, sd=sd)) 

iris_tr %>% summarise_if(.predicate = is.numeric,
                         .funs = c(mn=mean, sd=sd) )

iris_tr %>% summarise_if(.predicate = is.numeric, 
                         .funs = function(x){c(sm=sum(x), mn=mean(x), sd=sd(x))}) %>% 
  select(1:4)

iris_tr %>% group_by(size) %>% summarise_if(.predicate = is.numeric, 
                         .funs = function(x){c(sm=sum(x), mn=mean(x), sd=sd(x))})
 
iris_tr %>% mutate(av=mean(Sepal.Length)) %>% slice(1:5) %>% select(Sepal.Length, av)



iris_tr %>% rowwise() %>% mutate(sm1=sum(c(Sepal.Length,Petal.Length)), sm2=sum(c(Sepal.Width,Petal.Width)),
                                 mn1=mean(c(Sepal.Length,Petal.Length)), mn2=mean(c(Sepal.Width,Petal.Width)),
                                 rt1=Sepal.Length/Petal.Length, rt2=Sepal.Width/Petal.Width) %>% 
  select_if(is.numeric)


(iris_tr %>% rowwise() %>% mutate(sm1=sum(c(Sepal.Length,Petal.Length)), sm2=sum(c(Sepal.Width,Petal.Width)),
                                 mn1=mean(c(Sepal.Length,Petal.Length)), mn2=mean(c(Sepal.Width,Petal.Width)),
                                 rt1=Sepal.Length/Petal.Length, rt2=Sepal.Width/Petal.Width) %>% 
  select_if(is.numeric))[1:5,]

iris_tr %>% mutate(SL=Sepal.Length>5)

iris_tr %>% mutate(SL=Sepal.Length>5) %>% colnames() %>% as.data.frame()



iris_tr %>% group_by(Species) %>% summarise(SP_mean = mean(Sepal.Length),
  PL_mean = mean(Petal.Length)) 
       

iris_tr %>% group_by(Species) %>% summarise(SP_mean = mean(Sepal.Length)) %>%
  merge(y = iris_tr %>% group_by(Species) %>% summarise(PL_mean = mean(Petal.Length)), 
        by = "Species", all.x = T)


(iris_tr %>% rowwise() %>% mutate(av=mean(c_across(Sepal.Length:Petal.Width))) %>% 
  select_if(is.numeric))[1:5,]

iris_tr %>% mutate(group=case_when(Sepal.Length>5 ~ 1, Petal.Length<3 ~ 2, TRUE ~3)) %>% 
  select(Sepal.Length,Petal.Length,group) %>% slice(1:5)

iris_tr %>% distinct(Species)

iris_tr %>% arrange(Sepal.Length) %>% select(Sepal.Length) %>% slice(1:5)

iris_tr %>% count(Species)

iris_tr %>% mutate(group=case_when(Sepal.Length>5 ~ 1, Petal.Length<3 ~ 2, TRUE ~3)) %>% 
  count(group)

tb1 <- data.frame(x=c(1,1,2), code=c(1,2,1), y=c("a","b","c"))

tb2 <- data.frame(x=c(1,1,3), code=c(1,2,1), y=c("alpha","betta","gamma"))


tb1 %>% inner_join(y = tb2,by = "x")

tb1 %>% inner_join(tb2, by=c("x", "code"))

library(tidyr)

iris_tr<- iris_tr %>% mutate(code=1:150, .before = 1)

iris_tr_longer <- iris_tr %>% pivot_longer(cols = Sepal.Length:Petal.Width,
                                           names_to = "parameter",
                                           values_to = "length")

iris_tr_longer <- iris_tr %>% pivot_longer(cols = Sepal.Length:Petal.Width,
                                           names_to = "parameter",
                                           values_to = "length")


iris_tr_longer %>% pivot_wider(names_from = "parameter", values_from = "length")


iris_tr_longer <- iris_tr_longer %>% mutate(double_length=length * 2)

iris_tr_longer %>% pivot_wider(names_from = parameter, values_from = length )

(iris_tr_longer %>% pivot_wider(names_from = parameter, 
                               values_from = c(length, double_length) ))

View(iris_tr_longer %>% pivot_wider(names_from = parameter, 
                                    values_from = c(length, double_length) ))

library(writexl)
library(readxl)
library(haven)
library(odbc)

write_xlsx(iris_tr, "C:\\Users\\user\\Desktop\\iris_tr.xlsx")

Baza <- read_excel("C:\\Users\\user\\Desktop\\Baza.xlsx")

Baza %>% group_by(REG.N) %>% summarise(cnt=n(),mn=mean(D1), mnc=mean(C8))


####*Day11####

library(dplyr)

iris_tr %>% group_by(Species, size)%>% summarise(sm=sum(Sepal.Length, na.rm = T),
                      sd=sd(Sepal.Length, na.rm=T))


iris_tr %>% group_by(size)%>% summarise_if(.predicate = is.numeric, .funs = c(mn=mean, sdev=sd))


iris_tr %>% group_by(size)%>% 
  summarise_if(.predicate = is.numeric, .funs = c(mn= ~mean(.x, trim=0.2), 
                                                  sdev=sd , cnt= ~ n()))

View(iris_tr %>% group_by(size)%>% 
  summarise_if(.predicate = is.numeric, .funs = c(mn= ~mean(.x, trim=0.2), 
                                                  sdev=sd , cnt= ~ n())))
 
                                                  
iris_tr %>% group_by(size)  %>% summarise_at(.vars = c(3,4), mean)                                               
                                                  
summarize_at(.tbl = group_by(.data = iris_tr,Species), .vars = 3:4,.funs = mean)                                                  


iris_tr %>% filter(Sepal.Length>5) %>% select(2:4) %>% arrange(Sepal.Length) %>% 
  slice(1:5)


slice(.data = arrange(.data = select(.data = filter(.data = iris_tr, 
       Sepal.Length>5)  , 2:4), Sepal.Length ), 1:5)


iris_tr %>% pull(Sepal.Length) %>% as.data.frame()

tb1
tb2

tb1%>% merge(y = tb2, by="x")

                                                  
                                                  
iris_tr %>% group_by(Species)%>% 
  summarise(across(.cols=c(Sepal.Length,Petal.Length),
  .fns=c(mn=mean, sd=sd)), ssd=sd(Sepal.Width))
                                                
                                                  
iris_tr %>% mutate(var1=mean(c(Sepal.Length, Sepal.Width))) %>% slice(1:5)

iris_tr %>% rowwise()%>%  mutate(var1=mean(c(Sepal.Length, Sepal.Width))) %>% slice(1:5)

(iris_tr %>% rowwise()%>%  mutate(var1=mean(c(Sepal.Length, Sepal.Width))))[1:5,]

(iris_tr %>% rowwise()%>%  mutate(var1=mean(2:5)))[1:5,]

(iris_tr %>% rowwise()%>%  mutate(var1=mean(c_across(2:5))))[1:5,]

(iris_tr %>% rowwise()%>%  mutate(var1=mean(c_across(Sepal.Length:Petal.Width))))[1:5,]


iris_tr %>% mutate(size1=case_when(Sepal.Length>5~1, Petal.Length<4~2, TRUE~3))

table(iris_tr %>% mutate(size1=case_when(Sepal.Length>5~1, Petal.Length<4~2, TRUE~3)) %>% 
      select(size1))

table(iris_tr %>% mutate(size1=case_when(Sepal.Length>5~1, Petal.Length<4~2, FALSE~3)) %>% 
        select(size1) , useNA ="always")


####*Day12####

library(dplyr)

f1 <- function(x){x+1}
f1(4)

f2 <- function(x,y){x^2+y^3}
f2(x = 5, y = 4)
f3 <- function(x,y=2){x*y}
f3(x=5)

round1 <- function(x){round(x,digits = 1)}
round1(pi)

round2 <- function(x, dig=1){round(x,digits = dig)}
round2(pi)
round2(pi,2)
round2(pi,3)

sum(1:3,2:5)
median(1:3,2:5)
?sum
?median

med1 <- function(...){median(c(...))}
med1(1:3,100:101)

c(1,1,3) 
data.frame(x=c(1,1,3,3,5)) %>% 
  group_by(x) %>% 
  summarise(co=n()) %>%
  arrange(desc(co)) %>% 
  slice(1) %>% 
  select(x) %>% 
  as.numeric()
  
c(1,1,3,3,5) 

data.frame(x=c(1,1,3,3,5)) %>% 
  group_by(x) %>% 
  summarise(co=n()) %>%
  arrange(desc(co)) %>% 
 filter(co==max(co)) %>% 
  select(x) %>% 
  c() %>% 
  unlist()

f_mode <- function(vec){
  data.frame(x=vec) %>% 
    group_by(x) %>% 
    summarise(co=n()) %>%
    arrange(desc(co)) %>% 
    filter(co==max(co)) %>% 
    select(x) %>% 
    c() %>% 
    unlist()
}

f_mode(rpois(n = 100,lambda = 20))


func1 <- function(x){x^2}
func2 <- function(x){
  b <- x^2
  return(b)
  }

all(func1(1:100)==func2(1:100))

cbc_f <- function(x){x^3}
int <- integrate(f = cbc_f,lower = 0, upper = 1)
int$value


plot(rnorm(n = 10^2,mean = 10,sd = 1), 
     rexp(n = 10^2,rate = 1/10)
     )

plot(rnorm(n = 10^2,mean = 10,sd = 1), 
     rexp(n = 10^2,rate = 1/10), type = "l"
)

plot(seq(0,1,0.001), cbc_f(seq(0,1,0.001)))

plot(seq(0,1,0.001), cbc_f(seq(0,1,0.001)), type = "l")

plot(seq(0,1,0.001), cbc_f(seq(0,1,0.001)), type = "l",
     xlab = "სიღარიბის ქულა",
     ylab = "საშუალო ხარჯი" ,
     main = "ფუნქციის გრაფიკი"
     )

####*Day13####

library(dplyr)

library(moments)

library(VGAM)

skew_own <- function(x){
  mean(((x-mean(x))/sd(x))^3)
}

skew_own_1 <- function(x){
  mean(((x-mean(x))/sqrt(mean((x-mean(x))^2)))^3)
}

skew_own(c(1,2,100))
skewness(c(1,2,100))

skew_own_1(c(1,2,100))

skewness(rnorm(n = 10^6, mean = 100, 17))

skewness(rnorm(n = 1000, mean = 100, 17))


skewness(rpois(n = 10^6, lambda = 0.01))
table(rpois(n = 10^6, lambda = 0.01))


rand_Par <- rpareto(n = 10^6,scale = 20, shape = 10)
hist(rand_Par, breaks = 100)


skew_pareto <- function(x){
  2*(1+x)/(x-3)*sqrt((x-2)/x)
}

rpareto(n = 10^6,scale = 5, shape = 10)
skewness(rpareto(n = 10^6,scale = 5, shape = 10))

skew_pareto(10)


skewness(rpareto(n = 10^6,scale = 5, shape = 3.5))

skew_pareto(3.5)


rand_num <- rnorm(10^6)^2+rnorm(10^6)^2+rnorm(10^6)^2+rnorm(10^6)^2

hist(rand_num)
dns <- density(rand_num, bw = 1,kernel = "gaussian")

?density

dns$x

plot(dns)

rand_num1 <- rnorm(10^2)^2+rnorm(10^2)^2+rnorm(10^2)^2+rnorm(10^2)^2
dns1 <- density(rand_num1, bw = 1,kernel = "gaussian")

plot(dns1)

dens3 <- density(rnorm(1000))
dens3$y==max(dens3$y)
dens3$x[dens3$y==max(dens3$y)]



mode_kerner <- function(x){
  dens <- density(x)
  val <- dens$x[dens$y==max(dens$y)]
  return(val)
}


mode_kerner(rand_num)

skew_pearson <- function(x){
  (mean(x)-mode_kerner(x))/sd(x)
}

skew_pearson(rpois(n = 10^6,lambda = 100))
skewness(rpois(n = 10^6,lambda = 100))


skew_pearson(rand_num)
1/sqrt(2)



####*Day14####

for(i in 1:10){
  i+1
}

for(i in 1:10){
  print(i+1)
}

dt <- data.frame(x=NA)
for(i in 1:10){
  dt[i,1] <- i+1
}


for(i in 15:20){
  dt[i-14,1] <- i+1
}

dt <- data.frame(x=NA)
a <- 1
for(i in c(1,10,15,100)){
  dt[a,1] <- i+1
  a <- a+1
}


dt <- data.frame(x=NA)

for(i in c(1,10,15,100)){
  dt[which(i==c(1,10,15,100)),1] <- i+1
}

dt
#foris es Cawera aris kargi
for (variable in vector) {
  
}

#gamravlebis tabula

results <- matrix(NA,nrow = 10,ncol = 10)
for (i in 1:10) {
  for (j in 1:10) {
   results[i,j] <- i*j
  }
}

#shemtxveviti ricxvebis matrica

res <- as.data.frame(matrix(data = NA,nrow = 10^3,ncol = 10))

for (i in 1:10) {
  for (j in 1:10^3) {
    res[j,i] <- rpois(n = 1,lambda = i)
    
  }
}

for (i in 1:10) {
  for (j in 1:10^3) {
    res[j,i] <- rpois(n = 1,lambda = i)
    print(paste(i,"_",j))
  }
}


r_norm <- rnorm(n = 10^3,mean = 10,sd = 3)

data_norm <- data.frame(n=NA, mean=NA)
for (i in 1:10^3) {
  data_norm[i,1] <- i
  data_norm[i,2] <- mean(r_norm[1:i])
  print(i)
}

data_norm

plot(data_norm$n,data_norm$mean,type = "l")




r_cauchy <- rcauchy(n = 10^4,location = 0,scale = 2)

data_cauchy <- data.frame(n=NA, mean=NA)
for (i in 1:10^4) {
  data_cauchy[i,1] <- i
  data_cauchy[i,2] <- mean(r_cauchy[1:i])
  print(i)
}

data_cauchy

plot(data_cauchy$n,data_cauchy$mean,type = "l")


hist(rcauchy(n = 10^6,location = 0,scale = 2))

#damoukideblobis shemocmeba

rand_norm <- rnorm(n = 10^3,mean = 10,sd = 2)
rand_exp <- rexp(n = 10^3,rate = 1/10)
rand_db <- data.frame(norm=rand_norm,exp=rand_exp)

hist(rand_norm)
hist(rand_exp)

results <- data.frame(x=NA, y=NA, m=NA, res=NA)

a <- 1
for (n in 1:2) {
  for (m in 1:2) {
    for (i in seq(5,7,m)) {
      for (j in seq(0.1,2.1,n)) {
       results[a,1] <- i
       results[a,2] <- i+m
       results[a,3] <- j
       results[a,4] <- j+n
       results[a,4] <- (rand_db %>% filter(norm>i & norm<i+m & exp>j & exp<j+n) %>% nrow())/10^6-
         
         ((rand_db %>% filter(norm>i & norm<i+m) %>% nrow())/10^6)*
         ((rand_db %>% filter(exp>j & exp<j+n) %>% nrow())/10^6)
       a <- a+1
       print(paste0(n,"-",m,"-",i,"-",j))
      }
    } 
    
  }
}


results

max(abs(results$res))



####*Day15####

cor(iris_tr$Sepal.Length, iris_tr$Sepal.Width)
cor(iris_tr$Sepal.Length, iris_tr$Sepal.Width, method = "kendall")
cor(iris_tr$Sepal.Length, iris_tr$Sepal.Width, method = "spearman")

cor(iris_tr$Sepal.Length, iris_tr$Sepal.Width, method = "k")
cor(iris_tr$Sepal.Length, iris_tr$Sepal.Width, method = "s")

cor(iris_tr[,1:4],method = "s")
cor(iris_tr[,1:4],method = "k")

cor_ress <- list()
co <- 1
for (t in c("p","s","k")) {
  cor_ress[[co]] <- cor(iris_tr[,1:4],method = t)
  co <- co+1
}


cor_ress <- list()
co <- 1
for (t in c("p","s","k")) {
  cor_ress[[co]] <- cor(iris[,1:4],method = t)
  co <- co+1
}

cor_ress


rank(c(4,1,100))

cor(iris$Sepal.Length,iris$Sepal.Width, method = "s")
cor(rank(iris$Sepal.Length),rank(iris$Sepal.Width), method = "p")


hist(rchisq(n = 10^6,df = 1))
hist(rchisq(n = 10^6,df = 2), breaks = 100)
hist(rchisq(n = 10^6,df = 20), breaks = 100)
hist(rchisq(n = 10^6,df = 40), breaks = 100)


hist(rt(n = 10^6,df = 1), breaks = 100)
hist(rt(n = 10^6,df = 2), breaks = 100)
hist(rt(n = 10^6,df = 3), breaks = 100)
hist(rt(n = 10^6,df = 20), breaks = 100)
hist(rt(n = 10^6,df = 40), breaks = 100)


hist(rf(n = 10^6,df1 = 1, df2 = 1 ), breaks = 100)
hist(rf(n = 10^6,df1 = 10, df2 = 10 ), breaks = 100)
hist(rf(n = 10^6,df1 = 30, df2 = 30 ), breaks = 100)


dnorm(x = 0,mean = 0,sd = 1)
integrate(f = function(x){dnorm(x, mean = 10, sd = 2)},
          lower=11 , upper =12)

plot(seq(-3,3,0.01), dnorm(seq(-3,3,0.01)), type = "l")
lines(seq(-3,3,0.01), dt(seq(-3,3,0.01), df=1) , col="red")
lines(seq(-3,3,0.01), dt(seq(-3,3,0.01), df=7) , col="blue")

pnorm(q = 12,mean = 10,sd = 2) - pnorm(q = 11,mean = 10,sd = 2)

integrate(f = function(x){x*dchisq(x,df = 10)},
          lower=0 , upper =Inf)

integrate(f = function(x){1-pchisq(x,df = 10)},
          lower=0 , upper =Inf)


####Day16####

#empiriuli ganacilebis funqcia

cf <- ecdf(rnorm(10^2))
cf
cf(0)
plot(seq(-3,3,0.01),cf(seq(-3,3,0.01)), type = "l")
lines(seq(-3,3,0.01),pnorm(seq(-3,3,0.01)), col = "red")

#qvantilrbi

qnorm(p = 0.975, mean =0,sd = 1 )
qnorm(p = 0.995, mean =0,sd = 1 )

qt(p = 0.995, df = 5)
qt(p = 0.995, df = 300)

#centraluru zgvariti teotema

dt <- c(NA)
for (i in 1:1000) {
  r <- rexp(n = 10, rate = 1/10)
  dt[i] <- (mean(r) - 10)/sd(r)/sqrt(10)
  print(i)
}
hist(dt)


dt <- c(NA)
for (i in 1:1000) {
  r <- rexp(n = 1000, rate = 1/10)
  dt[i] <- (mean(r) - 10)/sd(r)/sqrt(1000)
  print(i)
}
hist(dt)

####Day17####

#ndobis intervali

vec <- rnorm(n = 10,mean = 100,sd = 15)
round(
c(mean(vec)-sd(vec)*qt(p = 1-0.05/2, df = length(vec)-1)/sqrt(length(vec)),
  mean(vec)+sd(vec)*qt(p = 1-0.05/2,df = length(vec)-1)/sqrt(length(vec)))
, digits = 2)

library(confintr)
ci_mean(x = vec,probs = c(0.25, 0.975),type = "t")
ci <- ci_mean(x = vec,probs = c(0.25, 0.975),type = "t")
ci$interval

res <- data.frame(x=NA)
for (i in 1:1000) {
  ver <- rnorm(n = 10^4,mean = 100,sd = 15)
  res[i,1] <- 
  100>=ci_mean(ver)$interval[1]&
    100<=ci_mean(ver)$interval[2] 
}
prop.table(table(res$x))


####Day18####

#hipotezis shemocmeba

vec <- rnorm(n = 10000,mean = 10,sd = 2)
mean(vec)

v <- (mean(vec)-10.1)/(sd(vec)/sqrt(length(vec)))
2*pt(q = v,df = length(vec)-1)

2*ifelse(
  v<0,
  2*pt(q = v,df = length(vec)-1),
  1-pt(q = v,df = length(vec)-1)
)

t.test(x = vec,mu = 10.1)

a <- rnorm(n = 10000,mean = 10,sd = 2)
b <- sort(a)
plot(a,b)  
cor(a,b)  
cor.test(x = a,y = b)  
cor.test(x = a,y = b,conf.level = 0.9) 

a <- rnorm(n = 100,mean = 10,sd = 2)
b <- rexp(n = 1000,rate = 1/10)
t.test(x = a,y = b,conf.level = 0.95)


vec <- rt(n = 1000,df = 2)
ks.test(x =vec,y = "pnorm" )

vec <- rt(n = 1000,df = 10)
ks.test(x =vec,y = "pnorm" )
  
vec <- rt(n = 1000,df = 30)
ks.test(x =vec,y = "pnorm" )

vec <- rt(n = 1000,df = 100)
ks.test(x =vec,y = "pnorm" )

vec <- rt(n = 1000,df = 100)
fun <- function(x){
    pnorm(q = x,mean = mean(vec), sd = sd(vec))-ecdf(vec)(x)
  
  }

plot(seq(-10,10,0.01), fun(seq(-10,10,0.01)), type = "l")

#wrfivi regresia

cor(iris[,1:4])
plot(iris$Petal.Length,iris$Petal.Width )
linear_mod <- lm(formula = Petal.Width ~ Petal.Length, data = iris)
summary(linear_mod)
abline(linear_mod, col="red")


