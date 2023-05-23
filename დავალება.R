#Task1
v1 <- c(1, 5, 12, 10) #c-vectori
v2 <- rep(c(0, 1), each = 10)
v3 <- seq(from = 0, to = 20, by = 0.5) #seq-aritmetikuli prog
v4 <- seq(from = 10, to = 30)
v4_1 <- v4
v4_1[c(TRUE, FALSE)] <- 0
v4_2 <- v4_1
v4_2[c(FALSE, TRUE)] <- 1

#Task2
v1 %% 2
which(v1 %% 2 == 1)

v2 %% 2
which(v2 %% 2 == 1)

#Task3
v5 <- c(v1, sum(v1))

#Task4
count <- 0
limit <- 10000000

for (i in 1:limit) {
  if (i %% 2 != 0 && i %% 3 != 0 && i %% 5 != 0 && i %% 7 != 0 && i %% 210 != 0) {
    count <- count + 1
  }
}

count

#Task5
n <- 100
matrix <- matrix(nrow = n, ncol = n)

for (i in 1:n) {
  for (j in 1:n) {
    matrix[i, j] <- i / j
  }
}

matrix

#Task6
v3_6 <- seq(from = 0, to = 20, by = 0.5)
v4_6 <- 10:30

v6 <- v3_6 * v4_6
v6
p1 <- sum(v6 < 100 | v6 > 300)
p2 <- sum(v6 >= 200 & v6 <= 400)
p3 <- sum(head(v6, 10)) + sum(tail(v6, 10))
p1
p2
p3

#Task7?





#task8
v1 <- c(1, 5, 12, 10)
unknown <- v1[3]

#task9
fifteenth <- floor(pi * 10^14) %% 10

fifteenth

