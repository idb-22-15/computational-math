runge_kutta4 <- function (f, a, b, y0, n) {
  m <- length(y0)
  h <- (b - a)/n
  x <- seq(a + h, b, by = h)
  y <- matrix(0, nrow = n, ncol = m)
  k1 <- h * f(a, y0)
  k2 <- h * f(a + h/2, y0 + k1/2)
  k3 <- h * f(a + h/2, y0 + k2/2)
  k4 <- h * f(a + h, y0 + k3)
  y[1, ] <- y0 + k1/6 + k2/3 + k3/3 + k4/6
  for (i in 1:(n - 1)) {
    k1 <- h * f(x[i], y[i, ])
    k2 <- h * f(x[i] + h/2, y[i, ] + k1/2)
    k3 <- h * f(x[i] + h/2, y[i, ] + k2/2)
    k4 <- h * f(x[i] + h, y[i, ] + k3)
    y[i + 1, ] <- y[i, ] + k1/6 + k2/3 + k3/3 + k4/6
  }
  x <- c(a, x)
  y <- rbind(y0, y)
  return(list(x = x, y = y))
}


H <- 500  # Глубина погружения 
rho0 <- 1000
rho1 <- 750
g <- 9.8
S_sech <- 15.5 * 10  # Площадь поперечного сечения подводной лодки
l <- 100.5  # Длина подводной лодки
V <- S_sech*l
alpha <- 0.01
eta <- 0.001
v <- 25  # Скорость лодки (в м/с)
k <- S_sech / l
# Функция для системы уравнений
submarine_sys <- function(t, y) {
  dy <- y[2]
  dz <- - (eta * k / (rho1*V)) * (1 + alpha * y[1] / H) * dy + g * (rho0 / rho1 - 1)
  return(c(dy, dz))
}
# Начальные условия
y0 <- c(-H, 0)  # y(0) = -H, v_y(0) = 0
a <- 0
b <- 100
n <-300
# Решение системы дифференциальных уравнений
solution <- runge_kutta4(submarine_sys, a, b, y0, n)
# Убираем значения, где y1 >= 0
solution$x <- solution$x[solution$y[, 1] <= 0]
solution$y <- solution$y[solution$y[, 1] <= 0, ]


# Функция для аппроксимации методом наименьших квадратов
custom_lm <- function(x, y, degree) {
  n <- length(x)
  X <- matrix(1, n, degree + 1)
  for (i in 1:degree) {
    X[, i + 1] <- x^i
  }
  coefficients <- solve(t(X) %*% X) %*% t(X) %*% y
  return(coefficients)
}

# Аппроксимация решения с помощью нашей функции
degree <- 2    # Степень полинома
coefs <- custom_lm(solution$x, solution$y[, 1], degree)

# Получение коэффициентов полинома
a <- coefs[3]
b <- coefs[2]
c <- coefs[1]

# Создание данных для аппроксимированной кривой
t_fit <- seq(0, max(solution$x), length.out = 100)
y_fit <- a * t_fit^2 + b * t_fit + c
fit_data <- data.frame(t = t_fit, y = y_fit)


# Построение графика траектории всплытия
library(ggplot2)
ggplot(data = data.frame(x = solution$x * v, y = solution$y[, 1]), aes(x = x, y = y)) +
  geom_point(color = "blue") +
  geom_line(data = fit_data, aes(x = t * v, y = y), color = "green") +
  ggtitle("Boat") +
  xlab("X") +
  ylab("Y") +
  theme_minimal()


# Определение времени всплытия и точки всплытия
solve_quadratic <- function(a, b, c, H) {
  roots <- polyroot(c(c - H, b, a))
  print(roots)
  real_roots <- Re(roots[abs(Im(roots)) < 1e-6])
  return(real_roots[real_roots > 0])
}
T <- solve_quadratic(a, b, c, 0)
L <- v * T
# Вывод значений времени всплытия и точки всплытия
cat("Время всплытия (T):", T, "\n")
cat("Абсцисса точки всплытия (L):", L, "\n")




