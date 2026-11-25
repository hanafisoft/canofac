# Expected values computed with the Python reference implementation
# (same algorithms, same deterministic starts), tolerance 1e-5.

.test_blocks <- function() {
  n <- 12
  mk <- function(k, p) {
    i <- matrix(seq_len(n), n, p); j <- matrix(seq_len(p), n, p, byrow = TRUE)
    sin(0.37 * i * j + k) + 0.5 * cos(1.3 * i + 0.7 * j * k) +
      0.3 * k * sin(0.9 * i) * cos(0.4 * j)
  }
  cf_preprocess(list(A = mk(1, 4), B = mk(2, 6), C = mk(3, 9)))
}

.expected <- list(
  list(q = 1, axes = "free", extraction = "simultaneous",
       d = matrix(c(0.48001226, 0.64827785, 0.57011564, 0.37301396, 0.51672248, 0.40477526), 3, 2, byrow = TRUE),
       crit = 2.99291746),
  list(q = 2, axes = "free", extraction = "simultaneous",
       d = matrix(c(0.67530083, 0.46244300, 0.57072265, 0.36667056, 0.49539201, 0.41734710), 3, 2, byrow = TRUE),
       crit = 1.54964821),
  list(q = 4, axes = "free", extraction = "simultaneous",
       d = matrix(c(0.71842075, 0.36216771, 0.48400785, 0.48089079, 0.43477046, 0.48921852), 3, 2, byrow = TRUE),
       crit = 0.48496333),
  list(q = 4, axes = "free", extraction = "sequential",
       d = matrix(c(0.70037663, 0.41928335, 0.54749134, 0.38898757, 0.47451026, 0.45183110), 3, 2, byrow = TRUE),
       crit = 0.47664031),
  list(q = 2, axes = "orthonormal", extraction = "sequential",
       d = matrix(c(0.67530083, 0.43575312, 0.57072265, 0.34485270, 0.49539201, 0.38824229), 3, 2, byrow = TRUE),
       crit = 1.48670503),
  list(q = 2, axes = "orthonormal", extraction = "simultaneous",
       d = matrix(c(0.69883675, 0.41514278, 0.53171955, 0.39327480, 0.47573204, 0.42854830), 3, 2, byrow = TRUE),
       crit = 1.50808170),
  list(q = 1, axes = "orthonormal", extraction = "simultaneous",
       d = matrix(c(0.65518407, 0.45574263, 0.56588052, 0.36247410, 0.50287340, 0.40285317), 3, 2, byrow = TRUE),
       crit = 2.94500788),
  list(q = 3, axes = "free", extraction = "simultaneous",
       d = matrix(c(0.71720171, 0.38111417, 0.48611990, 0.47714751, 0.43870784, 0.48117420), 3, 2, byrow = TRUE),
       crit = 0.84361866)
)
