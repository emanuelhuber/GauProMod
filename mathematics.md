---
layout: page
title: Maths
date: 2018-02-12
---



## Log determinant of positive definite matrices


<!--
$$\forall x \in R$$
-->
(Cholesky decomposition)

$$\mathbf{A} = \mathbf{L}\mathbf{L}'$$ 

(determinant of a lower triangular matrix)
$$|\mathbf{L}| = \prod_i L_{ii} $$

log determinant of a lower triangular matrix
$$ \log \prod_i x_i = \sum_i \log x_i$$

Thus to calculate the log determinant of a symmetric positive definite matrix:


```r
L <- chol(A);
logdetA <- 2*sum(log(diag(L)))
```

