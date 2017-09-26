# pmtest
R package that implements those modified tests

When dealing with scenarios where some samples from matched pairs design are missing, common statistical tests lack operating characteristics to robustly infer on parameters. The pmtest package provides users with five options when analyzing partially matched samples:


Modified T-Statistic by Kim et al

Corrected Z-Test by Looney and Jones

Maximum Likelihood Test by Lin and Stivers

Maximum Likelihood Test by Ekbohm

Weighted Z-test by Kuan and Huang

This document introduces you to situations where you will use these methods, and shows you how to apply them in partially matched samples situations.


Ideally, in a matched samples analysis, one would expect a total of 2n samples. However, there are many instances where samples are missing from one or both of the pairs. In partially matched samples situations, the data looks like

Case	Control
10	8
NA	4
12	NA
14	13
2	4
4	5
Note that each method takes in two one column vectors, each representing “Case” and “Control” respecitvely.

Partially matched samples can be viewed as data generated from two experimental designs where both designs intend to estimate the same parameter:

n1 matched pairs
Case	Control
10	8
14	13
2	4
4	5
independent groups with n2 and n3 per group
n2 would be the vector c(10,12,14,2,4) from the first column of the first table
n3 would be vector c(8,4,13,4,5) from the second column of the first table
