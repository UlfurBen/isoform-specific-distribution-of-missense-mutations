					     1.	     2.     3.    4.    5.    6.    7.   8.    9.   10.   11.   12.   13.   14.   15.  16
pathogenic_variants_in_EM_genes <-         c(349,   459,   123,  133,  261,  189,  432, 467,  101,  92,   42,   93,   394,  143,  35,  49)
benign_control_variants_in_EM_genes <-     c(1411,  756,   67,   168,  793,  477,  45,  297,  89,   169,  132,  222,  203,  543,  29,  213)
pathogenic_variants_in_non_EM_genes <-     c(7427,  3926,  1607, 2944, 2765, 2750, 433, 1796, 1058, 977,  289,  1115, 1791, 2405, 394, 1046)
benign_control_variants_in_non_EM_genes <- c(24712, 10221, 1574, 2782, 8483, 6451, 797, 4447, 1961, 2414, 2300, 2242, 2688, 9168, 721, 5003)

categories <- c("Nonpolar_to_Nonpolar", "Nonpolar_to_Polar", "Nonpolar_to_Negatively_charged", "Nonpolar_to_Positively_charged",
                "Polar_to_Nonpolar", "Polar_to_Polar", "Polar_to_Negatively_charged", "Polar_to_Positively_charged",
                "Negatively_charged_to_Nonpolar", "Negatively_charged_to_Polar", "Negatively_charged_to_Negatively_charged",
                "Negatively_charged_to_Positively_charged", "Positively_charged_to_Nonpolar", "Positively_charged_to_Polar",
                "Positively_charged_to_Negatively_charged", "Positively_charged_to_Positively_charged")

contingency_table <- data.frame(
  "EM genes" = c(pathogenic_variants_in_EM_genes[i], benign_control_variants_in_EM_genes[i]),
c(pathogenic_variants_in_non_EM_genes[i], benign_control_variants_in_non_EM_genes[i]),
row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )





1
Nonpolar_to_Nonpolar <- data.frame(
  "EM genes" = c(349, 1411),
  "non-EM genes" = c(7427, 24712),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Nonpolar_to_Nonpolar))

p-value = 0.001354
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.7278022 0.9287168
sample estimates:
odds ratio 
  0.822996

2
Nonpolar_to_Polar <- data.frame(
  "EM genes" = c(459, 756),
  "non-EM genes" = c(3926, 10221),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Nonpolar_to_Polar))

p-value = 4.461e-13
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.396298 1.787876
sample estimates:
odds ratio 
  1.580593 

3
Nonpolar_to_Negatively_charged <- data.frame(
  "EM genes" = c(123, 67),
  "non-EM genes" = c(1607, 1574),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Nonpolar_to_Negatively_charged))

p-value = 0.0001326
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.312944 2.479255
sample estimates:
odds ratio 
  1.797804 

4
Nonpolar_to_Positively_charged <- data.frame(
  "EM genes" = c(133, 168),
  "non-EM genes" = c(2944, 2782),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Nonpolar_to_Positively_charged))

p-value = 0.01526
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.5877914 0.9507179
sample estimates:
odds ratio 
 0.7481432

5
Polar_to_Nonpolar <- data.frame(
  "EM genes" = c(261, 793),
  "non-EM genes" = c(2765, 8483),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Polar_to_Nonpolar))

p-value = 0.9107
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.8688961 1.1707050
sample estimates:
odds ratio 
  1.009767

6
Polar_to_Polar <- data.frame(
  "EM genes" = c(189, 477),
  "non-EM genes" = c(2750, 6451),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Polar_to_Polar))

p-value = 0.4298
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.7765261 1.1091091
sample estimates:
odds ratio 
 0.9294832

7
Polar_to_Negatively_charged <- data.frame(
  "EM genes" = c(432, 45),
  "non-EM genes" = c(433, 797),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Polar_to_Negatively_charged))

p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 12.64986 25.08220
sample estimates:
odds ratio 
   17.6356

8
Polar_to_Positively_charged <- data.frame(
  "EM genes" = c(467, 297),
  "non-EM genes" = c(1796, 4447),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Polar_to_Positively_charged))

p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 3.323575 4.563824
sample estimates:
odds ratio 
   3.89243

9
Negatively_charged_to_Nonpolar <- data.frame(
  "EM genes" = c(101, 89),
  "non-EM genes" = c(1058, 1961),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Negatively_charged_to_Nonpolar))

p-value = 1.088e-06
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.549724 2.857494
sample estimates:
odds ratio 
  2.102925

10
Negatively_charged_to_Polar <- data.frame(
  "EM genes" = c(92, 169),
  "non-EM genes" = c(977, 2414),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Negatively_charged_to_Polar))

p-value = 0.03395
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.020407 1.764100
sample estimates:
odds ratio 
  1.344975

11
Negatively_charged_to_Negatively_charged <- data.frame(
  "EM genes" = c(42, 132),
  "non-EM genes" = c(289, 2300),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Negatively_charged_to_Negatively_charged))

p-value = 4.959e-06
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.707154 3.693101
sample estimates:
odds ratio 
  2.531143

12
Negatively_charged_to_Positively_charged <- data.frame(
  "EM genes" = c(93, 222),
  "non-EM genes" = c(1115, 2242),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Negatively_charged_to_Positively_charged))

p-value = 0.1884
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.6469903 1.0900665
sample estimates:
odds ratio 
 0.8423598

13
Positively_charged_to_Nonpolar <- data.frame(
  "EM genes" = c(394, 203),
  "non-EM genes" = c(1791, 2688),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Positively_charged_to_Nonpolar))

p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 2.426161 3.503809
sample estimates:
odds ratio 
  2.912332

14
Positively_charged_to_Polar <- data.frame(
  "EM genes" = c(143, 543),
  "non-EM genes" = c(2405, 9168),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Positively_charged_to_Polar))

p-value = 0.9614
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.8245454 1.2161333
sample estimates:
odds ratio 
  1.003912

15
Positively_charged_to_Negatively_charged <- data.frame(
  "EM genes" = c(35, 29),
  "non-EM genes" = c(394, 721),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Positively_charged_to_Negatively_charged))

p-value = 0.002982
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.288915 3.804765
sample estimates:
odds ratio 
  2.207009

16
Positively_charged_to_Positively_charged <- data.frame(
  "EM genes" = c(49, 213),
  "non-EM genes" = c(1046, 5003),
  row.names = c("pathogenic variants in Categorical-changes-of-interest", "control variants in Categorical-changes-of-interest"),
  stringsAsFactors = FALSE
  )

print(fisher.test(Positively_charged_to_Positively_charged))

p-value = 0.5597
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.7838165 1.5191459
sample estimates:
odds ratio 
  1.100315
