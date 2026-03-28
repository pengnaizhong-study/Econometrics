# 异质性稳健DID估计方法比较

## 概述

本项目使用Stata模拟数据，演示并比较六种主流的异质性稳健双重差分（DID）估计方法。在交错处理（staggered treatment）且处理效应存在异质性的设定下，传统TWFE估计会产生偏误，而新近发展的稳健估计方法可以正确识别因果效应。

## 比较的方法

| 方法 | Stata命令 | 参考文献 |
|------|-----------|----------|
| 传统TWFE | `reghdfe` | — |
| Callaway & Sant'Anna | `csdid` | Callaway & Sant'Anna (2021, JoE) |
| Borusyak, Jaravel & Spiess | `did_imputation` | Borusyak, Jaravel & Spiess (2024, ReStud) |
| de Chaisemartin & D'Haultfoeuille | `did_multiplegt` | de Chaisemartin & D'Haultfoeuille (2020, AER) |
| Gardner | `did2s` | Gardner (2022, ReStud) |
| Sun & Abraham | `eventstudyinteract` | Sun & Abraham (2021, JoE) |

## 模拟数据设计 (DGP)

- **面板结构**: 1000个个体, 20个时期
- **处理队列**: 4个交错处理队列 (t=8, 12, 15, 18) + 从未处理组
- **异质性处理效应**: `tau(g,t) = base_g + growth_g * (t - g)`
  - Cohort 8: 基础效应=5, 增长率=0.5/期
  - Cohort 12: 基础效应=3, 增长率=0.3/期
  - Cohort 15: 基础效应=2, 增长率=0.2/期
  - Cohort 18: 基础效应=1, 增长率=0.1/期

## 运行方式

```stata
do heterogeneous_did_comparison.do
```

### 前置安装

首次运行前，请取消do文件中安装命令的注释，或手动执行：

```stata
ssc install csdid
ssc install did_imputation
ssc install did_multiplegt
ssc install did2s
ssc install eventstudyinteract
ssc install reghdfe
ssc install ftools
ssc install event_plot
ssc install drdid
ssc install avar
```

## 主要参考文献

1. Callaway, B. & Sant'Anna, P.H.C. (2021). "Difference-in-Differences with Multiple Time Periods." *Journal of Econometrics*, 225(2), 200-230.
2. Sun, L. & Abraham, S. (2021). "Estimating Dynamic Treatment Effects in Event Studies with Heterogeneous Treatment Effects." *Journal of Econometrics*, 225(2), 175-199.
3. Borusyak, K., Jaravel, X. & Spiess, J. (2024). "Revisiting Event-Study Designs: Robust and Efficient Estimation." *Review of Economic Studies*, 91(6), 3253-3285.
4. de Chaisemartin, C. & D'Haultfoeuille, X. (2020). "Two-Way Fixed Effects Estimators with Heterogeneous Treatment Effects." *American Economic Review*, 110(9), 2964-2996.
5. Gardner, J. (2022). "Two-Stage Differences in Differences." *Review of Economic Studies*, 89(6), 3030-3060.
6. Goodman-Bacon, A. (2021). "Difference-in-Differences with Variation in Treatment Timing." *Econometrica*, 89(5), 2261-2290.
