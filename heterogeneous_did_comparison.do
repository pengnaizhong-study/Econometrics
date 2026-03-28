********************************************************************************
*  异质性稳健DID估计方法比较：基于模拟数据的实证演示
*  Heterogeneity-Robust DID Estimators Comparison
*
*  比较方法：
*    1. TWFE  — 传统双向固定效应（存在异质性偏误）
*    2. csdid — Callaway & Sant'Anna (2021)
*    3. did_imputation — Borusyak, Jaravel & Spiess (2024)
*    4. did_multiplegt — de Chaisemartin & D'Haultfoeuille (2020)
*    5. did2s — Gardner (2022) 两阶段DID
*    6. eventstudyinteract — Sun & Abraham (2021)
*
*  数据设定：
*    - 1000个个体，20期面板（t=1,...,20）
*    - 交错处理（staggered treatment）：4个处理队列 + 从未处理组
*    - 异质性处理效应：随队列和时间变化
********************************************************************************

clear all
set more off
set seed 20260328

* ============================================================================
*  第一部分：安装所需命令（如尚未安装请取消注释）
* ============================================================================

/*
ssc install csdid, replace
ssc install did_imputation, replace
ssc install did_multiplegt, replace
ssc install did2s, replace
ssc install eventstudyinteract, replace
ssc install avar, replace
ssc install reghdfe, replace
ssc install ftools, replace
ssc install event_plot, replace
ssc install drdid, replace
*/

* ============================================================================
*  第二部分：生成模拟数据 (DGP)
* ============================================================================

* --- 2.1 基本面板结构 ---
local N_units = 1000          // 个体数量
local T_periods = 20          // 时间期数

set obs `=`N_units' * `T_periods''
gen id = ceil(_n / `T_periods')
bysort id: gen time = _n
xtset id time

* --- 2.2 分配处理队列（cohort）---
*  队列设定：
*   cohort = 0  : 从未处理（never-treated），约200个单位
*   cohort = 8  : 第8期开始处理，约200个单位
*   cohort = 12 : 第12期开始处理，约200个单位
*   cohort = 15 : 第15期开始处理，约200个单位
*   cohort = 18 : 第18期开始处理，约200个单位

gen cohort = .
* 每个个体只分配一次
bysort id (time): gen first = (_n == 1)
gen rand_cohort = runiform() if first == 1
bysort id: egen rc = max(rand_cohort)

replace cohort = 0  if rc <= 0.20
replace cohort = 8  if rc > 0.20 & rc <= 0.40
replace cohort = 12 if rc > 0.40 & rc <= 0.60
replace cohort = 15 if rc > 0.60 & rc <= 0.80
replace cohort = 18 if rc > 0.80 & rc <= 1.00
drop first rand_cohort rc

* 处理状态
gen treat = (cohort > 0 & time >= cohort)

* 相对处理时间（event time）
gen rel_time = time - cohort if cohort > 0
replace rel_time = -1 if cohort == 0    // never-treated 设为 -1（基准）

label var id       "个体编号"
label var time     "时间期"
label var cohort   "处理队列（首次处理时间，0=从未处理）"
label var treat    "处理状态（0/1）"
label var rel_time "相对处理时间"

* --- 2.3 生成结果变量（含异质性处理效应）---

* 个体固定效应
bysort id (time): gen alpha_i = rnormal(0, 2) if _n == 1
bysort id: replace alpha_i = alpha_i[1]

* 时间固定效应
gen delta_t = 0.5 * time + 0.05 * time^2 / `T_periods'

* 关键：异质性处理效应 —— 随队列和相对时间变化
*   tau(g, t) = base_g + growth_g * (t - g)
*   早期处理队列效应更大，且效应随时间增长

gen true_te = 0
* cohort 8:  基础效应=5，每期增长0.5
replace true_te = 5 + 0.50 * (time - cohort) if cohort == 8 & treat == 1
* cohort 12: 基础效应=3，每期增长0.3
replace true_te = 3 + 0.30 * (time - cohort) if cohort == 12 & treat == 1
* cohort 15: 基础效应=2，每期增长0.2
replace true_te = 2 + 0.20 * (time - cohort) if cohort == 15 & treat == 1
* cohort 18: 基础效应=1，每期增长0.1
replace true_te = 1 + 0.10 * (time - cohort) if cohort == 18 & treat == 1

label var true_te "真实处理效应"

* 结果变量 = 个体FE + 时间FE + 处理效应 + 噪声
gen Y = alpha_i + delta_t + true_te + rnormal(0, 1)
label var Y "结果变量"

* --- 2.4 计算真实ATT（用于基准比较）---
quietly {
    * 总体ATT
    sum true_te if treat == 1
    local true_att = r(mean)

    * 各队列ATT
    forvalues g = 8(1)18 {
        sum true_te if cohort == `g' & treat == 1
        if r(N) > 0 {
            local true_att_`g' = r(mean)
        }
    }
}

di as result _n "=============================================="
di as result "  真实参数值 (True Parameter Values)"
di as result "=============================================="
di as result "  总体 ATT        = " %6.3f `true_att'
di as result "  ATT (cohort=8)  = " %6.3f `true_att_8'
di as result "  ATT (cohort=12) = " %6.3f `true_att_12'
di as result "  ATT (cohort=15) = " %6.3f `true_att_15'
di as result "  ATT (cohort=18) = " %6.3f `true_att_18'
di as result "=============================================="

* --- 2.5 数据概览 ---
tab cohort, matrow(C) matcell(F)
tab cohort treat

di _n ">>> 数据结构摘要 <<<"
di "个体数: `N_units'  |  时间期数: `T_periods'  |  观测数: `=_N'"


* ============================================================================
*  第三部分：估计与比较
* ============================================================================

* -----------------------------------------------
*  3.1 传统 TWFE（存在偏误的基准）
* -----------------------------------------------
di as result _n "=============================================="
di as result "  方法1: 传统TWFE双向固定效应"
di as result "=============================================="

reghdfe Y treat, absorb(id time) cluster(id)
est store twfe

di as text "  TWFE 估计的 ATT = " %6.3f _b[treat]
di as text "  真实 ATT         = " %6.3f `true_att'
di as text "  偏误             = " %6.3f (_b[treat] - `true_att')

* -----------------------------------------------
*  3.2 Callaway & Sant'Anna (2021) — csdid
* -----------------------------------------------
di as result _n "=============================================="
di as result "  方法2: Callaway & Sant'Anna (2021)"
di as result "=============================================="

* csdid需要gvar（首次处理时间，never-treated=0）
gen gvar = cohort
label var gvar "首次处理时间（csdid用）"

csdid Y, ivar(id) time(time) gvar(gvar) notyet
est store cs_did

* 汇总ATT
csdid_estat simple
est store cs_simple

csdid_estat group
est store cs_group

* 事件研究
csdid_estat event
est store cs_event

* -----------------------------------------------
*  3.3 Borusyak, Jaravel & Spiess (2024) — did_imputation
* -----------------------------------------------
di as result _n "=============================================="
di as result "  方法3: Borusyak, Jaravel & Spiess (2024)"
di as result "=============================================="

did_imputation Y id time cohort, allhorizons pretrend(5) minn(0)
est store bjs

* -----------------------------------------------
*  3.4 de Chaisemartin & D'Haultfoeuille (2020) — did_multiplegt
* -----------------------------------------------
di as result _n "=============================================="
di as result "  方法4: de Chaisemartin & D'Haultfoeuille (2020)"
di as result "=============================================="

did_multiplegt Y id time treat, robust_dynamic dynamic(5) placebo(5) breps(50) cluster(id)
est store dcdh

* -----------------------------------------------
*  3.5 Gardner (2022) — did2s 两阶段DID
* -----------------------------------------------
di as result _n "=============================================="
di as result "  方法5: Gardner (2022) 两阶段DID"
di as result "=============================================="

did2s Y, first_stage(i.id i.time) second_stage(treat) treatment(treat) cluster(id)
est store gardner

* -----------------------------------------------
*  3.6 Sun & Abraham (2021) — eventstudyinteract
* -----------------------------------------------
di as result _n "=============================================="
di as result "  方法6: Sun & Abraham (2021)"
di as result "=============================================="

* 生成相对时间虚拟变量
* 限定在 [-5, +7] 窗口，基准期 = -1
forvalues k = 5(-1)2 {
    gen lead`k' = (rel_time == -`k') if cohort > 0
    replace lead`k' = 0 if cohort == 0
}
forvalues k = 0/7 {
    gen lag`k' = (rel_time == `k') if cohort > 0
    replace lag`k' = 0 if cohort == 0
}
* 基准期 lead1 (rel_time == -1) 已排除

gen never_treat = (cohort == 0)

eventstudyinteract Y lead5 lead4 lead3 lead2 lag0-lag7, ///
    vce(cluster id) absorb(id time) cohort(cohort) ///
    control_cohort(never_treat)
est store sa

matrix sa_b = e(b_iw)
matrix sa_V = e(V_iw)


* ============================================================================
*  第四部分：结果汇总与比较
* ============================================================================

di as result _n "================================================================"
di as result "      各方法ATT估计结果比较"
di as result "================================================================"
di as result "  真实ATT                                = " %7.3f `true_att'
di as result "----------------------------------------------------------------"

* TWFE
est restore twfe
di as result "  TWFE                                   = " %7.3f _b[treat] ///
    "  (偏误=" %+6.3f (_b[treat] - `true_att') ")"

* Callaway & Sant'Anna — simple ATT
est restore cs_simple
di as result "  Callaway & Sant'Anna (simple ATT)      = " %7.3f _b[ATT]

* Gardner did2s
est restore gardner
di as result "  Gardner did2s                          = " %7.3f _b[treat]

di as result "================================================================"


* ============================================================================
*  第五部分：事件研究图（Event Study Plot）
* ============================================================================

* --- 5.1 BJS (did_imputation) 事件研究图 ---
est restore bjs
event_plot bjs, default_look ///
    graph_opt( ///
        title("Borusyak, Jaravel & Spiess (2024)") ///
        subtitle("did_imputation 事件研究") ///
        xtitle("相对处理时间") ytitle("处理效应估计") ///
        name(g_bjs, replace) ///
    ) stub_lag(tau#) stub_lead(pre#) together

* --- 5.2 Callaway & Sant'Anna 事件研究图 ---
est restore cs_event
event_plot cs_event, default_look ///
    graph_opt( ///
        title("Callaway & Sant'Anna (2021)") ///
        subtitle("csdid 事件研究") ///
        xtitle("相对处理时间") ytitle("处理效应估计") ///
        name(g_cs, replace) ///
    ) stub_lag(Tp#) stub_lead(Tm#) together


* --- 5.3 综合比较图（多方法叠加）---

* 提取 BJS 系数用于手动绘图
est restore bjs
local bjs_ncoef = e(df_m) + 1

* 构建比较数据集
preserve
clear

* 用BJS的事件研究系数
est restore bjs
local ncoef = 0
foreach v of local bjs_ncoef {
    local ++ncoef
}

* 生成真实效应的事件研究线
* 使用所有处理队列的平均值
set obs 13   // rel_time from -5 to +7
gen rel_t = _n - 6   // -5 to +7

* 真实平均效应（跨队列加权）
gen true_effect = 0 if rel_t < 0
replace true_effect = (5 + 0.5*rel_t + 3 + 0.3*rel_t + 2 + 0.2*rel_t + 1 + 0.1*rel_t)/4 if rel_t >= 0

twoway (line true_effect rel_t, lcolor(black) lwidth(thick) lpattern(solid)) ///
    , legend(label(1 "真实平均效应")) ///
      title("真实处理效应路径（跨队列平均）") ///
      xtitle("相对处理时间") ytitle("处理效应") ///
      xline(-0.5, lcolor(red) lpattern(dash)) ///
      name(g_true, replace)

restore


* ============================================================================
*  第六部分：各队列异质性展示
* ============================================================================

di as result _n "================================================================"
di as result "      各队列(Cohort)真实ATT"
di as result "================================================================"

preserve
collapse (mean) true_te Y, by(cohort time treat)

twoway (line true_te time if cohort == 8,  lcolor(blue) lwidth(medthick)) ///
       (line true_te time if cohort == 12, lcolor(red) lwidth(medthick)) ///
       (line true_te time if cohort == 15, lcolor(green) lwidth(medthick)) ///
       (line true_te time if cohort == 18, lcolor(orange) lwidth(medthick)) ///
    , legend(order(1 "Cohort 8" 2 "Cohort 12" 3 "Cohort 15" 4 "Cohort 18") ///
             pos(11) ring(0) col(2)) ///
      title("各队列真实处理效应路径") ///
      subtitle("异质性处理效应：早期队列效应更大、增速更快") ///
      xtitle("时间") ytitle("真实处理效应") ///
      xline(8 12 15 18, lcolor(gs12) lpattern(dash)) ///
      name(g_cohort_te, replace)

restore


* ============================================================================
*  第七部分：TWFE偏误来源的直觉解释
* ============================================================================

di as result _n "================================================================"
di as result "  TWFE偏误的直觉解释"
di as result "================================================================"
di as text ""
di as text "  在交错DID设定下，TWFE回归实际上是多个2×2 DID的"
di as text "  加权平均。问题在于："
di as text ""
di as text "  1) '已处理 vs 已处理' 比较：TWFE会使用早期处理组"
di as text "     作为晚期处理组的'对照组'，产生负权重。"
di as text ""
di as text "  2) 异质性处理效应：当处理效应随时间或队列变化时，"
di as text "     这些负权重会导致严重偏误。"
di as text ""
di as text "  3) 本例中，早期队列(cohort=8)的效应大且增长快，"
di as text "     被TWFE用作后期队列的隐含对照组时，会向下"
di as text "     拉动ATT估计，造成低估。"
di as text ""
di as text "  异质性稳健方法通过以下策略避免偏误："
di as text "  - CS(2021)  : 只使用'尚未处理'或'从未处理'组作对照"
di as text "  - BJS(2024) : 用未处理观测估计反事实，再插补处理效应"
di as text "  - dCDH(2020): 识别并报告负权重，使用替代估计量"
di as text "  - SA(2021)  : 用交互加权估计量纠正TWFE权重问题"
di as text "  - Gardner(2022): 两阶段法，先净化FE再估计处理效应"
di as text "================================================================"


* ============================================================================
*  第八部分：保存结果
* ============================================================================

* 保存估计结果表
est table twfe cs_simple gardner, b(%7.3f) se(%7.3f) ///
    title("ATT估计比较") stats(N r2)


di as result _n "================================================================"
di as result "  程序运行完毕！"
di as result "  所有图形已保存在内存中，可用 graph dir 查看"
di as result "================================================================"

log close _all
