# Problem Set 2
I finish five questions and summarize all answers, including equations, graphs and table, in result/HW2_Report.pdf.

In code folder, there are three Matlab code files:
   
   VFIstochstic.m: Code for Q2 and Q3 follows VFIdeterministic.m Prof. Kuhn uploded by doing VFI to plot value function, policy function and savings;
   
   Simulation.m: Code for Q4 adjusted A^h to get A^l by using long-run probability and do VFI to find policy function. Then simulate sequence of {A_t} and get {K_t+1} (also {K_t}) starting from K_0=30 to generate sequence {Y_t}. If sd(Y)<1.8% then keep [A^h, A^l]. Otherwise, go back to adjust A^h until sd(y)<1.8%. From my code, when A^h=1.002 and A^l=0.9936, sd(y)=1.75%.
   
   Extra_loop.m: Code for Q5 follows VFI idea but apply loop indexing row by i and colunm by j rather than vectorization. By using toc and tic on VFIstochastic.m and Extra_loop.m, vectorization (23.6378s) is four times faster than loop (135.3539).
