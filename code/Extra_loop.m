clear all
clc
close all
t0 = tic; %Timer start
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
pi_hh=0.977; pi_hl=1-pi_hh; pi_ll=0.926; pi_lh=1-pi_ll;
pi=[pi_hh, pi_hl; pi_lh, pi_ll]; %% pi matrix
A_h=1.1; A_l=0.678;
A=[A_h;A_l];
%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k
k = linspace(k_min, k_max, num_k);

%%%%%%%%%%%%loop%%%%%%%%%%%
%%%consumption function + reutrn function
ret_h = zeros(length(k));
ret_l = zeros(length(k));

 for i = 1:length(k) 
     for j = 1:length(k)
            cons_h(i,j) = A_h*k(i)^alpha+(1-delta)*k(i)-k(j);
            cons_l(i,j)= A_l*k(i)^alpha+(1-delta)*k(i)-k(j);
                
            ret_h(i,j)= cons_h(i,j)^(1-sigma)/(1-sigma);
            ret_l(i,j)= cons_l(i,j)^(1-sigma)/(1-sigma);
            if cons_h(i,j) < 0
                ret_h(i,j)= -Inf; 
            end
            if cons_l(i,j) < 0
                ret_l(i,j)= -Inf;
            end
       end
end; 

%%%%VFI%%%%%
v_guess = zeros(2,length(k));
dis = 1; tol = 1e-06;
while dis>tol
        for i = 1:length(k) 
            for j = 1:length(k)
                value_mat_h(i,j) = ret_h(i,j)+beta*(pi(1,1)*v_guess(1,j)+pi(1,2)*v_guess(2,j));
                value_mat_l(i,j) = ret_l(i,j)+beta*(pi(2,1)*v_guess(1,j)+pi(2,2)*v_guess(2,j));
        end; end         
    [vfn_h, pol_indx_h] = max(value_mat_h, [], 2); %max for each row
    vfn_h = vfn_h';
    
    [vfn_l, pol_indx_l] = max(value_mat_l, [], 2);
    vfn_l = vfn_l';
    
    dis = [max(abs(vfn_h - v_guess(1,:)));max(abs(vfn_l - v_guess(2,:)))];
    
    v_guess = [vfn_h;vfn_l];
end
    t = toc(t0);
    
    
% policy function
g_h = k(pol_indx_h); %High A
g_l = k(pol_indx_l); %Low A

%Saving
s_h=g_h-(1-delta)*k; %High A
s_l=g_l-(1-delta)*k; %Low A

%%%%%Plot Value function over k
figure (1)
suptitle('Value Function')
plot(k,vfn_h);
hold on;
plot(k,vfn_l);
hold off;
legend('for A^h','for A^l','location','northwest');
xlabel('k');
ylabel('V(k)');

%%%%%Plot Policy Function over k
figure (2)
plot(k,g_h);
hold on;
plot(k,g_l);
hold off;
legend('for A^h','for A^l','location','northwest')
suptitle('Policy Function')
xlabel('k');
ylabel('g(k)');

%%%%%Plot Saving Function over k
figure (3)
plot(k,s_h);
hold on;
plot(k,s_l);
hold off;
legend('for A^h','for A^l','location','northwest');
suptitle('Saving over k')
xlabel('k');
ylabel('Saving');