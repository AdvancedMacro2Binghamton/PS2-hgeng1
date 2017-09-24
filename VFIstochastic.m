clear all
clc
close all
%Timer
t0=tic;
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

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
    %Consumption function
cons_h = A(1)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; %at A_h
cons_l = A(2)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; %at A_l
    % return function
ret_h = cons_h .^ (1 - sigma) / (1 - sigma); % at high
ret_l = cons_l .^ (1 - sigma) / (1 - sigma); % at low
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret_h(cons_h < 0) = -Inf; % in the bracket is index like the dummary variable for location of negative consumption
ret_l(cons_l < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(2, num_k); %Firsr row is V_h, Second row is V_l.

while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_mat_h = ret_h + beta *(pi(1,1)* repmat(v_guess(1,:), [num_k 1])+pi(1,2)* repmat(v_guess(2,:), [num_k 1]));
    value_mat_l = ret_l + beta *(pi(2,1)* repmat(v_guess(1,:), [num_k 1])+pi(2,2)* repmat(v_guess(2,:), [num_k 1]));
    % find the optimal k' for every k:
    [vfn_h, pol_indx_h] = max(value_mat_h, [], 2); %max for each row
    vfn_h = vfn_h';
    
    [vfn_l, pol_indx_l] = max(value_mat_l, [], 2);
    vfn_l = vfn_l';
    
    % what is the distance between current guess and value function
    dis = [max(abs(vfn_h - v_guess(1,:)));max(abs(vfn_l - v_guess(2,:)))];
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = [vfn_h;vfn_l];
end; 
t=toc(t0); %timer end 

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



  
   
