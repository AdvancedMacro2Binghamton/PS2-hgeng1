clear all
clc
close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
pi_hh=0.977; pi_hl=1-pi_hh; pi_ll=0.926; pi_lh=1-pi_ll;
pi=[pi_hh, pi_hl; pi_lh, pi_ll]; %% pi matrix
A_h=1.1; A_l=0.678;
A=[A_h;A_l];



clear all; a  = 0.35; b = 0.99; d = 0.025; s = 2;
z  = [.678 1.1]';                        
p  = [.926 (1-.977)]';                   
kss= (a*inv(inv(b)-(1-d)))^inv(1-a);     
k  = linspace(0.001*kss, 1.3*kss, 1000); 

t0 = tic;
u = cat(3,zeros(length(k)),zeros(length(k)));
for m = 1:length(z) 
    for i = 1:length(k) 
        for j = 1:length(k)
            c = z(m)*k(i)^a+(1-d)*k(i)-k(j);
            if  c>0; u(i,j,m)= c^(1-s)/(1-s);
            else     u(i,j,m)= -Inf; end
end; end; end

v  = cat(3,zeros(length(k)),zeros(length(k))); 
v0 = zeros(length(z),length(k));
e  = 1; 
while e>1e-06
    for m=1:length(z)
        for i = 1:length(k) 
            for j = 1:length(k)
                v(i,j,m) = u(i,j,m)+b*(p(m)*v0(1,j)+...
                                   (1-p(m))*v0(2,j));
        end; end         
        [vfn(m,:),idx(m,:)] = max(v(:,:,m),[],2);
    end
    e = max(max(abs(vfn - v0)));
    v0 = vfn;
end; t = toc(t0);