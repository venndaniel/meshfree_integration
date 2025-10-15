% Script to generate Figure 4, with data generated using test4p2p1.m

% set plot text to LaTex
set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

% Number of points tested
Ns = [10
20
40
80
160
320
640
1280
2560
]; 

% recorded relative errors
errs = [1.0108E-01
5.7805E-03
3.3604E-03
1.8219E-04
3.7362E-06
3.1831E-06
9.2365E-08
6.7361E-10
1.1281E-12
];


% Create plot
loglog(Ns, errs, 'o-');
title('2D Singular Integral Convergence Test')
fontsize(18,'pixels')
hold on
loglog(Ns, 2*errs(1)*(Ns/Ns(1)).^-4, 'k--');
legend("Relative Error", "$\mathcal{O}(\tilde{N}^{-4}) = \mathcal{O}(h_{\textrm{max}}^{8})$")
grid on;
xlim([8, 2560*1.2])
xlabel("$\tilde{N}$");
ylabel('Relative Error')