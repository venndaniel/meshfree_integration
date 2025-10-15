% Script to generate Figure 1, with data generated using test4p1p1.m

set(groot,'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex')

% Number of points tested for q = 7/3
n73 = [400
800
1200
1600
2000
2400
2800
];

% errors for q = 7/3
err73 = [2.11E-02
5.65E-03
1.01E-03
3.83E-04
7.43E-05
7.94E-05
1.37E-05
];

% Number of points tested for q = 10/3
n103 = [400
800
1200
1600
2000
2400
2800
];


% errors for q = 10/3
err103 = [9.21E-03
1.77E-03
4.17E-04
3.61E-06
7.42E-05
1.19E-05
1.75E-05
];

% Number of points tested for triangulation and simple average
ntri= [400
800
1600
3200
6400
12800
25600
51200
];

% error for triangulation
errtri =[
 6.50E-02
5.04E-02
3.63E-03
1.72E-03
8.89E-04
4.64E-04
2.22E-04
1.14E-04
];

% error for simple average
erravg = [2.16E-02
3.02E-02
4.11E-04
1.02E-02
1.30E-02
1.46E-02
1.24E-02
1.56E-02
];

loglog(n73, err73, 'o-');
hold on;
plot(n103, err103, 'o-');
plot(ntri, errtri, 'o-');
plot(ntri, erravg, 'o-');

legend('Meshfree: $q=7/3$', 'Meshfree: $q=10/3$', 'Triangulation', 'Point Cloud Average');
grid on;
title('Surface Integral Convergence Test')
xlabel('Number of Points/Vertices')
ylabel('Relative Error');
fontsize(18, 'pixels')
