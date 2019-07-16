k1 = 45.8;
k2 = 27.2;
lk1 = log(k1);
lk2 = log(k2);
A = 161.6000;
B = log(12000);


a = A*lk1*lk2-A*(lk1+lk2)*B+A*B^2;
b = A*(lk1+lk2)-2*A*B-(k1*lk2-k2*lk1-(k1-k2)*B);
c = A-(k1-k2);

x = (-b + sqrt(b^2-4*a*c)) / (2*a);

B = log(20000);
f_cor = @(x,epsilon) (x./(epsilon.*log(x)+1-B*epsilon));


g = @(eps_1) ( (f_cor(2753.8,eps_1) - f_cor(2592.2,eps_1)) ./ (f_cor(45.8,eps_1) - f_cor(27.2,eps_1)) );

f_cor = @(x,epsilon) (x./(epsilon.*log(x)));
g = @(eps_1) ( (f_cor(10189,eps_1) - f_cor(9384,eps_1)) ./ (f_cor(7976,eps_1) - f_cor(7369,eps_1)) );