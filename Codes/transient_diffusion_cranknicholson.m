%Case Hardening Case Study
clc;
clear;

%Defining necessary boundary conditions, step size, diffusivity
D = input('Enter the diffusivity of the carbon particles in m^2/s: ');
C_o = input('Enter concentration at x = 0 mm and t = 0 s: ');
C_i = input('Enter concentration at x > 0 mm and t = 0 s: ');
t_f = input('Enter total time for diffusion in seconds: ');
C_x = input('Enter concentration at the end of the rod: ');
x_f = input('Enter the length of section of rod in m: ');
n_t = input('Enter the step size for time intervals: ');
n_x = input('Enter the step size for distance intervals: ');

delta_t = t_f/n_t; %We are taking divisions of the total time.
delta_x = x_f/n_x; 

x = (0:delta_x:x_f);

C_cn = zeros(n_x+1,n_t+1);
i = 1; %indexing for distance along the rod
l = 1; %indexing for time

C_cn(1,:) = C_o;
C_cn(2:n_x,1) = C_i;
C_cn(n_x+1,:) = C_x;

disp(C_cn);
lambda = (D*delta_t)/(delta_x)^2;

coeff_1 = 2*(lambda+1);
coeff_2 = -lambda;
coeff_mat = zeros(n_x-1,n_x-1);

%Preparation of Coefficient Matrix
while l<=n_x-1
   while i<=n_x-1
       if i == l 
           coeff_mat(i,l) = coeff_1;
       elseif ( ((l == i + 1) && (i ~= n_x-1)) || ((l == (i - 1)) && (i ~= 1 )))  
           coeff_mat(i,l) = coeff_2;
       end
       i = i + 1;
   end
   i = 1;
   l = l +1;
end

disp(coeff_mat);

%Initializing an RHS matrix
rhs_mat = zeros(n_x - 1,1);
l = 1;

coeff_3 = 2 * (1- lambda);

while l<=n_t
    while i <= n_x -1
        if i == 1
            rhs_mat(i,1) = (coeff_3 * C_cn(i+1,l))+ (lambda*(C_cn(i,l) + C_cn(i+2,l) + C_cn(i,l+1)));
        elseif i == n_x - 1 
            rhs_mat(i,1) = (coeff_3 * C_cn(i+1,l)) + (lambda*( C_cn(i,l) + C_cn(i+2,l+1) + C_cn(i+2,l))); 
        else
            rhs_mat(i,1) = lambda*(C_cn(i, l) + C_cn(i+2,l)) + (coeff_3* C_cn(i+1,l));
        end
        result = inv(coeff_mat) * rhs_mat;
        C_cn(2:n_x, l+1) = result;
        i = i + 1;
    end
i = 1;
l = l + 1;
end

l = 25;
figure; hold on
title('Crank Nicholson Method Solution');
xlabel('Distance');
ylabel('Concentration(wt%)');
while  l<= n_t + 1
    plot(x, C_cn(:, l));
    l = l + 1;
end


