%Case Hardening Case Study
clc;
clear;

%Defining necessary boundary conditions, step size, diffusivity
D = input('Enter the diffusivity of the carbon particles in m^2/s: ');
C_o = input('Enter concentration at x = 0 mm and t = 0 s: ');
C_i = input('Enter concentration at x > 0 mm and t = 0 s: ');
t_f = input('Enter total time for diffusion in seconds: ');
C_x = input('Enter concentration at the end of the rod: ');
x_f = input('Enter the length of section of rod in m:');
n_t = input('Enter the step size for time intervals: ');
n_x = input('Enter the step size for distance intervals: ');

delta_t = t_f/n_t; %We are taking divisions of the total time.
delta_x = x_f/n_x; 

C_e = zeros(n_x+1,n_t+1);
i = 2; %indexing for distance along the rod
l = 1; %indexing for time

x = (0:delta_x:x_f);

%Setting the boundary conditions
C_e(1,:) = C_o;
C_e(2:n_x,1) = C_i;
C_e(n_x+1,:) = C_x;

%Calculating value of lambda
lambda = (D*delta_t)/(delta_x)^2;

%Obtaining concentration profile of Carbon
while l <= n_t
    while i<= n_x
        C_e(i,l+1)= C_e(i,l) + (lambda * ( C_e(i+1,l) - (2*C_e(i,l)) + C_e(i-1,l)));
        i = i+1;     
    end
    l = l + 1;
    i = 2;
end

l = 25;
figure; hold on
title('Explicit Method Solution');
xlabel('Distance');
ylabel('Concentration(wt%)');
while  l<= n_t + 1
    plot(x, C_e(:, l));
    l = l + 1;
end



