%{ 
Author: Dongheon Han
Version: 1
I worked on the homework assignment alone, using only course materials
%}
%% RK 2nd: Heun's, Mod.Euler, Ralston
clc, clear;
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Input: 1st Order Differential equation %
dfx = @(x,fx) 1.9*x;
fprintf('CHECK your function!\n')
dfx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. RK_Heun's
fprintf('***1. RK_Heun''s\n')
n = input('How many steps: ');
stepSize_h = input('Step size(unitSize): ');
fprintf('***Initial Condition***\n');
x_i = input('Enter x_0: ');
fx_i = input('Enter fx_0: ');

for i = 1 : n
    k_1 = dfx(x_i,fx_i);
    k_2 = dfx(x_i+stepSize_h,fx_i+k_1*stepSize_h);
    fx_i = fx_i + (1/2)*(k_1+k_2)*(stepSize_h);
    fprintf('\n*********\nk_1 @ n=%d: %d\n',n , k_1);
    fprintf('k_2 @ n=%d: %d\n',n , k_2);
    fprintf('fx_%d: %d\n',n , fx_i);
end

%% 2. RK_Mod.Euler
fprintf('***2. RK_Mod.Euler\n')
n = input('How many steps: ');
stepSize_h = input('Step size(unitSize): ');
fprintf('***Initial Condition***\n');
x_i = input('Enter x_0: ');
fx_i = input('Enter fx_0: ');
for i = 1 : n
    k_1 = dfx(x_i,fx_i);
    k_2 = dfx(x_i+(1/2)*stepSize_h,fx_i+(1/2)*k_1*stepSize_h);
    fx_i = fx_i + (k_2)*(stepSize_h);
    fprintf('\n*********\nk_1 @ n=%d: %d\n',n , k_1);
    fprintf('k_2 @ n=%d: %d\n',n , k_2);
    fprintf('fx_%d: %d\n',n , fx_i);
end

%% 3. RK_Ralston
fprintf('***2. RK_Ralston\n')
n = input('How many steps: ');
stepSize_h = input('Step size(unitSize): ');
fprintf('***Initial Condition***\n');
x_i = input('Enter x_0: ');
fx_i = input('Enter fx_0: ');
for i = 1 : n
    k_1 = dfx(x_i,fx_i);
    k_2 = dfx(x_i+(3/4)*stepSize_h,fx_i+(3/4)*k_1*stepSize_h);
    fx_i = fx_i + ((1/3)*k_1 + (2/3)*k_2)*stepSize_h;
    fprintf('\n*********\nk_1 @ n=%d: %d\n',n , k_1);
    fprintf('k_2 @ n=%d: %d\n',n , k_2);
    fprintf('fx_%d: %d\n',n , fx_i);
end