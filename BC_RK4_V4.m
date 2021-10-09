%{ 
Author: Dongheon Han
Version: 1
I worked on the homework assignment alone, using only course materials
%}
clc, clear;
format long

dy1 = @(x,y1,y2) y2;
dy2 = @(x,y1,y2) 1.2*exp(-0.8*x)+2.7+1.6*sin(1.2*y2)-2.9*y1;
x_1 = 0.8;%reuse
fx_1 = -2.1;%reuse
x_b = 1.34;%dy_dx
fx_b = 0;%dy_dx

%dy/dx Guess in Boundary
a_im1 = (fx_b-fx_1)/(x_b-x_1);
%RK4

% use RK4
dfx_i = a_im1;
x_i = x_1;
fx_i = fx_1;
stepSize_h = 0.09;
n = round((x_b - x_1)/stepSize_h);
for i = 1 : n
    k_11 = dy1(x_i,fx_i,dfx_i);%dy1
    k_21 = dy2(x_i,fx_i,dfx_i);
    k_12 = dy1(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_11*stepSize_h, dfx_i+(1/2)*k_21*stepSize_h);%dy1
    k_22 = dy2(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_11*stepSize_h, dfx_i+(1/2)*k_21*stepSize_h);
    k_13 = dy1(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_12*stepSize_h, dfx_i+(1/2)*k_22*stepSize_h);%dy1
    k_23 = dy2(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_12*stepSize_h, dfx_i+(1/2)*k_22*stepSize_h);
    k_14 = dy1(x_i+stepSize_h, fx_i+k_13*stepSize_h, dfx_i+k_23*stepSize_h);%dy1
    k_24 = dy2(x_i+stepSize_h, fx_i+k_13*stepSize_h, dfx_i+k_23*stepSize_h);
    fx_i = fx_i + (1/6)*(k_11 + 2*k_12 + 2*k_13 + k_14)*stepSize_h;
    dfx_i = dfx_i + (1/6)*(k_21 + 2*k_22 + 2*k_23 + k_24)*stepSize_h;
    x_i = x_i + stepSize_h;
end
fa_im1 = fx_i - fx_b;


if (fx_i - fx_b < 0)
    a_i  = a_im1 + 25;
elseif (fx_i - fx_b > 0)
    a_i = a_im1 - 25 ;
else
    fprintf('Exit the code\n');
    return;
end

% use RK4
dfx_i = a_i;
x_i = x_1;
fx_i = fx_1;
stepSize_h = 0.09;
n = round((x_b - x_1)/stepSize_h);
for i = 1 : n
    k_11 = dy1(x_i,fx_i,dfx_i);%dy1
    k_21 = dy2(x_i,fx_i,dfx_i);
    k_12 = dy1(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_11*stepSize_h, dfx_i+(1/2)*k_21*stepSize_h);%dy1
    k_22 = dy2(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_11*stepSize_h, dfx_i+(1/2)*k_21*stepSize_h);
    k_13 = dy1(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_12*stepSize_h, dfx_i+(1/2)*k_22*stepSize_h);%dy1
    k_23 = dy2(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_12*stepSize_h, dfx_i+(1/2)*k_22*stepSize_h);
    k_14 = dy1(x_i+stepSize_h, fx_i+k_13*stepSize_h, dfx_i+k_23*stepSize_h);%dy1
    k_24 = dy2(x_i+stepSize_h, fx_i+k_13*stepSize_h, dfx_i+k_23*stepSize_h);
    fx_i = fx_i + (1/6)*(k_11 + 2*k_12 + 2*k_13 + k_14)*stepSize_h;
    dfx_i = dfx_i + (1/6)*(k_21 + 2*k_22 + 2*k_23 + k_24)*stepSize_h;
    x_i = x_i + stepSize_h;
end

fa_i = fx_i - fx_b;

if (fa_im1 * fa_i < 0)
    fprintf('\nSuccess Guess\n');
else 
    fprintf('\nGuess Fail\n');
    return;
end
a_ip1 = a_i - fa_i * (a_i - a_im1)/(fa_i - fa_im1);
iter = 0;
while (abs(fa_i) > 0.0001)
    % use RK4
    a_im1 = a_i;
    fa_im1 = fa_i;
    a_i = a_ip1;
    dfx_i = a_ip1;
    x_i = x_1;
    fx_i = fx_1;
    stepSize_h = 0.09;
    n = round((x_b - x_1)/stepSize_h);
    for i = 1 : n
        k_11 = dy1(x_i,fx_i,dfx_i);%dy1
        k_21 = dy2(x_i,fx_i,dfx_i);
        k_12 = dy1(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_11*stepSize_h, dfx_i+(1/2)*k_21*stepSize_h);%dy1
        k_22 = dy2(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_11*stepSize_h, dfx_i+(1/2)*k_21*stepSize_h);
        k_13 = dy1(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_12*stepSize_h, dfx_i+(1/2)*k_22*stepSize_h);%dy1
        k_23 = dy2(x_i+(1/2)*stepSize_h, fx_i+(1/2)*k_12*stepSize_h, dfx_i+(1/2)*k_22*stepSize_h);
        k_14 = dy1(x_i+stepSize_h, fx_i+k_13*stepSize_h, dfx_i+k_23*stepSize_h);%dy1
        k_24 = dy2(x_i+stepSize_h, fx_i+k_13*stepSize_h, dfx_i+k_23*stepSize_h);
        fx_i = fx_i + (1/6)*(k_11 + 2*k_12 + 2*k_13 + k_14)*stepSize_h;
        dfx_i = dfx_i + (1/6)*(k_21 + 2*k_22 + 2*k_23 + k_24)*stepSize_h;
        x_i = x_i + stepSize_h;
    end
    fa_i = fx_i - fx_b;
    a_ip1 = a_i - fa_i * (a_i - a_im1)/(fa_i - fa_im1);
    iter = iter + 1;
    fa_i
end
fprintf('iter : %d\n',iter);
dfx_i
fx_i
x_i 

