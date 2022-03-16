function [mf] = main_regression(x,z,e,h)
mf = pc_kernel_estimate(x,z,e,h);
mf=mf';

