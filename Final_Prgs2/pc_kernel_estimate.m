function p_hat = pc_kernel_estimate(x,z,e,h)
N=length(e);
for i=1:length(x)
    p(i) = 0;
    for j=1:length(e)  
       p(i) = p(i) + kernel((x(i)-e(j))/h); 
    end
    p_hat(i)= p(i)/(h*N);
end   
