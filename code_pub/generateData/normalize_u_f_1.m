function [u,f,data_str] = normalize_u_f(u,f,dx,normalizeOn)
% normalize u and f in L u=f
[case_num,tN,~] = size(u); 

if normalizeOn == 1
    for t=1:tN
       for n = 1:case_num
           temp = u(n,t,:); 
           L2norm = sum(temp.^2)*dx; 
           u(n,t,:) = temp/L2norm; 
           f(n,t,:) = f(n,t,:)/L2norm; 
       end
    end
    data_str = 'NormalizeOn';
else
    data_str = 'NormalizeOff';
end
end