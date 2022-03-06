function returnFlag = Bmat_detect_dimH(Bmat)
% check Bmat to see if the basis fucntions are linear independent in the Bmat norm fucntion space
[~,eigB] = eig(Bmat); 
eigB     = diag(eigB); returnFlag = 0; 
if isreal(eigB); Ind_zero = find(eigB<1e-20);  
    if length(Ind_zero)>=1
        n = length(Bmat(1,:));
        fprintf('\n Basis functions lead to matrix almost singular:\n ');
        fprintf('dimH= %i, number of eigenvalues <1e-20: %i \n',n,length(Ind_zero));
        % fprintf('\n Reduce basis numbers. \n'); returnFlag =1;
        %  return;
        returnFlag = 1;
    end
else
    returnFlag = 1;
    fprintf('\n Basis functions lead to matrix with complex eigenvalues. Change basis functions. \n'); returnFlag = 2;
end

end