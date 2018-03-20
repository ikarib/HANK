clear
load ss neqmiter wage netwage profit ra KNratio KYratio V gmat
V=reshape(V,2000,33); gmat=reshape(permute(reshape(gmat,40,50,33),[2 1 3]),2000,33);
fid = fopen('V.txt','w');
fprintf(fid,'%d\n',neqmiter);
fprintf(fid,'%26.15G\n',[wage netwage profit ra KNratio KYratio]);
fprintf(fid,[repmat('%26.15G',1,2000) '\n'],[V gmat]);
fclose(fid);