clear
load ss neqmiter wage netwage profit ra KNratio KYratio V gmat
V=reshape(V,2000,33); gmat=reshape(permute(reshape(gmat,40,50,33),[2 1 3]),2000,33);
fid = fopen('V.txt','w');
fprintf(fid,'%d\n',neqmiter);
fprintf(fid,'%26.15G\n',[wage netwage profit ra KNratio KYratio]);
fprintf(fid,[repmat('%26.15G',1,2000) '\n'],[V gmat]);
fclose(fid);



%% transition
clear
load trans
fid = fopen('TRANS/sol.dat','w');
fprintf(fid,[repmat('%26.15G',1,200) '\n'],[capital; labor; pi; ra]');
fclose(fid);


%% 
clear
load ../FortranOutputDir/BaselineOutputSubdir/IRF_Monetary/NOFS/STICKY/capital.txt
load ../FortranOutputDir/BaselineOutputSubdir/IRF_Monetary/NOFS/STICKY/labor.txt
load ../FortranOutputDir/BaselineOutputSubdir/IRF_Monetary/NOFS/STICKY/pi.txt
load ../FortranOutputDir/BaselineOutputSubdir/IRF_Monetary/NOFS/STICKY/ra.txt
fid = fopen('TRANS/sol.dat','w');
fprintf(fid,[repmat('%26.15G',1,4) '\n'],[capital labor pi ra]');
fclose(fid);
