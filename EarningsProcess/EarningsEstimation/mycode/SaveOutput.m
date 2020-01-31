ComputeMoments

if SaveSimulations
    f=fopen([OutputDir 'yannsim.txt'],'w');
    fprintf(f,'%20.14f,%20.14f,%20.14f,%20.14f,%20.14f\n',yannsim');
    fclose(f);
end

f=fopen([OutputDir 'parameters.txt'],'w');
fprintf(f,'lambda %g %g\n',guess(1,:));
fprintf(f,'zetaP %g %g\n',guess(2,:));
fprintf(f,'zetaN %g %g\n',guess(3,:));
fprintf(f,'rho %g %g\n',guess(4,:));
fprintf(f,'sigma %g %g\n',guess(5,:));
fprintf(f,'delta %g %g\n',guess(6,:));
fclose(f);

f=fopen([OutputDir 'moments.txt'],'w');
fprintf(f,'muy %g\n',muy);
fprintf(f,'mu2y %g\n',mu2y);
fprintf(f,'mu3y %g\n',mu3y);
fprintf(f,'mu4y %g\n',mu4y);
fprintf(f,'gam3y %g\n',gam3y);
fprintf(f,'gam4y %g\n',gam4y);

fprintf(f,'muylev %g\n',muylev);
fprintf(f,'mu2ylev %g\n',mu2ylev);
fprintf(f,'mu3ylev %g\n',mu3ylev);
fprintf(f,'mu4ylev %g\n',mu4ylev);
fprintf(f,'gam3ylev %g\n',gam3ylev);
fprintf(f,'gam4ylev %g\n',gam4ylev);

fprintf(f,'mu2dy1 %g\n',mu2dy1);
fprintf(f,'mu3dy1 %g\n',mu3dy1);
fprintf(f,'mu4dy1 %g\n',mu4dy1);
fprintf(f,'gam3dy1 %g\n',gam3dy1);
fprintf(f,'gam4dy1 %g\n',gam4dy1);

fprintf(f,'mu2dy5 %g\n',mu2dy5);
fprintf(f,'mu3dy5 %g\n',mu3dy5);
fprintf(f,'mu4dy5 %g\n',mu4dy5);
fprintf(f,'gam3dy5 %g\n',gam3dy5);
fprintf(f,'gam4dy5 %g\n',gam4dy5);
fclose(f);
