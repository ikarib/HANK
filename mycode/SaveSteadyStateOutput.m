if Display; fprintf('Saving output to disk\n'); end

% grids
save([OutputDir,'grids'],'agrid','bgrid','ygrid','ydist','adelta','bdelta')

% initial steady state summary stats
fid = fopen(fullfile(OutputDir,'InitialSteadyStateParameters.txt'),'w');
fprintf(fid,'gam %.15g\n',gam);
fprintf(fid,'rho %.15g\n',rho);
fprintf(fid,'deathrate %.15g\n',deathrate);
fprintf(fid,'kappa0_w %.15g\n',kappa0_w);
fprintf(fid,'kappa1_w %.15g\n',kappa1_w);
fprintf(fid,'kappa2_w %.15g\n',kappa2_w);
fprintf(fid,'kappa3 %.15g\n',kappa3);
fprintf(fid,'kappa4_w %.15g\n',kappa4_w);
fprintf(fid,'kappa0_d %.15g\n',kappa0_d);
fprintf(fid,'kappa1_d %.15g\n',kappa1_d);
fprintf(fid,'kappa2_d %.15g\n',kappa2_d);
fprintf(fid,'kappa4_d %.15g\n',kappa4_d);
fprintf(fid,'corptax %.15g\n',corptax);

fprintf(fid,'ra %.15g\n',equmINITSS.ra);
fprintf(fid,'rb %.15g\n',equmINITSS.rb);
fprintf(fid,'rborr %.15g\n',equmINITSS.rborr);
fprintf(fid,'rcapital %.15g\n',equmINITSS.rcapital);
fprintf(fid,'wage %.15g\n',equmINITSS.wage);
fprintf(fid,'netwage %.15g\n',equmINITSS.netwage);
fprintf(fid,'bond %.15g\n',equmINITSS.bond);
fprintf(fid,'capital %.15g\n',equmINITSS.capital);
fprintf(fid,'equity %.15g\n',equmINITSS.equity);
fprintf(fid,'labor %.15g\n',equmINITSS.labor);
fprintf(fid,'output %.15g\n',equmINITSS.output);
fprintf(fid,'investment %.15g\n',equmINITSS.investment);
fprintf(fid,'govexp %.15g\n',equmINITSS.govexp);
fprintf(fid,'lumptransfer %.15g\n',equmINITSS.lumptransfer);
fprintf(fid,'labtax %.15g\n',equmINITSS.labtax);
fprintf(fid,'taxrev %.15g\n',equmINITSS.taxrev);
fprintf(fid,'govbond %.15g\n',equmINITSS.govbond);
fprintf(fid,'worldbond %.15g\n',equmINITSS.worldbond);
fprintf(fid,'fundbond %.15g\n',equmINITSS.fundbond);
fprintf(fid,'KYratio %.15g\n',equmINITSS.KYratio);
fprintf(fid,'KNratio %.15g\n',equmINITSS.KNratio);
fprintf(fid,'mc %.15g\n',equmINITSS.mc);
fprintf(fid,'rb %.15g\n',equmINITSS.rb);
fprintf(fid,'tfp %.15g\n',equmINITSS.tfp);
fprintf(fid,'pi %.15g\n',equmINITSS.pi);
fprintf(fid,'rnom %.15g\n',equmINITSS.rnom);
fprintf(fid,'priceadjust %.15g\n',equmINITSS.priceadjust);
fprintf(fid,'profit %.15g\n',equmINITSS.profit);
fprintf(fid,'dividend %.15g\n',equmINITSS.dividend);
fprintf(fid,'divrate %.15g\n',equmINITSS.divrate);

fprintf(fid,'borrwedge %.15g\n',equmINITSS.borrwedge);
fprintf(fid,'rho %.15g\n',equmINITSS.rho);
fprintf(fid,'fundlev %.15g\n',equmINITSS.fundlev);
fprintf(fid,'deprec %.15g\n',deprec);
	
fprintf(fid,'Ea %.15g\n',statsINITSS.Ea);
fprintf(fid,'Eb %.15g\n',statsINITSS.Eb);
fprintf(fid,'Ec %.15g\n',statsINITSS.Ec);
fprintf(fid,'Ehours %.15g\n',statsINITSS.Ehours);
fprintf(fid,'Elabor %.15g\n',statsINITSS.Elabor);
fprintf(fid,'Ed %.15g\n',statsINITSS.Ed);
fprintf(fid,'Ewage %.15g\n',statsINITSS.Ewage);
fprintf(fid,'Enetlabinc %.15g\n',statsINITSS.Enetlabinc);
fprintf(fid,'Egrosslabinc %.15g\n',statsINITSS.Egrosslabinc);
fprintf(fid,'Enetprofinc %.15g\n',statsINITSS.Enetprofinc);
fprintf(fid,'Egrossprofinc %.15g\n',statsINITSS.Egrossprofinc);
fprintf(fid,'Einc %.15g\n',statsINITSS.Einc);
fprintf(fid,'Enw %.15g\n',statsINITSS.Enw);
fprintf(fid,'FRACa0 %.15g\n',statsINITSS.FRACa0);
fprintf(fid,'FRACb0 %.15g\n',statsINITSS.FRACb0);
fprintf(fid,'FRACb0a0 %.15g\n',statsINITSS.FRACb0a0);
fprintf(fid,'FRACnw0 %.15g\n',statsINITSS.FRACnw0);
fprintf(fid,'FRACb0aP %.15g\n',statsINITSS.FRACb0aP);
fprintf(fid,'FRACa0close %.15g\n',statsINITSS.FRACa0close);
fprintf(fid,'FRACb0close %.15g\n',statsINITSS.FRACb0close);
fprintf(fid,'FRACb0a0close %.15g\n',statsINITSS.FRACb0a0close);
fprintf(fid,'FRACnw0close %.15g\n',statsINITSS.FRACnw0close);
fprintf(fid,'FRACbN %.15g\n',statsINITSS.FRACbN);
fprintf(fid,'EbN %.15g\n',statsINITSS.EbN);
fprintf(fid,'EbP %.15g\n',statsINITSS.EbP);
fprintf(fid,'Eadjcost %.15g\n',statsINITSS.Eadjcost);
fprintf(fid,'GINIa %.15g\n',statsINITSS.GINIa);
fprintf(fid,'GINIb %.15g\n',statsINITSS.GINIb);
fprintf(fid,'GINInw %.15g\n',statsINITSS.GINInw);
fprintf(fid,'GINIc %.15g\n',statsINITSS.GINIc);
fprintf(fid,'GINIinc %.15g\n',statsINITSS.GINIinc);
fprintf(fid,'Ec_bN %.15g\n',statsINITSS.Ec_bN);
fprintf(fid,'Ec_b0close %.15g\n',statsINITSS.Ec_b0close);
fprintf(fid,'Ec_b0far %.15g\n',statsINITSS.Ec_b0far);
	
fclose(fid);

% initial steady state distributions and policy functions
save([OutputDir,'solnINITSS'],'-struct','solnINITSS')
save([OutputDir,'cumINITSS'],'-struct','cumINITSS')
save([OutputDir,'equmINITSS'],'-struct','equmINITSS')
save([OutputDir,'statsINITSS'],'-struct','statsINITSS')