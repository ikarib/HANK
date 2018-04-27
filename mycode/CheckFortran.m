max(max(max(abs(agrid-load('agrid.txt')))))
max(max(max(abs(bgrid-reshape(load('bgrid.txt'),1,ngpb,1)))))
max(max(max(abs(ygrid-reshape(load('ygrid.txt'),1,1,ngpy)))))
max(max(max(abs(V-reshape(load('V.txt'),ngpa,ngpb,ngpy)))))
max(max(max(abs(c-reshape(load('c.txt'),ngpa,ngpb,ngpy)))))
max(max(max(abs(h-reshape(load('h.txt'),ngpa,ngpb,ngpy)))))
max(max(max(abs(s-reshape(load('s.txt'),ngpa,ngpb,ngpy)))))
max(max(max(abs(d-reshape(load('d.txt'),ngpa,ngpb,ngpy)))))
max(max(max(abs(u-reshape(load('u.txt'),ngpa,ngpb,ngpy)))))
% max(max(max(abs(bdot-reshape(load('bdot.txt'),ngpa,ngpb,ngpy)))))
max(max(max(abs(V-reshape(load('V.txt'),ngpa,ngpb,ngpy)))))
max(max(max(abs(Vnew-reshape(load('Vnew.txt'),ngpa,ngpb,ngpy)))))

% ACOOrow = load('ACOOrow.txt'); 
% ACOOcol = load('ACOOcol.txt');
% ACOOval = load('ACOOval.txt');
% ACOOdif = inf(ngpy,1);
% iba = reshape(1:nab,ngpb,ngpa)';
% for iy=1:ngpy
%     iab = ACOOrow(iy,:)~=0 & ACOOcol(iy,:)~=0;
%     Acoo{iy} = sparse(ACOOrow(iy,iab),ACOOcol(iy,iab),ACOOval(iy,iab),nab,nab,5*nab);
%     ACOOdif(iy)=max(max(abs(ACOO{iy}-Acoo{iy}(iba,iba))));
% end
% max(ACOOdif)

AUCOOrow = load('AUCOOrow.txt'); 
AUCOOcol = load('AUCOOcol.txt');
AUCOOval = load('AUCOOval.txt');
AUCOOdif = inf(ngpy,1);
iba = reshape(1:nab,ngpb,ngpa)';
for iy=1:ngpy
    iab = AUCOOrow(iy,:)~=0 & AUCOOcol(iy,:)~=0;
    AUCOO{iy} = sparse(AUCOOrow(iy,iab),AUCOOcol(iy,iab),AUCOOval(iy,iab),nab,nab,5*nab);
    AUCOOdif(iy)=max(max(abs(AUCOO{iy}(iba,iba)-AU{iy})));
end
max(AUCOOdif)


max(max(max(abs(gmat-reshape(permute(reshape(load('gmat.txt'),ngpb,ngpa,ngpy),[2 1 3]),nab,ngpy)))))
max(max(max(abs(ydist-reshape(load('ydist.txt'),1,1,ngpy)))))
max(max(max(abs(adelta-reshape(load('adelta.txt'),ngpa,1,1)))))
max(max(max(abs(bdelta-reshape(load('bdelta.txt'),1,ngpb,1)))))
max(max(abs(squeeze(gamarg)-load('gamarg.txt'))))
max(max(abs(squeeze(gbmarg)-load('gbmarg.txt'))))
max(max(abs(gabmarg-load('gabmarg.txt'))))
max(max(abs(gabcum-load('gabcum.txt'))))
max(max(max(abs(lmargdist-permute(reshape(load('lmargdist.txt'),ngpy,ngpb,ngpa),[3 2 1])))))
max(abs(lnw_a-load('lnw_a.txt')'))
max(abs(lnw_b-load('lnw_b.txt')'))
max(abs(lnw_c-load('lnw_c.txt')'))
max(abs(lnw_h-load('lnw_h.txt')'))
max(abs(lnw_inc-load('lnw_inc.txt')'))
max(abs(lnwmargcum-load('lnwmargcum.txt')'))
max(abs(lcgrid-load('lcgrid.txt')'))
% max(abs(lcmargdist-load('lcmargdist.txt')'))
max(abs(lcdelta-load('lcdelta.txt')'))
% max(abs(lcmargcum-load('lcmargcum.txt')'))
max(abs(lincgrid-load('lincgrid.txt')'))
max(abs(lincdelta-load('lincdelta.txt')'))
max(abs(linc_a-load('linc_a.txt')'))
max(abs(linc_b-load('linc_b.txt')'))
max(abs(linc_c-load('linc_c.txt')'))
max(abs(linc_h-load('linc_h.txt')'))
max(abs(linc_nw-load('linc_nw.txt')'))
max(abs(lincmargcum-load('lincmargcum.txt')'))
max(max(abs([Ea_nwQ;Eb_nwQ;Ec_nwQ;Einc_nwQ;Ea_incQ;Eb_incQ;Ec_incQ;Einc_incQ]-load('E.txt'))))

ccum4_=permute(reshape(load('lccumvec.txt'),ngpb,ngpa,ngpy),[2 1 3]);
dcum4_=permute(reshape(load('ldcumvec.txt'),ngpb,ngpa,ngpy),[2 1 3]);
max(max(max(abs(ccum4./ccum4_-1))))
max(max(max(abs(dcum4./dcum4_-1))))

% DiscountedMPC
load('lsubeff1ass.txt')
max(max(abs(subeff1ass-lsubeff1ass(reshape(1:nab,ngpb,ngpa)',:))))
load('lsubeff2ass.txt')
max(max(abs(subeff2ass-lsubeff2ass(reshape(1:nab,ngpb,ngpa)',:))))
load('lwealtheff1ass.txt')
max(max(abs(wealtheff1ass-lwealtheff1ass(reshape(1:nab,ngpb,ngpa)',:))))
load('lwealtheff2ass.txt')
max(max(abs(wealtheff2ass-lwealtheff2ass(reshape(1:nab,ngpb,ngpa)',:))))


max(abs(deltatransvec-load('deltatransvec.txt')'))

%% check IRF
for ipe = 0:15
    if ipe; OutputFileIRF = [OutputDir,'IRF_',IRFDir,'/PE',num2str(ipe)];
      else; OutputFileIRF = [OutputDir,'IRF_',IRFDir,'/NOFS']; end
    S = load(OutputFileIRF);
    for f=fields(S.equmSTICKY)'
        f{1}
        temp = load([OutputFileIRF,'/STICKY/',f{1},'.txt']);
        max(max(abs(S.equmSTICKY.(f{1})-temp')))
    end
%     S.statsSTICKY
%     S.solnSTICKY
end
