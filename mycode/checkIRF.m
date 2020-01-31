function checkIRF(irf,OutputFileIRF)
equm = struct;
for f = setdiff(fields(irf.equmSTICKY),{'labtax'},'stable')'
    equm.(f{1}) = load([OutputFileIRF,f{1},'.txt'])';
    dif = max(abs(irf.equmSTICKY.(f{1})-equm.(f{1})));
    fprintf('%s = %g\n',f{1},dif)
    if dif>1e-4
        plot([irf.equmSTICKY.(f{1}); equm.(f{1})]'),title(f{1}),legend('mine','theirs')
        error('%s differs',f{1})
    end
end
stats = struct;
for f = setdiff(fields(irf.statsSTICKY),{'Enetprofinc','Egrossprofinc','FRACb0a0close'},'stable')'
    stats.(f{1}) = load([OutputFileIRF,f{1},'.txt']);
    dif = max(max(abs(vertcat(irf.statsSTICKY(:).(f{1}))-stats.(f{1}))));
    fprintf('%s = %g\n',f{1},dif)
%     if dif>1e-4
%         figure,plot([vertcat(irf.statsSTICKY(:).(f{1})) stats.(f{1})]),title(f{1}),legend('mine','theirs')
%         error('%s differs',f{1})
%     end
end
soln = struct;
for f = setdiff(fields(irf.solnSTICKY),{'AU','gamarg','gbmarg'},'stable')'
    switch f{1}
        case 'c'; name='con';
        case 'd'; name='dep';
        case 'gmat'; name='gjoint';
        case 'h'; name='hour';
        case 's'; name='sav';
        otherwise; name=f{1};
    end
    soln.(f{1})=[];
    for i=1:33
        soln.(f{1})(:,:,i) = load([OutputFileIRF,name,'_T1_y',num2str(i),'.txt']);
    end
    dif = max(max(max(abs(irf.solnSTICKY.(f{1}){1}-soln.(f{1})))));
    fprintf('%s = %g\n',name,dif)
%     if dif>1e-4
%         surface(irf.solnSTICKY.(f{1}){1}),title(f{1}),legend('mine','theirs')
%         soln.(f{1})
%         error('%s differs',f{1})
%     end
end
