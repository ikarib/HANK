clear
SetParameters
% distcomp.feature('LocalUseMpiexec',false);parpool('local',ngpy);
Procedures
if CalibrateCostFunction; Calibration; end
InitialSteadyState
SaveSteadyStateOutput
FinalSteadyState
if DoImpulseResponses; ImpulseResponses; end