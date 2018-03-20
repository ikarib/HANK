PROGRAM Main

USE Parameters
USE Globals

CALL SetParameters
CALL MakeGuess
CALL Estimation
CALL CombinedProcess
CALL SaveOutput

END PROGRAM Main