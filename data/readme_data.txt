data.mat (MATLAB data file)

Contains: BZAL, BZPL, EPAL, EPPL, SCAL, SCPL
SC is participant 1
EP is participant 2
BZ is participant 3
AL = active learning condition
PL = passive learning condition
The data type of these are matlab struct


Relevant struct fields:

RevealPosX, RevealPosY:
In the AL data: participant's revealing positions on each trial
In the PL data: computer generated revealing position in each trial
Rows index trials; columns index revealing number.

AnswerChoice: 
participant's categorization

AnswerReal:
the true category of the underlying image

ImageID:
identifiers of the image used for each trial.
col 1: 1->patch; 2->stripy
col 2: inverse length scale used in GP for x: 20->patch; 6->stripy horizontal; 30->stripy vertical
col 3: inverse length scale used in GP for y
col 4: an image id
The actual images are in the data on <https://archive.org/download/images_201907>.
This is a big MATLAB datafile (~13GB).
These images are grouped into 6 structs named DIMBZAL, DIMBZPL, DIMEPAL, DIMEPPL, DIMSCAL, DIMSCPL to correspond to the BZAL, BZPL, EPAL, EPPL, SCAL, SCPL structs.
The ith row in a struct is the underlying image (770x770 pixels) used in the ith trial.
For example, the 10the row in DIMSCAL is the underlying image used for the 10th trial in SCAL.

MaxRevealingTrial:
The maximum number of revealings allowed for each trial

RevealType (only in the PL data):
The type of revealing generated by the computer
0 = random revealings
1 = using BAS when the underlying image is patch
2 = using BAS when the underlying image is stripy horizontal
3 = using BAS when the underlying image is stripy vertical