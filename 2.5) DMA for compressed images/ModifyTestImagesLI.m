clear
clc
load('SPTest15Gsmall')
load('SPTest20Gsmall')
load('SWFreshOnlyGsmall')
load('SWOnlyGsmall')

SPTest15Gdark=SPTest15Gsmall.*1.1;
SPTest20Gdark=SPTest20Gsmall.*1.1;
SWFreshOnlyGdark=SWFreshOnlyGsmall.*1.1;
SWOnlyGdark=SWOnlyGsmall.*1.1;
save('SPTest15Gdark','SPTest15Gdark')
save('SPTest20Gdark','SPTest20Gdark')
save('SWFreshOnlyGdark','SWFreshOnlyGdark')
save('SWOnlyGdark','SWOnlyGdark')