//
//  main.c
//  HCV Model
//
//  Created by Neil Bretana on 15/08/2014.
//  Copyright (c) 2014 Neil Bretana. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "individual.h"

#define TCOLCTR 2   //Array of locations containing:
                    //[row][0] IDU+;
                    //[row][1] IDU-;
#define COLCTR 12  //Array of locations containing:
                    //[row][0] IDU+; HCV+; ATSI
                    //[row][1] IDU+; HCV+; NON-ATSI
                    //[row][2] IDU+; HCV-; ATSI; PREV. EXPOSED
                    //[row][3] IDU+; HCV-; ATSI; SUSCEPTIBLE
                    //[row][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED
                    //[row][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE
                    //[row][6] IDU-; HCV+; ATSI
                    //[row][7] IDU-; HCV+; NON-ATSI
                    //[row][8] IDU-; HCV-; ATSI; PREV. EXPOSED
                    //[row][9] IDU-; HCV-; ATSI; SUSCEPTIBLE
                    //[row][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED
                    //[row][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE
#define ROWPRIS 4 //# prisons +1 community
#define NEVENTS 27 //number of events
#define RGROUPS 13 //Risk groups: 0: non-injecting; 1: injecting less than daily; opioid; not sharing 2: injecting less than daily; opioid; sharing less daily 3: injecting less than daily; opioid; sharing daily more 4: injecting less than daily; non-opioid; no sharing 5: injecting less than daily; non-opioid; sharing less daily 6: injecting less than daily; non-opioid; sharing daily or more 7: ...

const gsl_rng_type *T;

int generateRand();
int draw_multinom(gsl_rng **r, int nEvents, double probsInput[]);
int distributePop(int input);
void traverse (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countHCVpris (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int countHCVCom (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count0(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count1 (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count2(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count3(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count4 (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count5(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count6(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count7 (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count8(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count9(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count10 (sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count11(sIndiv **pHeadCopy, sIndiv **pTailCopy);
int count12(sIndiv **pHeadCopy, sIndiv **pTailCopy);
void newIndiv (int *idGlobal, sIndiv **pHeadCopy, sIndiv **pTailCopy, gsl_rng **r, int nEvents, int (*pLocArray2)[ROWPRIS][COLCTR], int prison, int category);
double probNtrDeath(gsl_rng **r, sIndiv **pTarget); //0 
double probStopInj(gsl_rng **r, sIndiv **pTarget); //1 
double probInjLdOpNoSh(gsl_rng **r, sIndiv **pTarget); //2
double probInjLdOpShLd(gsl_rng **r, sIndiv **pTarget); //3
double probInjLdOpShDm(gsl_rng **r, sIndiv **pTarget); //4
double probInjLdNoNoSh(gsl_rng **r, sIndiv **pTarget); //5 
double probInjLdNoShLd(gsl_rng **r, sIndiv **pTarget); //6
double probInjLdNoShDm(gsl_rng **r, sIndiv **pTarget); //7
double probInjDmOpNoSh(gsl_rng **r, sIndiv **pTarget); //8
double probInjDmOpShLd(gsl_rng **r, sIndiv **pTarget); //9
double probInjDmOpShDm(gsl_rng **r, sIndiv **pTarget); //10
double probInjDmNoNoSh(gsl_rng **r, sIndiv **pTarget); //11
double probInjDmNoShLd(gsl_rng **r, sIndiv **pTarget); //12
double probInjDmNoShDm(gsl_rng **r, sIndiv **pTarget); //13 
double probMoveCom(gsl_rng **r, sIndiv **pTarget); //14
double probMoveP1(gsl_rng **r, sIndiv **pTarget); //15
double probMoveP2(gsl_rng **r, sIndiv **pTarget); //16
double probMoveP3(gsl_rng **r, sIndiv **pTarget); //17
double probInfectIndiv(sIndiv **pTarget); //18
double probNtrClear(gsl_rng **r, sIndiv **pTarget, int currDay); //19
double probProgress(gsl_rng **r, sIndiv **pTarget); //20
double probHCVdeath(gsl_rng **r, sIndiv **pTarget); //21
double probOSTstart(gsl_rng **r, sIndiv **pTarget); //22
double probOSTstop(gsl_rng **r, sIndiv **pTarget); //23
double probDAAstart(gsl_rng **r, sIndiv **pTarget); //24
double probDAAstop(gsl_rng **r, sIndiv **pTarget); //25
double probDAAclear(gsl_rng **r, sIndiv **pTarget, int currDay); //26

void age(sIndiv **pTargetCopy);
void removeIndiv(sIndiv **pTarget, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR]);
void changeRisk(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int newRisk);
void moveLocation(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int newLoc);
int infect(gsl_rng **r, int nEvents, sIndiv **pTargetCopy, sIndiv **pHeadCopy, sIndiv **pTailCopy, int (*pLocArray2)[ROWPRIS][COLCTR], int currDay);
void clearHCV(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR]);
void progressHCV(sIndiv **pTargetCopy);
int deathHCV(sIndiv **pTargetCopy);
void startOST(sIndiv **pTargetCopy, int currDay);
void stopOST(sIndiv **pTargetCopy);
void startDAA(sIndiv **pTargetCopy, int currDay);
void stopDAA(sIndiv **pTargetCopy);
void clearDAA(sIndiv **pTargetCopy, int (*pLocArray2)[ROWPRIS][COLCTR]);

int main(int argc, const char * argv[])
{
    //int newInfected, newCleared; //new cases per week/month
    FILE *fp, *fw, *fr, *fm;
    int idG=0, totalInfectedCommunity=0, daysSim, i, j, k, l, normCount, nShare, totalIndivDAA=0, totalIndivOST=0, totalClrN=0, totalClrD=0, fi, ni, newCases=0, success, currDay=1, totalPrisonPop=0, dead=0, deadHCV=0, deadHCVCom=0, deadHCVPris=0, prisonHCV=0, communityHCV=0, releasedHCV=0, releasedHCVCom=0, releasedHCVPris=0, released=0, eIndex=999, eInf=999, deadFlag=0, eFlag, mFlag=0, dFlag=0, cnFlag=0, cdFlag=0, cEvents, cIDUp=0, cIDUn=0, NSWpop=6774000, inPop; //total prison population, number of days to simulate; NSW pop=2005
    int E0=0, E1=0, E2=0, E3=0, E4=0, E5=0, E6=0, E7=0, E8=0, E9=0, E10=0, E11=0, E12=0, E13=0, E14=0, E15=0, E16=0, E17=0, E18=0, E19=0, E20=0, E21=0, E22=0, E23=0, E24=0, E25=0, E26=0, R0=0, R1=0, R2=0, R3=0, R4=0, R5=0, R6=0, R7=0, R8=0, R9=0, R10=0, R11=0, R12=0;//Event counters
    gsl_rng *r;
    //int locArray[ROWPRIS][COLCTR]; //Array of locations containing: //[row][0] IDU+ //[row][1] IDU-
    //int (*pLocArray)[ROWPRIS][COLCTR]=&locArray;
    //double probLocMat[ROWPRIS][ROWPRIS]; //[row][col] for movement between locations
    //double probRiskMat[RGROUPS][RGROUPS]; //[row][col] for movement between locations
    double riskNormDenominator, movesNormDenominator, normDenominator, probEvents[NEVENTS], probInf[2];
    //arrays for events
	double eventProbs[27], eventDraw[2]; 
	int eventDecision[27], eventOrder[27], outArray[2];
	int eCtr, iFlag, iCrawl, iCount; 
    int tempLocArray[ROWPRIS][TCOLCTR]; //Array of locations containing: //[row][0] IDU+ //[row][1] IDU-
    int (*pTempLocArray)[ROWPRIS][TCOLCTR]=&tempLocArray;
    int locArray[ROWPRIS][COLCTR]; //Array of locations containing: [row][0] IDU+; HCV+; ATSI, [row][1] IDU+; HCV+; NON-ATSI, [row][2] IDU+; HCV-; ATSI; PREV. EXPOSED, [row][3] IDU+; HCV-; ATSI; SUSCEPTIBLE, [row][4] IDU+; HCV-; NON-ATSI; PREV. EXPOSED, [row][5] IDU+; HCV-; NON-ATSI; SUSCEPTIBLE, [row][6] IDU-; HCV+; ATSI, [row][7] IDU-; HCV+; NON-ATSI, [row][8] IDU-; HCV-; ATSI; PREV. EXPOSED, [row][9] IDU-; HCV-; ATSI; SUSCEPTIBLE, [row][10] IDU-; HCV-; NON-ATSI; PREV. EXPOSED, [row][11] IDU-; HCV-; NON-ATSI; SUSCEPTIBLE
    int (*pLocArray)[ROWPRIS][COLCTR]=&locArray;
    double fd;
    float rProp, lb, ub; //For population
    char option1[2]={"b"}, option2[2]={"d"}; //basic or detailed input
    
    //Create a generator chosen by the environment variable gsl_rng_type *
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    gsl_rng_set(r, (unsigned long)time(NULL));
    
    //printf("Days:");
    //scanf("%d", &daysSim);
    
    srand((unsigned int)time(NULL));
    
    //Set linked list parameters
    sIndiv *nHead, *nTail, *target;
    nHead=NULL;
    nTail=NULL;
    
    //printf("int: %f\n", gsl_ran_exponential(r, 0.0004593));
    //printf("int: %f\n", gsl_ran_beta(r, 1.28527, 2797.038));
    //printf("int: %f\n", gsl_ran_exponential(r, 0.0004593));
    //printf("int: %f\n", gsl_ran_beta(r, 1.28527, 2797.038));
    //printf("int: %f\n", gsl_ran_exponential(r, 0.0004593));
    //printf("int: %f\n", gsl_ran_beta(r, 1.28527, 2797.038));
    //exit(0);
    
    for(i=0;i<ROWPRIS;i++){ //read through prisons
        for(j=0;j<COLCTR;j++){
            locArray[i][j]=0;
        }
    }
    
    if ( argc != 4) // There should be one option
    {
        printf( "Missing arguments\n");
        printf( "Input file name, Output file name, Run option\n");
        exit(0);
    }
    
    printf( "Input read: done");
    //printf( "Input: %s\n", argv[1] );
    //printf( "Output: %s\n", argv[2] );
    //printf( "Option: %s\n", argv[3] ); //basic or detailed; use d
    //exit(0);
     
    //use files
    //fp=fopen("/Users/neilbretana/Dropbox/HCV Model/HCV New Model/HCV Model/input.txt", "r");
    fp=fopen(argv[1], "r");
    //fw=fopen("/Users/neilbretana/Dropbox/HCV Model/HCV New Model/HCV Model/output.xls", "w");
    fw=fopen(argv[2], "w");
    
    //argv[3]="d";

    if(fw==NULL){
        printf("Cannot open file!\n");
        exit(1);
    }
    fprintf(fw, "TIME\tC IDU+HCV+ATSI+\tC IDU+HCV+ATSI-\tC IDU+HCV-ATSI+PE\tC IDU+HCV-ATSI+S\tC IDU+HCV-ATSI-PE\tC IDU+HCV-ATSI-S\tC IDU-HCV+ATSI+\tC IDU-HCV+ATSI-\tC IDU-HCV-ATSI+PE\tC IDU-HCV-ATSI+S\tC IDU-HCV-ATSI-PE\tC IDU-HCV-ATSI-S\tP1 IDU+HCV+ATSI+\tP1 IDU+HCV+ATSI-\tP1 IDU+HCV-ATSI+PE\tP1 IDU+HCV-ATSI+S\tP1 IDU+HCV-ATSI-PE\tP1 IDU+HCV-ATSI-S\tP1 IDU-HCV+ATSI+\tP1 IDU-HCV+ATSI-\tP1 IDU-HCV-ATSI+PE\tP1 IDU-HCV-ATSI+S\tP1 IDU-HCV-ATSI-PE\tP1 IDU-HCV-ATSI-S\tP2 IDU+HCV+ATSI+\tP2 IDU+HCV+ATSI-\tP2 IDU+HCV-ATSI+PE\tP2 IDU+HCV-ATSI+S\tP2 IDU+HCV-ATSI-PE\tP2 IDU+HCV-ATSI-S\tP2 IDU-HCV+ATSI+\tP2 IDU-HCV+ATSI-\tP2 IDU-HCV-ATSI+PE\tP2 IDU-HCV-ATSI+S\tP2 IDU-HCV-ATSI-PE\tP2 IDU-HCV-ATSI-S\tP3 IDU+HCV+ATSI+\tP3 IDU+HCV+ATSI-\tP3 IDU+HCV-ATSI+PE\tP3 IDU+HCV-ATSI+S\tP3 IDU+HCV-ATSI-PE\tP3 IDU+HCV-ATSI-S\tP3 IDU-HCV+ATSI+\tP3 IDU-HCV+ATSI-\tP3 IDU-HCV-ATSI+PE\tP3 IDU-HCV-ATSI+S\tP3 IDU-HCV-ATSI-PE\tP3 IDU-HCV-ATSI-S\tNEW CASES\tPRISON HCV\tCOMMUNITY HCV\tRELEASED\tRELEASED HCV\tRELEASED HCV COM\tRELEASED HCV PRIS\tCLEAR HCV\tDEATHS\tHCV DEATHS\tHCV COM DEATH\tHCV PRISON DEATH\tTOTAL PRISON POP\tE0\tE1\tE2\tE3\tE4\tE5\tE6\tE7\tE8\tE9\tE10\tE11\tE12\tE13\tE14\tE15\tE16\tE17\tE18\tE19\tE20\tE21\tE22\tE23\tE24\tE25\tE26\tR0\tR1\tR2\tR3\tR4\tR5\tR6\tR7\tR8\tR9\tR10\tR11\tR12\n");
    
    if(fp==NULL){
        printf("Cannot open file\n");
        exit(0);
    }else{
        //check if correct file format
        
        //if file correct
        fscanf(fp, "%d,", &daysSim);
        //printf("Days: %d\n", daysSim);
        
        //option a
        //    for(i=0;i<4;i++){ //read through all 4 populations
        //        fscanf(fp, "%d,", &fi); //read input
        //        printf("No. of people in location %d group %d: %d\n", i, j, fi); //print input
        //    }
        
        if(strcmp(argv[3],option1)==0){ //basic populations USE COLCTR to read
            //fp=fopen("/Users/neilbretana/Dropbox/HCV Model/HCV New Model/HCV Model/inputII.txt", "r");
            for(i=0;i<ROWPRIS;i++){ //read through prisons
                for(j=0;j<TCOLCTR;j++){ //0 and 1, if 0 IDU+, 1 IDU-
                    fscanf(fp, "%d,", &fi);
                    //printf("No. of people in input location %d group %d: %d\n", i, j, fi); //
                    for(k=0;k<COLCTR-6;k++){ //subpopulation array index
                        if(fi>0){
                            nShare=distributePop(fi);//randomise
                            fi=fi-nShare;
                            if(i!=0){ //prisons 1-3
                                //printf("No. of people in location %d group %d: %d\n", i, k, nShare); //
                                for(l=0; l<nShare; l++){
                                    //newIndiv(&idG, &nHead, &nTail, &r, ROWPRIS, probLocMat, pLocArray, i, k);
                                }
                            }else{ //community
                                if(j==0){ //IDU+
                                    locArray[i][k]=nShare;
                                    //printf("No. of people in location %d group %d: %d\n", i, k, locArray[i][k]);
                                }else if(j==1){ //IDU-
                                    locArray[i][k+6]=nShare;
                                    //printf("No. of people in location %d group %d: %d\n", i, k+6, locArray[i][k+6]);
                                }
                            }
                        }
                    }
                }
            }
        }else if(strcmp(argv[3],option2)==0){ //detailed populations USE COLCTR2 to read
            //fp=fopen("/Users/neilbretana/Dropbox/HCV Model/HCV New Model/HCV Model/input.txt", "r");
            for(i=0;i<ROWPRIS;i++){ //read through prisons
                for(j=0;j<COLCTR;j++){ //subpopulation array index
                    fscanf(fp, "%d,", &fi);
                    if(i!=0){ //prisons 1-3
                        printf("No. of people in location %d group %d: %d\n", i, j, fi); //
                        for(l=0; l<fi; l++){
                            newIndiv(&idG, &nHead, &nTail, &r, ROWPRIS, pLocArray, i, j); //creating individual agents
                        }
                    }else{ //community
                        locArray[i][j]=fi;
                        printf("No. of people in location %d group %d: %d\n", i, j, locArray[i][j]);
                    }
                }
            }

        }
        
        printf( "Initialise populations: done \n");
        fclose(fp);
        
        //else if incorrect
            //exit
    }
    
    //simulation starts here
    //loop depending on number of days being simulated
    while(currDay!=daysSim+1){
        //printf("DAY %d\n", currDay);
        
        //Move non-individual subpopulations from community to prisons
        //i=0; //community
        //function creates new individuals from the community, choose from 0-23
        //lb=0;//Min influx from DCS
        //ub=23;//Max influx from DCS
        //lb=1;//Min Max from influx frquencies Influx0515 file MONTHLY though
        //ub=703;
        //inPop=gsl_ran_flat(r, lb, ub); // uniform distribution between min and max (min= 4.793333e-06; max=6.74e-06) <- decimal percent
        inPop=gsl_ran_exponential(r, 2.67); //See influxtrnds.xls 2.67 is the mean influx in days based on data from DCS
        //printf("rProp %f\n", rProp);
        //inPop=NSWpop*rProp; // get proportion from NSW population
        //printf("inPop %d\n", inPop);
        for(k=0;k<inPop;k++){
            newIndiv(&idG, &nHead, &nTail, &r, ROWPRIS, pLocArray, 99, 99);//new prisoner to be moved from com to a random prison based on **
        }
       
        /*
        i=0; //community
            for(j=0;j<COLCTR;j++){ //read through non-individual subpoulations
                if(locArray[i][j]>0){ //make sure current loc has enough people
                    //ni=1; //Determine n number of non-injectors to be moved to this location nloc; ; has to be 0-n with n<=total current numbr in that bucket
                    lb=0;//Min influx from DCS
                    ub=23;//Max influx from DCS
                    inPop=gsl_ran_flat(r, lb, ub); // uniform distribution between min and max (min= 4.793333e-06; max=6.74e-06) <- decimal percent
                    //printf("rProp %f\n", rProp);
                    //inPop=NSWpop*rProp; // get proportion from NSW population
                    //printf("inPop %d\n", inPop);
                    for(k=0;k<inPop;k++){
                        newIndiv(&idG, &nHead, &nTail, &r, ROWPRIS, pLocArray, 99, j);//new prisoner to be moved from com to a random prison based on **
                    }
                }
            }
         */
        
        //Traverse one by one the individuals in each prison
        //Minimum
        target=nHead; //Target points to the individual at the beginning of the list
        while(target!=NULL){
            //age if 365
            printf("Modifying individal \n");
            
			if(currDay%365==0){
                age(&target);
            }
            
            //Fill up probability array of events
            if(target->group>=0&&target->group<=5){ //if at-risk, set probabilities accordingly
                probEvents[0]=probNtrDeath(&r, &target);//*natural death rate
                probEvents[1]=probStopInj(&r, &target);//*rate of stopping from injecting
                probEvents[2]=probInjLdOpNoSh(&r, &target);//*rate of injecting less daily, opioid, no sharing
                probEvents[3]=probInjLdOpShLd(&r, &target);//*rate of injecting less daily, opioid, sharing less than daily
                probEvents[4]=0.0;//probInjLdOpShDm(&r, &target);//*rate of injecting less daily, opioid, sharing daily or more
                probEvents[5]=probInjLdNoNoSh(&r, &target);//*rate of injecting less daily, non-opioid, no sharing
                probEvents[6]=probInjLdNoShLd(&r, &target);//*rate of injecting less daily, non-opioid, sharing less than daily
                probEvents[7]=0.0;//probInjLdNoShDm(&r, &target);//*rate of injecting less daily, non-opioid, sharing daily or more
                probEvents[8]=probInjDmOpNoSh(&r, &target);//*rate of injecting daily or more, opioid, no sharing
                probEvents[9]=probInjDmOpShLd(&r, &target);//*rate of injecting daily or more, opioid, sharing less than daily
                probEvents[10]=0.0;//probInjDmOpShDm(&r, &target);//*rate of injecting daily or more, opioid, sharing daily or more
                probEvents[11]=probInjDmNoNoSh(&r, &target);//*rate of injecting daily or more, non-opioid, no sharing
                probEvents[12]=probInjDmNoShLd(&r, &target);//*rate of injecting daily or more, non-opioid, sharing less than daily
                probEvents[13]=0.0;//probInjDmNoShDm(&r, &target);//*rate of injecting daily or more, non-opioid, sharing daily or more
                probEvents[14]=probMoveCom(&r, &target);//*rate of release to community
                probEvents[15]=probMoveP1(&r, &target);//*rate of movement to Prison 1
                probEvents[16]=probMoveP2(&r, &target);//*rate of movement to Prison 2
                probEvents[17]=probMoveP3(&r, &target);//*rate of movement to Prison 3
                probEvents[18]=probInfectIndiv(&target);//*rate of infecting susceptible individuals //PROBABILITY OF SHARING
                probEvents[19]=probNtrClear(&r, &target, currDay);//*rate of natural clearance
                probEvents[20]=probProgress(&r, &target);//*rate of HCV progression
                probEvents[21]=probHCVdeath(&r, &target);//*rate of HCV-related death
                probEvents[22]=0.0;//probOSTstart(&r, &target);//*rate of commencing OST
                probEvents[23]=0.0;//probOSTstop(&r, &target);//*rate of stopping OST
                probEvents[24]=0.0;//probDAAstart(&r, &target);//*rate of commencing DAA treatment
                probEvents[25]=0.0;//probDAAstop(&r, &target);//*rate of stopping DAA treatment
                probEvents[26]=0.0;//probDAAclear(&r, &target, currDay);//*rate of clearance via DAA
            }else{ //otherwise, set at-risk-related events to 0
                probEvents[0]=probNtrDeath(&r, &target);//*natural death rate
                if(target->group==9||target->group==11){ //target->group>=6&&target->group<=11 all non-IDU ever injected //target->group==9||target->group==11 only nonIDU never exposed to HCV never injected
                    probEvents[1]=probStopInj(&r, &target);//*rate of stopping from injecting
                    probEvents[2]=probInjLdOpNoSh(&r, &target)/5;//*rate of injecting less daily, opioid, no sharing //0.0
                    probEvents[3]=probInjLdOpShLd(&r, &target)/5;//*rate of injecting less daily, opioid, sharing less than daily
                    probEvents[4]=0.0;//probInjLdOpShDm(&r, &target)/3;//*rate of injecting less daily, opioid, sharing daily or more
                    probEvents[5]=probInjLdNoNoSh(&r, &target)/5;//*rate of injecting less daily, non-opioid, no sharing
                    probEvents[6]=probInjLdNoShLd(&r, &target)/5;//*rate of injecting less daily, non-opioid, sharing less than daily
                    probEvents[7]=0.0;//probInjLdNoShDm(&r, &target)/3;//*rate of injecting less daily, non-opioid, sharing daily or more
                    probEvents[8]=probInjDmOpNoSh(&r, &target)/5;//*rate of injecting daily or more, opioid, no sharing
                    probEvents[9]=probInjDmOpShLd(&r, &target)/5;//*rate of injecting daily or more, opioid, sharing less than daily
                    probEvents[10]=0.0;//probInjDmOpShDm(&r, &target)/3;//*rate of injecting daily or more, opioid, sharing daily or more
                    probEvents[11]=probInjDmNoNoSh(&r, &target)/5;//*rate of injecting daily or more, non-opioid, no sharing
                    probEvents[12]=probInjDmNoShLd(&r, &target)/5;//*rate of injecting daily or more, non-opioid, sharing less than daily
                    probEvents[13]=0.0;//probInjDmNoShDm(&r, &target)/3;//*rate of injecting daily or more, non-opioid, sharing daily or more
                }else{
                    probEvents[1]=probStopInj(&r, &target);//*rate of stopping from injecting
                    probEvents[2]=probInjLdOpNoSh(&r, &target)/5;//*rate of injecting less daily, opioid, no sharing
                    probEvents[3]=probInjLdOpShLd(&r, &target)/5;//*rate of injecting less daily, opioid, sharing less than daily
                    probEvents[4]=0.0;//probInjLdOpShDm(&r, &target);//*rate of injecting less daily, opioid, sharing daily or more
                    probEvents[5]=probInjLdNoNoSh(&r, &target)/5;//*rate of injecting less daily, non-opioid, no sharing
                    probEvents[6]=probInjLdNoShLd(&r, &target)/5;//*rate of injecting less daily, non-opioid, sharing less than daily
                    probEvents[7]=0.0;//probInjLdNoShDm(&r, &target);//*rate of injecting less daily, non-opioid, sharing daily or more
                    probEvents[8]=probInjDmOpNoSh(&r, &target)/5;//*rate of injecting daily or more, opioid, no sharing
                    probEvents[9]=probInjDmOpShLd(&r, &target)/5;//*rate of injecting daily or more, opioid, sharing less than daily
                    probEvents[10]=0.0;//probInjDmOpShDm(&r, &target);//*rate of injecting daily or more, opioid, sharing daily or more
                    probEvents[11]=probInjDmNoNoSh(&r, &target)/5;//*rate of injecting daily or more, non-opioid, no sharing
                    probEvents[12]=probInjDmNoShLd(&r, &target)/5;//*rate of injecting daily or more, non-opioid, sharing less than daily
                    probEvents[13]=0.0;//probInjDmNoShDm(&r, &target);//*rate of injecting daily or more, non-opioid, sharing daily or more
                }
                probEvents[14]=probMoveCom(&r, &target);//*rate of release to community
                probEvents[15]=probMoveP1(&r, &target);//*rate of movement to Prison 1
                probEvents[16]=probMoveP2(&r, &target);//*rate of movement to Prison 2
                probEvents[17]=probMoveP3(&r, &target);//*rate of movement to Prison 3
                probEvents[18]=0.0;//*rate of infecting susceptible individuals
                probEvents[19]=0.0;//*rate of natural clearance
                probEvents[20]=0.0;//*rate of HCV progression
                probEvents[21]=0.0;//*rate of HCV-related death
                probEvents[22]=0.0;//*rate of commencing OST
                probEvents[23]=0.0;//*rate of stopping OST
                probEvents[24]=0.0;//*rate of commencing DAA treatment
                probEvents[25]=0.0;//*rate of stopping DAA treatment
                probEvents[26]=0.0;//*rate of clearance via DAA
            }
            
            //NORMALIZE
            //Normalise each risk change first over all risk change
            //change risk behaviour indeces: 1, 2, 3, 5, 6, 8, 9, 11, 12
            //riskNormDenominator=probEvents[1]+probEvents[2]+probEvents[3]+probEvents[5]+probEvents[6]+probEvents[8]+probEvents[9]+probEvents[11]+probEvents[12];
            //probEvents[1]=probEvents[1]/riskNormDenominator;
            //probEvents[2]=probEvents[2]/riskNormDenominator;
            //probEvents[3]=probEvents[3]/riskNormDenominator;
            //probEvents[5]=probEvents[5]/riskNormDenominator;
            //probEvents[6]=probEvents[6]/riskNormDenominator;
            //probEvents[8]=probEvents[8]/riskNormDenominator;
            //probEvents[9]=probEvents[9]/riskNormDenominator;
            //probEvents[11]=probEvents[11]/riskNormDenominator;
            //probEvents[12]=probEvents[12]/riskNormDenominator;
            
            normDenominator=probEvents[0]+probEvents[1]+probEvents[2]+probEvents[3]+probEvents[4]+probEvents[5]+probEvents[6]+probEvents[7]+probEvents[8]+probEvents[9]+probEvents[10]+probEvents[11]+probEvents[12]+probEvents[13]+probEvents[14]+probEvents[15]+probEvents[16]+probEvents[17]+probEvents[18]+probEvents[19]+probEvents[20]+probEvents[21]+probEvents[22]+probEvents[23]+probEvents[24]+probEvents[25]+probEvents[26];
            
            for(normCount=0;normCount<27;normCount++){
                probEvents[normCount]=probEvents[normCount]/normDenominator;
                //printf("probEvents %f\n", probEvents[normCount]);
            }
            
            //1. fill up eventDraw 1 or 2; will a particular event occur or not? 
            eventDraw[0]=0; eventDraw[1]=0; outArray[0]=0; outArray[1]=0;
			for(eCtr=0; eCtr<27; eCtr++){
            	eventDecision[eCtr]=0;
            	eventOrder[eCtr]=eCtr;
			}
			printf("event draw array prepared \n");
            
            //2. Fill-up 1/0s for eventDecision, based on eventDraw
            for(eCtr=0;eCtr<27;eCtr++){
            	eventDraw[0]=1-probEvents[eCtr];
            	eventDraw[1]=probEvents[eCtr];
            	//gsl_ran_multinomial(r, 2, 1, eventDraw, outArray);
            	eventDecision[eCtr]=draw_multinom(&r, 2, eventDraw);
            //	printf("event decision %d: %d\n", eCtr, eventDecision[eCtr]);
            	//iCrawl=0;
            	//iFlag=0;
            	//while(iFlag==0){
            	//	if(outArray[iCrawl]==1){
            	//		eventDecision[eCtr]=iCrawl;
            	//		iFlag=1;
				//	}
				//	iCrawl++;
				//}
			}
			printf("event decision array prepared \n");
			
			//3. Randomise eventOrder array
			gsl_ran_shuffle(r, eventOrder, 27, sizeof(int));
			iCount=0;
			eFlag=0; //0 to keep drawing events; 1 if individual is dead 
			
			//for(eCtr=0;eCtr<27;eCtr++){
			//	printf("eventOrder %d: %d\n", eCtr, eventOrder[eCtr]);
			//}
			
			//for(eCtr=0;eCtr<27;eCtr++){
			//	printf("eventDecision %d: %d\n", eCtr, eventDecision[eCtr]);
			//}
			
            while(iCount<27 && eFlag==0){
            	iFlag=0; //flag for event search
            	iCrawl=0; //index of event to initiate
                printf("Looking for event order %d\n", iCount);
            	while(iFlag==0){
            		if(eventOrder[iCrawl]==iCount){
            			printf("to perform event %d: y or n?\n", iCrawl);
            			//if(iCrawl>27){
            			//	exit(1);
						//}
						if(eventDecision[iCrawl]==1){
							printf("PERFORM!\n");
            				 switch(iCrawl){
			                   case 0: //natural death
		                        	printf("case 0\n");
									if(target->placeInfected==1||target->placeInfected==2||target->placeInfected==3){
        		                	    deadHCVPris++;
                		        	}else if(target->placeInfected==4){
                        	    		deadHCVCom++;
                        			}
                        			removeIndiv(&target, &nHead, &nTail, pLocArray);
                        			dead++; deadFlag=1;
                        			E0++; eFlag=1;
                        			printf("natural death\n");
                        			break;
			                    case 1://stopping from injecting
			                    	printf("case 1\n");
            			            changeRisk(&target, pLocArray, 0);
                        			E1++;
                        			printf("stop injecting\n");
                        			break;
                    			case 2://inject less daily, opioid, no sharing
                        			printf("case 2\n");
									changeRisk(&target, pLocArray, 1);
                        			E2++;
                        			printf("inject less daily, no opioid, no sharing\n");
                        			break;
                    			case 3://inject less daily, opioid, sharing less than daily
                        			printf("case 3\n");
									changeRisk(&target, pLocArray, 2);
                        			E3++;
                        			printf("inject less daily, opioid, sharing less daily\n");
                        			break;
                    			case 4://inject less daily, opioid, sharing daily or more
                        			printf("case 4\n");
									E4++;
                        			printf("null event\n");
                        			break;
                    			case 5://inject less daily, non-opioid, no sharing
                        			printf("case\n");
									changeRisk(&target, pLocArray, 4);
                        			E5++;
                        			printf("inject less daily, no opioid, no sharing\n");
                        			break;
                    			case 6://inject less daily, non-opioid, sharing less than daily
                        			printf("case\n");
									changeRisk(&target, pLocArray, 5);
                        			E6++;
                        			printf("inject less daily, no opioid, sharing less daily\n");
                        			break;
                    			case 7://inject less daily, non-opioid, sharing daily or more
                       				printf("case\n");
									E7++;
                       				printf("null event\n");
                        			break;
                    			case 8://inject daily or more, opioid, no sharing
                        			printf("case\n");
									changeRisk(&target, pLocArray, 7);
                        			E8++;
                        			printf("inject daily more, opioid, no sharing\n");
                        			break;
                    			case 9://inject daily or more, opioid, sharing less than daily
                        			printf("case\n");
									changeRisk(&target, pLocArray, 8);
                        			E9++;
                        			printf("inject daily more, opioid, sharing less daily\n");
                        			break;
                    			case 10://inject daily or more, opioid, sharing daily or more
                        			printf("case\n");
									E10++;
                        			printf("null event\n");
                        			break;
                    			case 11: //injecting daily or more, non-opioid, no sharing
                        			printf("case\n");
									changeRisk(&target, pLocArray, 10);
                        			E11++;
                        			printf("inject daily more, no opioid, no sharing\n");
                        			break;
                    			case 12: //injecting daily or more, non-opioid, sharing less than daily
                        			printf("case\n");
									changeRisk(&target, pLocArray, 11);
                        			E12++;
                        			printf("inject daily more, no opioid, sharing less daily\n");
                        			break;
                    			case 13: //inject daily or more, non-opioid, sharing daily or more
                        			printf("case\n");
									E13++;
                        			printf("null event\n");
                        			break;
                    			case 14://move to community
                        			printf("case\n");
									moveLocation(&target, pLocArray, 0); //Community
                            		if(target->metavir<5){ //infected
                                		releasedHCV++;
                                		if(target->placeInfected==1||target->placeInfected==2||target->placeInfected==3){
                                    		releasedHCVPris++;
                                		}else if(target->placeInfected==4){
                                    		releasedHCVCom++;
                                		}
                            		}
                        			removeIndiv(&target, &nHead, &nTail, pLocArray); //remove from list and release to community
                        			released++;
                        			deadFlag=1; //totalIndivAgents--; eFlag=1; //to move to the next indiv
                        			E14++; eFlag=1;
                        			printf("release to community\n");
                        			break;
                    			case 15://move to prison 1
                        			printf("case\n");
									if(target->location!=1){
                            			moveLocation(&target, pLocArray, 1); //Prison 1
                        			}
                        			E15++;
                        			printf("move prison 1\n");
                        			break;
			                    case 16://move to prison 2
            			            printf("case\n");
									if(target->location!=2){
                        			    moveLocation(&target, pLocArray, 2); //Prison 2
                        			}
                        			E16++;
                        			printf("move prison 2\n");
                        			break;
                    			case 17://move to prison 3
                        			printf("case\n");
									if(target->location!=3){
                            			moveLocation(&target, pLocArray, 3); //Prison 3
                        			}
                        			E17++;
                        			printf("move prison 3\n");
                        			break;
			                    case 18://SHARING EVENT. If called, only 0.0057 (0.0032-0.0105) will lead to an infection (Boelen, L., et al., Per-event probability paper
            			            printf("case\n");
									probInf[1]=gsl_ran_exponential(r, 0.0057);
                        			//probInf[1]=gsl_ran_exponential(r, 0.0032);
                        			probInf[0]=1-probInf[1];
                        			eInf=draw_multinom(&r, 2, probInf);
                        			if(eInf==1){
                            			success=infect(&r, NEVENTS, &target, &nHead, &nTail, pLocArray, currDay);
                            			if(success==1){
                                			newCases++;
                            			}
                        			}
                        			E18++;
                        			printf("sharing event\n");
                        			break;
                    			case 19: //NATURAL CLEARANCE HCV
                        			printf("case\n");
									clearHCV(&target, pLocArray);
                        			totalClrN++;
                        			E19++;
                        			printf("natural clearance\n");
                        			break;
                    			case 20: //HCV PROGRESS
                        			printf("case\n");
									progressHCV(&target);
                        			E20++;
                        			printf("HCV progress\n");
                        			break;
                    			case 21: //DEATH BY HCV
                        			printf("case\n");
									if(target->placeInfected==1||target->placeInfected==2||target->placeInfected==3){
                            			deadHCVPris++;
                        			}else if(target->placeInfected==4){
                            			deadHCVCom++;
                        			}
                        			removeIndiv(&target, &nHead, &nTail, pLocArray);
                        			dead++;
                        			deadHCV++; deadFlag=1;
                        			E21++; eFlag=1;
                        			printf("death by HCV\n");
                        			break;
                    			case 22: //START OST THERAPY
                        			printf("case\n");
									E22++;
                        			printf("OST event\n");
                        			//totalIndivOST++;
                        			break;
                    			case 23: //STOP OST
                        			printf("case\n");
									E23++;
                        			totalIndivOST--;
                        			printf("Stop OST event\n");
                        			break;
                    			case 24: //START DAA THERAPY
                        			//totalIndivDAA++;
                        			printf("case\n");
									E24++;
                        			printf("DAA event\n");
                        			break;
                    			case 25: //STOP DAA THERAPY
                        			//totalIndivDAA--;
                        			printf("case\n");
									E25++;
                        			printf("Stop DAA event\n");
                        			break;
                    			case 26: //CLEAR VIA DAA
                        			//clearDAA(&target, pLocArray);
                        			printf("case\n");
									E26++;
                        			printf("Clear DAA event\n");
                        			break;	
								 }
						}
						iCount++;
						iFlag=1;
					}
					iCrawl++;
				}
				//exit(1);
			}
            //exit(1);
            //printf("prob0: %f\n", probEvents[0]);
            //printf("TARGET> ID: %d, ATSI: %d, Location: %d, metavir: %d, injFreq: %d, sharing: %d, shaFreq: %d, injOpd: %d\n", target->ID, target->atsi, target->location, target->metavir, target->injFreq, target->sharing, target->shaFreq, target->injOpd);
            //while loop for initiating events; keep drawing events until eFlag becomes 1
/*            eFlag=0; //0 to keep drawing events, 1 to stop
            cEvents=0; //Set counter of events to 0
            while(eFlag==0){
                //multinomial function to select event; get index of event to do (index with 1 on it not 0)
                eIndex=draw_multinom(&r, NEVENTS, probEvents);
                    
                switch(eIndex){
                   case 0: //natural death
                        //printf("Index 0 ");
                        //printf("DIE (naturally)!!!!\n");
                        if(target->placeInfected==1||target->placeInfected==2||target->placeInfected==3){
                            deadHCVPris++;
                        }else if(target->placeInfected==4){
                            deadHCVCom++;
                        }
                        removeIndiv(&target, &nHead, &nTail, pLocArray);
                        dead++; deadFlag=1;
                        //deadHCV deadHCVCom deadHCVPris
                        E0++;
                        break;
                    case 1://stopping from injecting
                        //printf("Index 1 ");
                        //printf("Stop injecting!!!!\n");
                        changeRisk(&target, pLocArray, 0);
                        E1++;
                        break;
                    case 2://inject less daily, opioid, no sharing
                        //printf("Index 2 ");
                        //printf("Inject less daily, opioid, no sharing!!!!\n");
                        changeRisk(&target, pLocArray, 1);
                        E2++;
                        break;
                    case 3://inject less daily, opioid, sharing less than daily
                        //printf("Index 3 ");
                        //printf("Inject less daily, opioid, sharing less than daily!!!!\n");
                        changeRisk(&target, pLocArray, 2);
                        E3++;
                        break;
                    case 4://inject less daily, opioid, sharing daily or more
                        //printf("Index 4 ");
                        //printf("Inject less daily, opioid, sharing daily or more!!!!\n");
                        //changeRisk(&target, pLocArray, 3);
                        E4++;
                        break;
                    case 5://inject less daily, non-opioid, no sharing
                        //printf("Index 5 ");
                        //printf("Inject less daily, non-opioid, no sharing!!!!\n");
                        changeRisk(&target, pLocArray, 4);
                        E5++;
                        break;
                    case 6://inject less daily, non-opioid, sharing less than daily
                        //printf("Index 6 ");
                        //printf("Inject less daily, non-opioid, sharing less than daily!!!!\n");
                        changeRisk(&target, pLocArray, 5);
                        E6++;
                        break;
                    case 7://inject less daily, non-opioid, sharing daily or more
                        //printf("Index 7 ");
                        //printf("Inject less daily, non-opioid, sharing daily or more!!!!\n");
                        //changeRisk(&target, pLocArray, 6);
                        E7++;
                        break;
                    case 8://inject daily or more, opioid, no sharing
                        //printf("Index 8 ");
                        //printf("Inject daily or more, opioid, no sharing!!!!\n");
                        changeRisk(&target, pLocArray, 7);
                        //E8++;
                        break;
                    case 9://inject daily or more, opioid, sharing less than daily
                        //printf("Index 9 ");
                        //printf("Inject daily or more, opioid, sharing less than daily!!!!\n");
                        changeRisk(&target, pLocArray, 8);
                        E9++;
                        break;
                    case 10://inject daily or more, opioid, sharing daily or more
                        //printf("Index 10 ");
                        //printf("inject daily or more, opioid, sharing daily or more!!!!\n");
                        //changeRisk(&target, pLocArray, 9);
                        E10++;
                        break;
                    case 11: //injecting daily or more, non-opioid, no sharing
                        //printf("Index 11 ");
                        //printf("injecting daily or more, non-opioid, no sharing!!!!\n");
                        changeRisk(&target, pLocArray, 10);
                        E11++;
                        break;
                    case 12: //injecting daily or more, non-opioid, sharing less than daily
                        //printf("Index 12 ");
                        //printf("injecting daily or more, non-opioid, sharing less than daily!!!!\n");
                        changeRisk(&target, pLocArray, 11);
                        E12++;
                        break;
                    case 13: //inject daily or more, non-opioid, sharing daily or more
                        //printf("Index 13 ");
                        //printf("injecting daily or more, non-opioid, sharing less than daily!!!!\n");
                        //changeRisk(&target, pLocArray, 12);
                        E13++;
                        break;
                    case 14://move to community
                        //printf("Index 14 \n ");
                        //printf("MOVE LOCATION???\n");
                        moveLocation(&target, pLocArray, 0); //Community
                            if(target->metavir<5){ //infected
                                releasedHCV++;
                                if(target->placeInfected==1||target->placeInfected==2||target->placeInfected==3){
                                    releasedHCVPris++;
                                }else if(target->placeInfected==4){
                                    releasedHCVCom++;
                                }
                            }
                        removeIndiv(&target, &nHead, &nTail, pLocArray); //remove from list and release to community
                        released++;
                        deadFlag=1; //totalIndivAgents--; eFlag=1; //to move to the next indiv
                        E14++;
                        break;
                    case 15://move to prison 1
                        //printf("Index 15 \n ");
                        //printf("MOVE LOCATION???\n");
                        //if already in prison 1, don't move; stay.
                        if(target->location!=1){
                            moveLocation(&target, pLocArray, 1); //Prison 1
                        }
                        E15++;
                        break;
                    case 16://move to prison 2
                        //printf("Index 16 \n ");
                        //printf("MOVE LOCATION???\n");
                        if(target->location!=2){
                            moveLocation(&target, pLocArray, 2); //Prison 2
                        }
                        E16++;
                        break;
                    case 17://move to prison 3
                        //printf("Index 17 \n ");
                        //printf("MOVE LOCATION???\n");
                        if(target->location!=3){
                            moveLocation(&target, pLocArray, 3); //Prison 3
                        }
                        E17++;
                        break;
                    case 18://SHARING EVENT. If called, only 0.0057 (0.0032-0.0105) will lead to an infection (Boelen, L., et al., Per-event probability paper
                        //printf("Index 18 \n");
                        //printf("INFECT SOMEONE!!!\n");
                        //if sharing, what is the probability of infecting? use exponential mean per-event probability
                        //probInf[1]=gsl_ran_flat(r, 0.0032, 0.0105);
                        probInf[1]=gsl_ran_exponential(r, 0.0057);
                        //probInf[1]=gsl_ran_exponential(r, 0.0032);
                        probInf[0]=1-probInf[1];
                        eInf=draw_multinom(&r, 2, probInf);
                        if(eInf==1){
                            success=infect(&r, NEVENTS, &target, &nHead, &nTail, pLocArray, currDay);
                            if(success==1){
                                newCases++;
                            }
                        }
                        E18++;
                        //traverse(&nHead, &nTail);
                        break;
                    case 19: //NATURAL CLEARANCE HCV
                        //printf("Index 19 \n ");
                        //printf("CLEAR NATURALLY!!!\n");
                        clearHCV(&target, pLocArray);
                        totalClrN++;
                        E19++;
                        break;
                    case 20: //HCV PROGRESS
                        //printf("Index 20 \n ");
                        //printf("HCV Progress\n");
                        progressHCV(&target);
                        E20++;
                        break;
                    case 21: //DEATH BY HCV
                        //printf("Index 21 \n ");
                        //printf("DEATH BY HCV!!!\n");
                        if(target->placeInfected==1||target->placeInfected==2||target->placeInfected==3){
                            deadHCVPris++;
                        }else if(target->placeInfected==4){
                            deadHCVCom++;
                        }
                        removeIndiv(&target, &nHead, &nTail, pLocArray);
                        dead++;
                        deadHCV++; deadFlag=1;
                        E21++;
                        break;
                    case 22: //START OST THERAPY
                        //printf("Index 22 \n ");
                        //printf("COMMENCE OST\n");
                        //startOST(&target, currDay);
                        E22++;
                        //totalIndivOST++;
                        break;
                    case 23: //STOP OST
                        //printf("Index 23 \n ");
                        //printf("STOP OST\n");
                        //stopOST(&target);
                        //totalIndivOST--;
                        E23++;
                        break;
                    case 24: //START DAA THERAPY
                        //printf("Index 24 \n ");
                        //printf("COMMENCE DAA\n");
                        //startDAA(&target, currDay);
                        //totalIndivDAA++;
                        E24++;
                        break;
                    case 25: //STOP DAA THERAPY
                        //printf("Index 25 \n ");
                        //printf("STOP DAA\n");
                        //stopDAA(&target);
                        //totalIndivDAA--;
                        E25++;
                        break;
                    case 26: //CLEAR VIA DAA
                        //printf("Index 20 \n ");
                        //printf("CLEAR VIA DAA\n");
                        //clearDAA(&target, pLocArray);
                        E26++;
                        break;
                }
                //cEvents++;
                //if NEVENTS is == total number of events
                //if(cEvents==NEVENTS){
                    //eFlag=1;
                //} for implementing more than one event at a time
                
                eFlag=1; //implement only 1 event at a time
            }
*/            
            
            //Move to the next individual (if current individual is still alive)
            if(deadFlag==0){
                target=target->nextIndiv;
            }
            deadFlag=0;

        }
        
        //while loop for events ends here
        //print counts and individuals
        fprintf(fw, "%d\t", currDay);
        for(i=0; i<ROWPRIS; i++) //4
        {
            for(j=0; j<COLCTR; j++){ //12
                //printf("LOCATION %d GROUP %d >>> %d \n", i, j, locArray[i][j]);
                fprintf(fw, "%d\t", locArray[i][j]);
                if(i!=0){
                    totalPrisonPop=totalPrisonPop+locArray[i][j];
                }
            }
        }
        
        //count totalPrisPop
        
        
        //count prevalence in prison
        if(nTail!=NULL){
            prisonHCV=countHCVpris(&nHead, &nTail);
        }
        
        //count prevalence in prison obtained in community
        if(nTail!=NULL){
            communityHCV=countHCVCom(&nHead, &nTail);
        }
        
        //count Rs
        if(nTail!=NULL){
            R0=count0(&nHead, &nTail);
            R1=count1(&nHead, &nTail);
            R2=count2(&nHead, &nTail);
            R3=count3(&nHead, &nTail);
            R4=count4(&nHead, &nTail);
            R5=count5(&nHead, &nTail);
            R6=count6(&nHead, &nTail);
            R7=count7(&nHead, &nTail);
            R8=count8(&nHead, &nTail);
            R9=count9(&nHead, &nTail);
            R10=count10(&nHead, &nTail);
            R11=count11(&nHead, &nTail);
            R12=count12(&nHead, &nTail);
        }
        
        fprintf(fw, "%d\t", newCases); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", prisonHCV); //prevalence of those infected in prisons
        fprintf(fw, "%d\t", communityHCV); //prevalence of those infected in community
        fprintf(fw, "%d\t", released); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", releasedHCV); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", releasedHCVCom); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", releasedHCVPris); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", totalClrN); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", dead); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", deadHCV); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", deadHCVCom); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", deadHCVPris); //print HCV-infected individuals released per day
        fprintf(fw, "%d\t", totalPrisonPop);
        fprintf(fw, "%d\t", E0);
        fprintf(fw, "%d\t", E1);
        fprintf(fw, "%d\t", E2);
        fprintf(fw, "%d\t", E3);
        fprintf(fw, "%d\t", E4);
        fprintf(fw, "%d\t", E5);
        fprintf(fw, "%d\t", E6);
        fprintf(fw, "%d\t", E7);
        fprintf(fw, "%d\t", E8);
        fprintf(fw, "%d\t", E9);
        fprintf(fw, "%d\t", E10);
        fprintf(fw, "%d\t", E11);
        fprintf(fw, "%d\t", E12);
        fprintf(fw, "%d\t", E13);
        fprintf(fw, "%d\t", E14);
        fprintf(fw, "%d\t", E15);
        fprintf(fw, "%d\t", E16);
        fprintf(fw, "%d\t", E17);
        fprintf(fw, "%d\t", E18);
        fprintf(fw, "%d\t", E19);
        fprintf(fw, "%d\t", E20);
        fprintf(fw, "%d\t", E21);
        fprintf(fw, "%d\t", E22);
        fprintf(fw, "%d\t", E23);
        fprintf(fw, "%d\t", E24);
        fprintf(fw, "%d\t", E25);
        fprintf(fw, "%d\t", E26);
        fprintf(fw, "%d\t", R0);
        fprintf(fw, "%d\t", R1);
        fprintf(fw, "%d\t", R2);
        fprintf(fw, "%d\t", R3);
        fprintf(fw, "%d\t", R4);
        fprintf(fw, "%d\t", R5);
        fprintf(fw, "%d\t", R6);
        fprintf(fw, "%d\t", R7);
        fprintf(fw, "%d\t", R8);
        fprintf(fw, "%d\t", R9);
        fprintf(fw, "%d\t", R10);
        fprintf(fw, "%d\t", R11);
        fprintf(fw, "%d\n", R12);
        newCases=0;
        prisonHCV=0; communityHCV=0;
        released=0; releasedHCV=0; releasedHCVCom=0; releasedHCVPris=0; //set counter back to zero
        totalClrN=0;
        deadHCV=0; deadHCVCom=0; deadHCVPris=0;
        dead=0; totalPrisonPop=0;
        E0=0; E1=0; E2=0; E3=0; E4=0; E5=0; E6=0; E7=0; E8=0; E9=0; E10=0; E11=0; E12=0; E13=0; E14=0; E15=0; E16=0; E17=0; E18=0; E19=0; E20=0; E21=0; E22=0; E23=0; E24=0; E25=0; E26=0;
        R0=0; R1=0; R2=0; R3=0; R4=0; R5=0; R6=0; R7=0; R8=0; R9=0; R10=0; R11=0; R12=0;
        //printf("Total : %d \n", locArray[1][0]+locArray[2][0]+locArray[3][0]);
        
        if(nTail!=NULL){
            traverse(&nHead, &nTail);
        }
        
        currDay++;
        //printf("\n");
    }

	printf("DONE!\n");
    fclose(fw);
    gsl_rng_free(r);
    return 0;
}

