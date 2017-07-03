//
//  individual.h
//  HCV Model
//
//  Created by Neil Bretana on 1/10/2014.
//  Copyright (c) 2014 Neil Bretana. All rights reserved.
//

#ifndef HCV_Model_individual_h
#define HCV_Model_individual_h

struct individual
{
    int ID;
    int gender;
    int age; //18-100
    int atsi; //0-no; 1-yes
    int injecting; //0-no; 1-yes;
    int sharing; //0-no; 1-yes;
    int injFreq; //0 - not applicable; 1 less than daily; 2 daily or more;
    int shaFreq; //0 - not applicable; 1 less than daily; 2 daily or more;
    int location; //0, 1, 2, 3
    int group; //0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 This is for the population groups
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
    int placeInfected; //0 (N/A), 1 (min), 2 (med), 3 (max), 4 (community)
    int risk; //0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 This is for the injecting risk group
    //Risk groups: 0: non-injecting; 1: injecting less than daily; opioid; not sharing 2: injecting less than daily; opioid; sharing less daily 3: injecting less than daily; opioid; sharing daily more 4: injecting less than daily; non-opioid; no sharing 5: injecting less than daily; non-opioid; sharing less daily 6: injecting less than daily; non-opioid; sharing daily or more 7: ...
    int timeOfImprisonment; //Record time of imprisonment
    
    int infectionNumber; //Tells the number of times an individual has been infected
    int metavir; //0-4; 5-clearer; 6-non-infected//never-infected
    int cirrhosis; //0 and 1 link to metavir 4
    int viralLoad;
    int timeOfInfection; //Record timepoint when it got infected
    
    int injOpd; //Injecting Opioids 0 1
    int her; //0 1
    int herDose;
    int herDaysMissed;
    int met; //0 1
    int metDose;
    int metDaysMissed;
    int bup;
    int bupDose;
    int bupDaysMissed;
    int OST;
    int timeStartOST; //record time of OST start
    int bleach; //0 1
    int syringeEx; //0 1
    int DAA; //0 1
    int timeStartDAA; //record time of DAA start
    
    struct individual *prevIndiv;
    struct individual *nextIndiv; //nextIndiv
};

typedef struct individual sIndiv;

//Methods

#endif
