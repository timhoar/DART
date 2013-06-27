/*
 * generate the merge and query fortran namelists
 * in a loop for all of sept.   merges original ncep
 * obs files with new gps files.
 *
 * nsc 23jul08
 */

#if 0
! example of what text needs to be generated:
&query_obs_seq_nml
   filename_seq = '/ptmp/nancy/ncep/obs_seq20060906'  /

&merge_obs_seq_nml
   num_input_files = 5,
   filename_seq = '../obs_seq.gpsro_2006091806', 
                  '../obs_seq.gpsro_2006091812', 
                  '../obs_seq.gpsro_2006091818', 
                  '../obs_seq.gpsro_2006091824', 
                  '/ptmp/nancy/ncep/obs_seq20060918', 
   filename_out = '/ptmp/nancy/newgps/obs_seq20060918',
   first_obs_days           = 148183,
   first_obs_seconds        = 10801,
   last_obs_days            = 148184,
   last_obs_seconds         = 10800,
   /
#endif

#include "stdio.h"
#include "stdlib.h"

void addme(char *thing);

main(int argc, char **argv)
{
    int i;
    int gregbase = 148408;
    int startd = 1;
    int endd = 31;
    int month = 5;
    int year = 2006;
    char buf[128];

    /* set start, end dates, plus gregorian base */
    if (argc > 4) {
        month = atoi(argv[1]);
        gregbase = atoi(argv[2]);
        startd = atoi(argv[3]);
        endd = atoi(argv[4]);
    }

    printf("base=%d, year=%d, month=%d, start=%d, end=%d\n", 
            gregbase, year, month, startd, endd);

    for (i=startd; i<=endd; i++) {
        /* get a fresh copy of the namelist each time, 
         * then append a merge and query section and run.
         */ 

        system("cp -f input.template input.nml");

        addme("&merge_obs_seq_nml");
        addme("   num_input_files = 5,");

        sprintf(buf, "   filename_seq = \'../obs_seq.gpsro_%4d%02d%02d06\', ", 
                year, month, i);
        addme(buf);
        sprintf(buf, "                  \'../obs_seq.gpsro_%4d%02d%02d12\', ",
                year, month, i);
        addme(buf);
        sprintf(buf, "                  \'../obs_seq.gpsro_%4d%02d%02d18\', ", 
                year, month, i);
        addme(buf);
        sprintf(buf, "                  \'../obs_seq.gpsro_%4d%02d%02d24\', ",
                year, month, i);
        addme(buf);
        sprintf(buf, "                  \'/ptmp/nancy/ncep/obs_seq%4d%02d%02d\', ", 
                year, month, i);
        addme(buf);
        sprintf(buf, "   filename_out = \'/ptmp/nancy/newgps/obs_seq%4d%02d%02d\', ", 
                year, month, i);
        addme(buf);
        sprintf(buf, "   first_obs_days           = %d,", gregbase+i-1);
        addme(buf);
        addme("   first_obs_seconds        = 10801,");
        sprintf(buf, "   last_obs_days            = %d, ", gregbase+i);
        addme(buf);
        addme("   last_obs_seconds         = 10800 ");
    
        addme("   /");
        addme("");
        addme("&query_obs_seq_nml");
        sprintf(buf, "   filename_seq = \'/ptmp/nancy/newgps/obs_seq%4d%02d%02d\'", 
                year, month, i);
        addme(buf);
        addme("   /");

        sprintf(buf, "echo day %d", i);
        system(buf);
        system("./merge_obs_seq");
        system("./query_obs_seq");

    }    

    exit(0);
}

void addme(char *thing)
{
    char cmd[256];

    printf("echo \"%s\"\n", thing);
    sprintf(cmd, "echo \"%s\" >> input.nml", thing);
    system(cmd);
}

