	  t  O   k820309    4          12.1        2�@Q                                                                                                           
       ../../../reg_factor/reg_factor_mod.f90 REG_FACTOR_MOD              COMP_REG_FACTOR                      @                              
       R8                      @                              
       GET_UNIT OPEN_FILE REGISTER_MODULE ERROR_HANDLER E_ERR E_MSG NMLFILEUNIT FIND_NAMELIST_IN_FILE CHECK_NAMELIST_READ DO_NML_FILE DO_NML_TERM                      @                              
       TIME_TYPE WRITE_TIME GET_TIME                �  !@                                '                    #SECONDS    #DAYS                 � D                                                              � D                                                                                                                                                      %         @                                                            %         @                                	                          #OPEN_FILE%LEN 
   #OPEN_FILE%PRESENT    #OPEN_FILE%TRIM    #OPEN_FILE%MIN    #FNAME    #FORM    #ACTION                  @                            
     LEN               @                                 PRESENT               @                                 TRIM               @                                 MIN           
   @                                                 1           
  @                                                 1           
  @                                                 1 #         @                                                   #REGISTER_MODULE%TRIM    #SRC    #REV    #RDATE                  @                                 TRIM           
   @                                                 1           
   @                                                 1           
   @                                                 1 #         @                                                	   #ERROR_HANDLER%PRESENT    #ERROR_HANDLER%TRIM    #LEVEL    #ROUTINE    #TEXT    #SRC    #REV    #RDATE    #AUT    #TEXT2     #TEXT3 !                 @                                 PRESENT               @                                 TRIM           
   @                                                   
   @                                                 1           
   @                                                 1           
  @                                                 1           
  @                                                 1           
  @                                                 1           
  @                                                 1           
  @                                                  1           
  @                             !                    1                                              "                                                      2                                             #                                                       0               !                            $            #         @                                 %                  #FIND_NAMELIST_IN_FILE%PRESENT &   #FIND_NAMELIST_IN_FILE%TRIM '   #FIND_NAMELIST_IN_FILE%ADJUSTL (   #NAMELIST_FILE_NAME )   #NML_NAME *   #IUNIT +   #WRITE_TO_LOGFILE_IN ,                 @                            &     PRESENT               @                            '     TRIM               @                            (     ADJUSTL           
   @                             )                    1           
   @                             *                    1             @                              +                      
  @                              ,           #         @                                 -                  #CHECK_NAMELIST_READ%LEN .   #CHECK_NAMELIST_READ%PRESENT /   #CHECK_NAMELIST_READ%TRIM 0   #CHECK_NAMELIST_READ%INDEX 1   #IUNIT 2   #IOSTAT_IN 3   #NML_NAME 4   #WRITE_TO_LOGFILE_IN 5                 @                            .     LEN               @                            /     PRESENT               @                            0     TRIM               @                            1     INDEX           
   @                              2                     
   @                              3                     
   @                             4                    1           
  @                              5           %         @                                6                            %         @                                7                            #         @                                  8                  #WRITE_TIME%PRESENT 9   #WRITE_TIME%TRIM :   #FILE_UNIT ;   #TIME <   #FORM =   #IOS_OUT >                 @                            9     PRESENT               @                            :     TRIM           
   @                              ;                     
   @                              <                   #TIME_TYPE              
  @                             =                    1             @                              >            #         @                                 ?                  #GET_TIME%PRESENT @   #GET_TIME%HUGE A   #TIME B   #SECONDS C   #DAYS D                 @                            @     PRESENT               @                            A     HUGE           
   @                              B                   #TIME_TYPE                @                              C                        @                              D            %         @                                E                   
       #COMP_REG_FACTOR%SUM F   #NUM_GROUPS G   #REGRESS H   #OBS_TIME I   #OBS_INDEX J   #STATE_INDEX K   #OBS_STATE_IND L   #OBS_STATE_MAX M                 @                            F     SUM           
   @                              G                    
   @                             H                    
    p          5 � p        r G       5 � p        r G                               
  @@                              I                   #TIME_TYPE              
   @                              J                     
   @                              K                     
  @                              L                     
  @                              M              �   >      fn#fn $   �       b   uapp(REG_FACTOR_MOD    �   C   J  TYPES_MOD    A  �   J  UTILITIES_MOD !     ^   J  TIME_MANAGER_MOD +   j  g       TIME_TYPE+TIME_MANAGER_MOD ;   �  H   %   TIME_TYPE%SECONDS+TIME_MANAGER_MOD=SECONDS 5     H   %   TIME_TYPE%DAYS+TIME_MANAGER_MOD=DAYS    a  p       R8+TYPES_MOD '   �  P       GET_UNIT+UTILITIES_MOD (   !  �       OPEN_FILE+UTILITIES_MOD 0   �  <      OPEN_FILE%LEN+UTILITIES_MOD=LEN 8     @      OPEN_FILE%PRESENT+UTILITIES_MOD=PRESENT 2   _  =      OPEN_FILE%TRIM+UTILITIES_MOD=TRIM 0   �  <      OPEN_FILE%MIN+UTILITIES_MOD=MIN .   �  L   e   OPEN_FILE%FNAME+UTILITIES_MOD -   $  L   e   OPEN_FILE%FORM+UTILITIES_MOD /   p  L   e   OPEN_FILE%ACTION+UTILITIES_MOD .   �         REGISTER_MODULE+UTILITIES_MOD 8   ;  =      REGISTER_MODULE%TRIM+UTILITIES_MOD=TRIM 2   x  L   e   REGISTER_MODULE%SRC+UTILITIES_MOD 2   �  L   e   REGISTER_MODULE%REV+UTILITIES_MOD 4     L   e   REGISTER_MODULE%RDATE+UTILITIES_MOD ,   \  �       ERROR_HANDLER+UTILITIES_MOD <   5	  @      ERROR_HANDLER%PRESENT+UTILITIES_MOD=PRESENT 6   u	  =      ERROR_HANDLER%TRIM+UTILITIES_MOD=TRIM 2   �	  @   e   ERROR_HANDLER%LEVEL+UTILITIES_MOD 4   �	  L   e   ERROR_HANDLER%ROUTINE+UTILITIES_MOD 1   >
  L   e   ERROR_HANDLER%TEXT+UTILITIES_MOD 0   �
  L   e   ERROR_HANDLER%SRC+UTILITIES_MOD 0   �
  L   e   ERROR_HANDLER%REV+UTILITIES_MOD 2   "  L   e   ERROR_HANDLER%RDATE+UTILITIES_MOD 0   n  L   e   ERROR_HANDLER%AUT+UTILITIES_MOD 2   �  L   e   ERROR_HANDLER%TEXT2+UTILITIES_MOD 2     L   e   ERROR_HANDLER%TEXT3+UTILITIES_MOD $   R  q       E_ERR+UTILITIES_MOD $   �  q       E_MSG+UTILITIES_MOD *   4  @       NMLFILEUNIT+UTILITIES_MOD 4   t  �       FIND_NAMELIST_IN_FILE+UTILITIES_MOD D   l  @      FIND_NAMELIST_IN_FILE%PRESENT+UTILITIES_MOD=PRESENT >   �  =      FIND_NAMELIST_IN_FILE%TRIM+UTILITIES_MOD=TRIM D   �  @      FIND_NAMELIST_IN_FILE%ADJUSTL+UTILITIES_MOD=ADJUSTL G   )  L   e   FIND_NAMELIST_IN_FILE%NAMELIST_FILE_NAME+UTILITIES_MOD =   u  L   e   FIND_NAMELIST_IN_FILE%NML_NAME+UTILITIES_MOD :   �  @   e   FIND_NAMELIST_IN_FILE%IUNIT+UTILITIES_MOD H     @   e   FIND_NAMELIST_IN_FILE%WRITE_TO_LOGFILE_IN+UTILITIES_MOD 2   A        CHECK_NAMELIST_READ+UTILITIES_MOD :   E  <      CHECK_NAMELIST_READ%LEN+UTILITIES_MOD=LEN B   �  @      CHECK_NAMELIST_READ%PRESENT+UTILITIES_MOD=PRESENT <   �  =      CHECK_NAMELIST_READ%TRIM+UTILITIES_MOD=TRIM >   �  >      CHECK_NAMELIST_READ%INDEX+UTILITIES_MOD=INDEX 8   <  @   e   CHECK_NAMELIST_READ%IUNIT+UTILITIES_MOD <   |  @   e   CHECK_NAMELIST_READ%IOSTAT_IN+UTILITIES_MOD ;   �  L   e   CHECK_NAMELIST_READ%NML_NAME+UTILITIES_MOD F     @   e   CHECK_NAMELIST_READ%WRITE_TO_LOGFILE_IN+UTILITIES_MOD *   H  P       DO_NML_FILE+UTILITIES_MOD *   �  P       DO_NML_TERM+UTILITIES_MOD ,   �  �       WRITE_TIME+TIME_MANAGER_MOD <   �  @      WRITE_TIME%PRESENT+TIME_MANAGER_MOD=PRESENT 6   �  =      WRITE_TIME%TRIM+TIME_MANAGER_MOD=TRIM 6   
  @   e   WRITE_TIME%FILE_UNIT+TIME_MANAGER_MOD 1   J  W   e   WRITE_TIME%TIME+TIME_MANAGER_MOD 1   �  L   e   WRITE_TIME%FORM+TIME_MANAGER_MOD 4   �  @   e   WRITE_TIME%IOS_OUT+TIME_MANAGER_MOD *   -  �       GET_TIME+TIME_MANAGER_MOD :   �  @      GET_TIME%PRESENT+TIME_MANAGER_MOD=PRESENT 4   �  =      GET_TIME%HUGE+TIME_MANAGER_MOD=HUGE /   <  W   e   GET_TIME%TIME+TIME_MANAGER_MOD 2   �  @   e   GET_TIME%SECONDS+TIME_MANAGER_MOD /   �  @   e   GET_TIME%DAYS+TIME_MANAGER_MOD       �       COMP_REG_FACTOR $   �  <      COMP_REG_FACTOR%SUM +   )  @   a   COMP_REG_FACTOR%NUM_GROUPS (   i  �   a   COMP_REG_FACTOR%REGRESS )     W   a   COMP_REG_FACTOR%OBS_TIME *   t  @   a   COMP_REG_FACTOR%OBS_INDEX ,   �  @   a   COMP_REG_FACTOR%STATE_INDEX .   �  @   a   COMP_REG_FACTOR%OBS_STATE_IND .   4  @   a   COMP_REG_FACTOR%OBS_STATE_MAX 