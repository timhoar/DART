#!/bin/bash

ele=`expr $1 - 1`

v1=( -77.64109263 -78.26319591 -78.88378501 -79.50144611 -80.11477639 -80.72238923 -81.32291942 -81.91502833 -82.49740887 -83.06879045 -83.62794375 -84.17368535 -84.70488220 -85.22045591 -85.71938687 -86.20071810 -86.66355892 -87.10708838 -87.53055842 -87.93329679 -88.31470966 -88.67428402 -89.01158970 -89.32628117 -89.61809898 -89.88687091 -90.13251283 -90.35502918 -90.55451321 -90.73114684 -90.88520020 -91.01703090 -91.12708291 -91.21588521 -91.28405006 -91.33227097 -91.36132045 -91.37204739 -91.36537416 -91.34229351 -91.30386515 -91.25121207 -91.18551671 -91.10801678 -91.02000102 -90.92280463 -90.81780460 -90.70641488 -90.59008135 -90.47027675 -90.34849540 -90.22624793 -90.10505584 -89.98644611 -89.87194565 -89.76307590 -89.66134723 -89.56825358 -89.48526695 -89.41383210 -89.35536122 -89.31122878 -89.28276644 -89.27125815 -89.27793533 -89.30397228 -89.35048178 -89.41851082 -89.50903664 -89.62296292 -89.76111628 -89.92424296 -90.11300587 -90.32798182 -90.56965912 -90.83843539 -91.13461579 -91.45841145 -91.80993829 -92.18921613 -92.59616817 -93.03062072 -93.49230336 -93.98084932 -94.49579631 -95.03658756 -95.60257330 -96.19301244 -96.80707464 -97.44384270 -98.10231520 -98.78140944 -99.47996470 -100.19674570 -100.93044643 -101.67969410 -102.44305340 -103.21903097 -104.00608008 -104.80260543 -105.60696823 -106.41749136 -107.23246468 -108.05015049 -108.86878902 -109.68660411 -110.50180888 -111.31261141 -112.11722057 -112.91385169 -113.70073240 -114.47610825 -115.23824845 -115.98545137 -116.71605010 -117.42841780 -118.12097298 -118.79218455 -119.44057686 -120.06473436 -120.66330625 -121.23501080 -121.77863947 -122.29306082 -122.77722414 -123.23016281 -123.65099739 -124.03893845 -124.39328903 -124.71344686 -124.99890626 -125.24925969 -125.46419898 -125.64351626 -125.78710455 -125.89495795 -125.96717164 -126.00394136 -126.00556272 -125.97243010 -125.90503523 -125.80396544 -125.66990166 -125.50361601 -125.30596919 -125.07790752 -124.82045968 -124.53473326 -124.22191102 -123.88324682 -123.52006148 -123.13373828 -122.72571835 -122.29749582 -121.85061286 -121.38665450 -120.90724338 -120.41403434 -119.90870896 -119.39296997 -118.86853561 -118.33713404 -117.80049757 -117.26035705 -116.71843616 -116.17644583 -115.63607861 -115.09900327 -114.56685934 -114.04125183 -113.52374612 -113.01586297 -112.51907361 -112.03479520 -111.56438628 -111.10914260 -110.67029307 -110.24899603 -109.84633568 -109.46331889 -109.10087218 -108.75983905 -108.44097760 -108.14495841 -107.87236277 -107.62368122 -107.39931241 -107.19956224 -107.02464340 -106.87467518 -106.74968360 -106.64960194 -106.57427149 -106.52344270 -106.49677661 -106.49384657 -106.51414034 -106.55706238 -106.62193650 -106.70800878 -106.81445072 -106.94036267 -107.08477752 -107.24666457 -107.42493364 -107.61843942 -107.82598591 -108.04633113 -108.27819192 -108.52024885 -108.77115137 -109.02952288 -109.29396603 -109.56306801 -109.83540590 -110.10955202 -110.38407932 -110.65756675 -110.92860455 -111.19579955 -111.45778035 -111.71320244 -111.96075320 -112.19915677 -112.42717878 -112.64363094 -112.84737540 -113.03732899 -113.21246717 -113.37182788 -113.51451500 -113.63970171 -113.74663347 -113.83463084 -113.90309194 -113.95149468 -113.97939860 -113.98644658 -113.97236602 -113.93696994 -113.88015755 -113.80191468 -113.70231375 -113.58151354 -113.43975855 -113.27737807 -113.09478500 -112.89247425 -112.67102094 -112.43107824 -112.17337493 -111.89871273 -111.60796330 -111.30206497 -110.98201933 -110.64888745 -110.30378593 -109.94788281 -109.58239316 -109.20857459 -108.82772254 -108.44116549 -108.05025995 -107.65638542 -107.26093923 -106.86533129 -106.47097878 -106.07930087 -105.69171332 -105.30962317 -104.93442341 -104.56748765 -104.21016498 -103.86377472 -103.52960146 -103.20889002 -102.90284073 -102.61260472 -102.33927941 -102.08390426 -101.84745658 -101.63084766 -101.43491909 -101.26043931 -101.10810044 -100.97851531 -100.87221488 -100.78964582 -100.73116846 -100.69705499 -100.68748805 -100.70255952 -100.74226973 -100.80652692 -100.89514707 -101.00785402 -101.14427996 -101.30396617 -101.48636415 -101.69083707 -101.91666147 -102.16302932 -102.42905037 -102.71375479 -103.01609611 -103.33495442 -103.66913981 -104.01739614 -104.37840495 -104.75078971 -105.13312018 -105.52391701 -105.92165656 -106.32477582 -106.73167751 -107.14073536 -107.55029939 -107.95870139 -108.36426046 -108.76528854 -109.16009609 -109.54699766 -109.92431761 -110.29039566 -110.64359250 -110.98229535 -111.30492336 -111.60993300 -111.89582326 -112.16114082 -112.40448498 -112.62451242 -112.81994190 -112.98955856 -113.13221823 -113.24685130 -113.33246651 -113.38815441 -113.41309054 -113.40653836 -113.36785192 -113.29647813 -113.19195885 -113.05393256 -112.88213579 -112.67640414 -112.43667309 -112.16297833 -111.85545591 -111.51434192 -111.13997193 -110.73278005 -110.29329769 -109.82215197 -109.32006384 -108.78784585 -108.22639965 -107.63671317 -107.01985752 -106.37698362 -105.70931853 -105.01816161 -104.30488029 -103.57090580 -102.81772853 -102.04689327 -101.25999429 -100.45867015 -99.64459853 -98.81949074 -97.98508631 -97.14314731 -96.29545277 -95.44379287 -94.58996328 -93.73575930 -92.88297017 -92.03337330 -91.18872859 -90.35077276 -89.52121387 -88.70172581 -87.89394302 -87.09945524 -86.31980256 -85.55647050 -84.81088533 -84.08440963 -83.37833804 -82.69389323 -82.03222212 -81.39439239 -80.78138925 -80.19411244 -79.63337362 -79.09989396 -78.59430211 -78.11713245 -77.66882368 -77.24971772 -76.86005892 -76.49999366 -76.16957023 -75.86873907 -75.59735329 -75.35516960 -75.14184951 -74.95696084 -74.79997959 -74.67029209 -74.56719744 -74.48991027 -74.43756374 -74.40921287 -74.40383809 -74.42034899 -74.45758844 -74.51433674 -74.58931618 -74.68119559 -74.78859521 -74.91009163 -75.04422293 -75.18949386 -75.34438117 -75.50733906 -75.67680459 -75.85120325 -76.02895446 -76.20847714 -76.38819520 -76.56654310 -76.74197124 -76.91295139 -77.07798191 -77.23559298 -77.38435165 -77.52286672 -77.64979350 -77.76383842 -77.86376336 -77.94838984 -78.01660299 -78.06735525 -78.09966985 -78.11264402 -78.10545192 -78.07734729 -78.02766586 -77.95582738 -77.86133744 -77.74378885 -77.60286288 -77.43833000 -77.25005044 -77.03797437 -76.80214173 -76.54268180 -76.25981240 -75.95383883 -75.62515242 -75.27422882 -74.90162604 -74.50798207 -74.09401233 -73.66050681 -73.20832693 -72.73840215 -72.25172636 -71.74935402 -71.23239611 -70.70201586 -70.15942432 -69.60587570 -69.04266269 -68.47111147 -67.89257677 -67.30843674 -66.72008772 -66.12893904 -65.53640772 -64.94391310 -64.35287160 -63.76469133 -63.18076690 -62.60247411 -62.03116488 -61.46816214 -60.91475489 -60.37219338 -59.84168444 -59.32438691 -58.82140732 -58.33379576 -57.86254185 -57.40857111 -56.97274138 -56.55583959 -56.15857879 -55.78159543 -55.42544685 -55.09060922 -54.77747558 -54.48635437 -54.21746808 -53.97095242 -53.74685561 -53.54513811 -53.36567264 -53.20824449 -53.07255222 -52.95820858 -52.86474183 -52.79159733 -52.73813943 -52.70365367 -52.68734928 -52.68836194 -52.70575682 -52.73853187 -52.78562142 -52.84589992 -52.91818597 -53.00124658 -53.09380153 -53.19452804 -53.30206548 -53.41502033 -53.53197119 -53.65147397 -53.77206713 -53.89227697 -54.01062311 -54.12562382 -54.23580152 -54.33968821 -54.43583090 -54.52279699 -54.59917962 -54.66360291 -54.71472713 -54.75125374 -54.77193032 -54.77555534 -54.76098276 -54.72712644 -54.67296440 -54.59754281 -54.49997978 -54.37946891 -54.23528259 -54.06677503 -53.87338499 -53.65463826 -53.41014987 -53.13962588 -52.84286500 -52.51975982 -52.17029772 -51.79456145 -51.39272943 -50.96507569 -50.51196944 -50.03387441 -49.53134776 -49.00503871 -48.45568689 -47.88412031 -47.29125307 -46.67808273 -46.04568742 -45.39522265 -44.72791790 -44.04507282 -43.34805337 -42.63828752 -41.91726093 -41.18651227 -40.44762841 -39.70223947 -38.95201361 -38.19865184 -37.44388254 -36.68945596 -35.93713868 -35.18870788 -34.44594567 -33.71063338 -32.98454575 -32.26944528 -31.56707647 -30.87916019 -30.20738812 -29.55341722 -28.91886440 -28.30530121 -27.71424875 -27.14717273 -26.60547867 -26.09050734 -25.60353035 -25.14574608 -24.71827569 -24.32215951 -23.95835363 -23.62772681 -23.33105763 -23.06903195 -22.84224073 -22.65117806 -22.49623959 -22.37772126 -22.29581838 -22.25062497 -22.24213352 -22.27023505 -22.33471950 -22.43527644 -22.57149616 -22.74287102 -22.94879718 -23.18857663 -23.46141951 -23.76644674 -24.10269295 -24.46910973 -24.86456902 -25.28786695 -25.73772777 -26.21280807 -26.71170124 -27.23294209 -27.77501171 -28.33634240 -28.91532292 -29.51030370 -30.11960227 -30.74150877 -31.37429155 -32.01620278 -32.66548416 -33.32037266 -33.97910618 -34.63992928 -35.30109887 -35.96088977 -36.61760026 -37.26955755 -37.91512309 -38.55269777 -39.18072696 -39.79770547 -40.40218216 -40.99276455 -41.56812306 -42.12699512 -42.66818900 -43.19058738 -43.69315072 -44.17492023 -44.63502072 -45.07266298 -45.48714602 -45.87785887 -46.24428217 -46.58598934 -46.90264750 -47.19401805 -47.45995691 -47.70041442 -47.91543494 -48.10515612 -48.26980783 -48.40971078 -48.52527486 -48.61699711 -48.68545944 -48.73132608 -48.75534068 -48.75832317 -48.74116643 -48.70483258 -48.65034918 -48.57880507 -48.49134615 -48.38917085 -48.27352551 -48.14569955 -48.00702053 -47.85884906 -47.70257367 -47.53960554 -47.37137312 -47.19931687 -47.02488376 -46.84952190 -46.67467517 -46.50177781 -46.33224913 -46.16748823 -46.00886888 -45.85773439 -45.71539274 -45.58311171 -45.46211431 -45.35357423 -45.25861159 -45.17828888 -45.11360705 -45.06550190 -45.03484070 -45.02241902 -45.02895791 -45.05510127 -45.10141354 -45.16837774 -45.25639369 -45.36577666 -45.49675624 -45.64947558 -45.82399094 -46.02027151 -46.23819962 -46.47757124 -46.73809675 -47.01940213 -47.32103036 -47.64244315 -47.98302300 -48.34207555 -48.71883215 -49.11245276 -49.52202913 -49.94658815 -50.38509553 -50.83645964 -51.29953562 -51.77312962 -52.25600326 -52.74687829 -53.24444135 -53.74734885 -54.25423203 -54.76370208 -55.27435532 -55.78477846 -56.29355395 -56.79926525 -57.30050220 -57.79586631 -58.28397606 -58.76347213 -59.23302255 -59.69132776 -60.13712559 -60.56919607 -60.98636615 -61.38751418 -61.77157433 -62.13754068 -62.48447124 -62.81149161 -63.11779855 -63.40266316 -63.66543393 -63.90553941 -64.12249067 -64.31588345 -64.48539997 -64.63081056 -64.75197478 -64.84884246 -64.92145418 -64.96994165 -64.99452759 -64.99552542 -64.97333852 -64.92845923 -64.86146755 -64.77302946 -64.66389496 -64.53489587 -64.38694322 -64.22102445 -64.03820028 -63.83960135 -63.62642459 -63.39992929 -63.16143309 -62.91230756 -62.65397375 -62.38789747 -62.11558437 -61.83857498 -61.55843949 -61.27677250 -60.99518763 -60.71531205 -60.43878099 -60.16723212 -59.90230000 -59.64561048 -59.39877512 -59.16338561 -58.94100834 -58.73317892 -58.54139690 -58.36712050 -58.21176157 -58.07668061 -57.96318201 -57.87250941 -57.80584135 -57.76428701 -57.74888229 -57.76058605 -57.80027665 -57.86874871 -57.96671021 -58.09477980 -58.25348447 -58.44325750 -58.66443673 -58.91726309 -59.20187959 -59.51833046 -59.86656080 -60.24641639 -60.65764398 -61.09989182 -61.57271058 -62.07555452 -62.60778307 -63.16866271 -63.75736908 -64.37298953 -65.01452589 -65.68089749 -66.37094455 -67.08343182 -67.81705237 -68.57043178 -69.34213243 -70.13065809 -70.93445866 -71.75193507 -72.58144446 -73.42130533 -74.26980301 -75.12519509 -75.98571704 -76.84958785 )


echo ${v1[$ele]}
