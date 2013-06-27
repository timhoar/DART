%  read AIRS Level-2 retrievals of T and moisture (mixing ratio)
clf;
clear;
%
iyear=2006;
imonth=8;
beg_day = 16;
end_day = 31;
%
ayear = int2str(iyear);
amonth= int2str(imonth);
%
for iday=beg_day : end_day;

totnum = 0;
aday  = int2str(iday);

if (imonth < 10)
bmonth = strcat('0', amonth);
end
inname = strcat('ls -1 data/AIRS.', ayear, '.', bmonth, '.', aday, '* >list');
catname = strcat('cat out', ayear, '.', bmonth, '.', aday, '* >satemp_obs.',ayear,bmonth,aday);

[s,w] = system('rm -f out*');
[s,w] = system(inname);

fid = fopen('list', 'r');
allstr = fscanf(fid, '%s', [1 inf]);   % read in all filenames in a single row
fclose(fid);
%
str_length = 74;                       % may vary at different paths.
nfiles = size(allstr,2)/str_length;
fff=reshape(allstr,str_length, nfiles);
fname_all=fff';
%
for nn=1:size(fname_all,1)
fname=fname_all(nn, :)
outname=fname_all(nn, 11:24);
tt  = read_L12_swath_file(fname, 2);
airs= read_L12_swath_file(fname, 3);
%
year   = tt.start_year(1:1);
month  = tt.start_month(1:1);
day    = tt.start_day(1:1);
hour   = tt.start_hour(1:1);
minute = tt.start_minute(1:1) + 3.0 ;    % 6 minute window for each granule
second = tt.start_sec(1:1);
time  = double(hour) + double(minute)/60.0 + double(second)/3600.0;
%
lat  = airs.Latitude;       % (deg)
lon  = airs.Longitude;
land = airs.landFrac;       % (0 -> 1.0)
qa   = airs.RetQAFlag;
cld  = airs.CldFrcStd;      % (0 -> 1.0)
%
tair = airs.TAirStd;        % (K)
qair = airs.H2OMMRStd;      % (g/kg)
%
%  the following 12 of 28 levels data are extracted.
pres = [1000. 925. 850. 700. 600. 500. 400. 300. 250. 200. 150];

out = zeros(37800, 6);
kk = 0;
%
nprofile = 0;
for i = 1:size(lat, 1)
for j = 1:size(lat, 2)
   if(qa(i,j) == 0)
   if(land(i,j) < 0.001)
%
   totcld = cld(1,1,1,i,j) + cld(1,2,1,i,j) + cld(1,3,1,i,j) + ...
            cld(1,1,2,i,j) + cld(1,2,2,i,j) + cld(1,3,2,i,j) + ...
            cld(1,1,3,i,j) + cld(1,2,3,i,j) + cld(1,3,3,i,j) ;
%
   if(totcld > 0.0001)
   if(totcld <= 0.01)
%    output temperature 
      nprofile = nprofile + 1; 
     for k = 1:11;
       kk = (nprofile-1)*22 + 2*k-1;
      totnum = totnum + 1;
      out(kk, 1) = lon(i,j);
      out(kk, 2) = lat(i,j);
      out(kk, 3) = pres(k);          %  hPa
      out(kk, 4) = double(tair(k+1,i,j));
      out(kk, 5) = totnum + 1000000;
      out(kk, 6) = time;

%   for mixing ratio
       kk = (nprofile-1)*22 + 2*k;
      totnum = totnum + 1;
      out(kk, 1) = lon(i,j);
      out(kk, 2) = lat(i,j);
      out(kk, 3) = pres(k);
      out(kk, 4) = double(qair(k+1,i,j));
      out(kk, 5) = totnum + 5000000;
      out(kk, 6) = time;
     end
   end 
   end 
   end 
   end 
end
end
%

if (kk > 0)
out2 = out(1:kk, 1:6);
oname = strcat('out', outname);
save (oname, 'out2', '-ASCII');         % print out each single granule file
end

end                                     %  end of the day files loop
%
[ss,ww] = system(catname);              % print out to daily files.
[s,w] = system('rm -f out*');
end                                     %  end of the day loop
%
