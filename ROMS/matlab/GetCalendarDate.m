function h = GetCalendarDate(seconds,days,calendartype)
%% h = GetCalendarDate(seconds,days [,calendartype] )
%
% seconds, days are the DART times 
%
% EXAMPLE:
%
% mydate = GetCalendarDate(82761,148520);

%% DART software - Copyright 2004 - 2011 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% <next few lines under version control, do not edit>
% $URL$
% $Id$
% $Revision$
% $Date$

if (nargin < 2) 
   error('Must supply at least two arguments, seconds and days')
elseif (nargin ==2)
   calendartype = 'Gregorian';
end

switch lower(calendartype)
   case 'gregorian'
      mytime = datenum(1601,1,1) + days + seconds/86400;
      h = datestr(mytime);
      fprintf('DART time (%d s, %d d) is %s %s\n', ...
                     seconds,days,h,calendartype)
   case 'noleap'
      error('noleap not supported yet')
   case 'thirty_day_months'
      error('thirty_day_months not supported yet')
   case 'julian'
      error('julian not supported yet')
   case 'no_calendar'
      error('no_calendar not supported yet')
   case 'gregorian_mars'
      error('gregorian_mars not supported yet')
   otherwise
end
