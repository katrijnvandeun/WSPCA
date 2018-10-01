function YSTAR=ystar(XHAT,RES,wWSQ)
%YSTAR calculates the fixed part of the majorizing WLS function (Kiers,
%1997 and 2002)
%YSTAR=XHAT+wWSQ.*RES
%K. Van Deun, Dept. Psychology, KU Leuven
%version 1: March 2013
YSTAR=XHAT+wWSQ.*RES;