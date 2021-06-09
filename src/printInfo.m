function printInfo(iterCount,val,grad,relerr,step)

%print reporting information

if(nargin < 1)
    %print the header
    fprintf(1,'\n\nStarting solver: %s\n\n',mfilename);
    fprintf(1,'iterate  |  objective value  | gradient value |  relative improvement | step size \n');
    fprintf(1,'---------------------------------------------------------------------\n');
else
    fprintf(1,'%4i     |  %0.6e    |  %0.6e    |  %0.6e    |  %0.3e \n',iterCount,val,grad,relerr,step);
end