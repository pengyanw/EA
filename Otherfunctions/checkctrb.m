function flag = checkctrb(A,B)
if rank(ctrb(A,B))==size(A,2)
    flag=1;
    disp('Controllable!')
else
    flag=0;
     disp('Not Controllable!')
end