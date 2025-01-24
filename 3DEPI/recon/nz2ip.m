function IP=nz2ip(nz,Rz)
IP = -nz/2:Rz:nz/2-1;
negIP=IP(IP<0);
posIP=IP(IP>=0);
IP_co=zeros(size(IP)); %center out
negIP=fliplr(negIP);

IP_co(1)=negIP(1);
if length(negIP)>length(posIP)
    %length(IP) is odd in this case
    for i=2:length(IP)
        if mod(i,2)==0
            IP_co(i)=negIP(1+i/2);
        else
            IP_co(i)=posIP(floor(i/2));
        end
    end
elseif length(negIP)==length(posIP)
    %length(IP) is even in this case
    for i=2:length(IP)
        if mod(i,2)==0
            IP_co(i)=posIP(i/2);
        else
            IP_co(i)=negIP(1+floor(i/2));
        end
    end
end
IP=IP_co;
IP=IP+nz/2+1;
end