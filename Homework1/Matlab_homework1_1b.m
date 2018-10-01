clear all;

trialstot = 100000;

n = 100;
p = [12,20,40,60,80,100]; 
errlist = zeros(1,length(p));
for k = 1:length(p)
    choice = p(k);
    errcount = 0;
    currp = zeros(choice,n);
    for trials = 1:trialstot
        w = zeros(n,n);
        for x = 1:choice
            temp = randi([0 1],1,n)*2-1;
            currp(x, :) = temp;
            w = w + temp'*temp/n;
        end

%         w = w - diag(diag(w));
        pnum = randi(choice, 1);
        nnum = randi(n, 1);
        
        nout = w(nnum, :)*currp(pnum, :)';

        if nout < 0
            nout = -1;
        else
            nout = 1;
        end
        if nout == currp(pnum, nnum)
            continue
        else
            errcount = errcount + 1;
        end
    end
    errperc = errcount / trialstot;
    errlist(k) = errperc;
end