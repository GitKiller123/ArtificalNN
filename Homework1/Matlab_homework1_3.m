clear all;

trialstot = 100000;
mtot = 100;

betha = 2;
g = @(b)(1/(1+exp(-2*betha*b)));
n = 200;
p = 40; 
m = zeros(mtot,1);
for mcount = 1:mtot
for k = 1:length(p)
    choice = p(k);
    S_org = zeros(choice,n);
    
    w = zeros(n,n);
    for x = 1:choice
        temp = randi([0 1],1,n)*2-1;
        S_org(x, :) = temp;
        w = w + temp'*temp/n;
    end
    w = w - diag(diag(w));
    pnum = 1;
    nnum = 1;
    M = zeros(trialstot,1);
    S = S_org(pnum,:);
    for trials = 1:trialstot
%         nnum = randi(n, 1);
        b = w(nnum, :)*S(pnum, :)';
        vrand = rand;
        
        if vrand <= g(b)
            nout = 1;
        else
            nout = -1;
        end
        if nout ~= S(nnum)
            S(nnum) = nout;
        end
%         if nout == S(pnum, nnum)
%             continue
%         else
%             errcount = errcount + 1;
%         end
        M(trials) = S_org(pnum,:)*S'/n;
        if nnum < n
            nnum = nnum + 1;
        else
            nnum = 1;
        end
    end
end
m(mcount) = sum(M)/trialstot;
end
avgm = sum(m)/mtot;