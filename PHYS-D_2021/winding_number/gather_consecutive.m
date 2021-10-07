function out = gather_consecutive(v)
out = NaN*ones(length(v),length(v));
global_tick = 0;
ticker = 1;
while ticker <= length(v)
    while ticker <= length(v) && v(ticker)==0
        ticker = ticker+1;
    end
    slide = 0;
    while ticker + slide <= length(v) && v(ticker+slide)~=0
        slide = slide + 1;
    end
    if ticker+slide-1<=length(v) && slide>=1
        out(1:slide,global_tick+1) = v(ticker:ticker+slide-1);
        global_tick = global_tick + 1;
        ticker = ticker + slide;
    end
end
out = out(:,1:global_tick);
end