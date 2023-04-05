% find the mask this function gets
function [probability,X2,Y2,boundaries ,BW ]= findWhiskers(frame_size,framei,dlc,GuassFilterSize,numberofwhiskers,XYbiasDLC)
    if numberofwhiskers ==1
        dlc.L3=dlc.L2;% make the 3rd whisker the second one in the absence of 3 whiskers
        dlc.x3=dlc.x2;
        dlc.y3 = dlc.y2;
    end
    dlc_this_frame = dlc(dlc.frames==framei,:);
    if numberofwhiskers ==1
        [probability,X2,Y2] = read_dlc(dlc_this_frame,3,XYbiasDLC,frame_size);
    else
        [probability,X2,Y2] = read_dlc(dlc_this_frame,numberofwhiskers,XYbiasDLC,frame_size);
    end
    if sum(probability<0.01)<2
        thisim = zeros(frame_size);
        ind = sub2ind(frame_size,floor((Y2)),floor((X2)));
        thisim(ind(~isnan(ind)))=1;
        B = imgaussfilt(thisim,GuassFilterSize);
        BW = imbinarize(B);
        boundaries = bwboundaries(BW);
        theselenghts=cellfun(@length,boundaries);
        [~,Maxmatind] = max(theselenghts);
        boundaries  = boundaries{Maxmatind};
    else
        BW= zeros(frame_size);
        boundaries=[];
    end
end
function [probability,X2,Y2] = read_dlc(dlc_this_frame,numberofwhiskers,XYbiasDLC,frame_size)
    probability = arrayfun(@(i)dlc_this_frame.(['L' num2str(i)]),1:numberofwhiskers);
    X2 = arrayfun(@(i)dlc_this_frame.(['L' num2str(i)]),1:numberofwhiskers);
    Y2 = arrayfun(@(i)dlc_this_frame.(['L' num2str(i)]),1:numberofwhiskers);
    Y2(Y2<1)=1;
    X2(X2<1)=1;
    X2=X2+XYbiasDLC(1);
    Y2=Y2+XYbiasDLC(2);
    Y2(probability<0.01) = NaN;
    X2(probability<0.01) = NaN;
    Y2(Y2>(frame_size(1)))=frame_size(1);
    X2(X2>(frame_size(2)))=frame_size(2);
end