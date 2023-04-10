
function [TLeft,TRight]=Smooth_whiskerDLCTLeft(TLeft,TRight,allnumberofwhiskers)
smoothwindowsize = 30;
for jj=1:numel(allnumberofwhiskers)
    numberofwhiskers = allnumberofwhiskers(jj);
switch numberofwhiskers
    case 1
        A= TLeft.x1;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        TLeft.x1 = C;
        
    case 2
        A= TLeft.x2;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
                TLeft.x2 = C;

    case 3
        A= TLeft.x3;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        TLeft.x3 = C;
    case 4
        A= TLeft.x4;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        TLeft.x4 = C;
        
end
switch numberofwhiskers
    case 1
        A= TRight.x1;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        
        TRight.x1 = C;
    case 2
        A= TRight.x2;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        TRight.x2 = C;
        
    case 3
        A= TRight.x3;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        TRight.x3 = C;
        
    case 4
        
        A= TRight.x4;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        TRight.x4 = C;
        
end
%%
switch numberofwhiskers
    case 1
        
        A= TLeft.y1;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        TLeft.y1 = C;
    case 2
        
        A= TLeft.y2;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        
        TLeft.y2 = C;
    case 3
        
        A= TLeft.y3;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        
        TLeft.y3 = C;
    case 4
        
        A= TLeft.y4;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        
        TLeft.y4 = C;
end


switch numberofwhiskers
    case 1
        
        A= TRight.y1;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        
        TRight.y1 = C;
    case 2
        
        A= TRight.y2;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        
        TRight.y2 = C;
    case 3
        
        A= TRight.y3;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        
        TRight.y3 = C;
    case 4
        
        A= TRight.y4;
        [B, TF] = rmoutliers(A);
        C = NaN(size(A));
        C (~TF)=B; C = nanfastsmooth(C,smoothwindowsize,1,1);
        C(isnan(C))=A(isnan(C));
        
        TRight.y4 = C;
end
end