function OD = Smooth(OD,thresh) 
%     SmoothOD = zeros(size(OD));
%     for i = 2:size(OD,1)-1
%         for j = 2:size(OD,2)-1             
%             difs = [OD(i,j)-OD(i+1,j),OD(i,j)-OD(i+1,j-1),OD(i,j)-OD(i+1,j+1)...
%                     ,OD(i,j)-OD(i,j-1),OD(i,j)-OD(i,j+1)...
%                     ,OD(i,j)-OD(i-1,j),OD(i,j)-OD(i-1,j-1),OD(i,j)-OD(i-1,j+1)];            
%             if sum(abs(difs) > thresh)                
%                 SmoothOD(i,j) = OD(i,j) + (sum(difs) - 8*OD(i,j))/8.0;
%             else
%                 SmoothOD(i,j) = OD(i,j);
%             end
%         end
%     end   
kernel = ones(3,3)/8;
kernel(2,2) = 0;
movingAverage = conv2(OD, kernel);
movingAverage = movingAverage(2:end-1,2:end-1);
d = movingAverage - OD;
problem_index = abs(d) > thresh;
OD(problem_index) = movingAverage(problem_index);
end