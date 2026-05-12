function plot_handle = ciplot(lower,upper,x,colour,alpha);
if nargin<5
    alpha=0.5;
end
if nargin<4
    colour='b';
end
if nargin<3
    x=1:length(lower);
end
% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end
plot_handle = fill([x fliplr(x)],[upper fliplr(lower)],colour,'FaceAlpha',alpha,"EdgeColor","none")
