function h=m_quiver(long,lat,u,v,varargin)
global MAP_PROJECTION MAP_VAR_LIST
[X,Y]=m_ll2xy(long,lat,'clip','point');
[XN ,YN ]=m_ll2xy([long(:) long(:)]',[lat(:) lat(:)+.001]','clip','off');
[XE ,YE ]=m_ll2xy([long(:) long(:)+(.001)./cos(lat(:)*pi/180)]',[lat(:) lat(:)]','clip','off');
mU=u.*reshape(diff(XE),size(lat))*1000 + v.*reshape(diff(XN),size(lat))*1000;
mV=u.*reshape(diff(YE),size(lat))*1000 + v.*reshape(diff(YN),size(lat))*1000;
ii=isfinite(X(:)) & isfinite(u(:));
h=quiver(X(ii),Y(ii),mU(ii),mV(ii),varargin{:});
set(h,'tag','m_quiver');
if nargout==0
 clear h
end
