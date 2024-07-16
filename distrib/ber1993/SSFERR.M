function zern=ssferr(f)

global capt xi

zern=exp(f/(capt-xi-f))-capt/(capt-xi-f) ;
