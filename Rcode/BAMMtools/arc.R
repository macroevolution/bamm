##################################
#	Internal function called by plot.dtrates(...)
#	Arguments:
#		x,y = coordinates of center of curvature of arc, e.g. (0,0)
#		theta1 = initial theta of arc (radians)
#		theta2 = ending theta of arc (radians)
#		rad = radius of arc
arc = function(x,y,theta1,theta2,rad,border,...)
{
	step = (theta2-theta1)/100;
	xv = x+rad*cos(seq(theta1,theta2,step));
	yv = y+rad*sin(seq(theta1,theta2,step));
	if(step) lines(xv,yv,lend=2,col=border,...);
}

arc = Vectorize(arc);