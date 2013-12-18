##################################
#	Internal function called by plot.dtrates(...)
#	Arguments:
#		x,y = coordinates of center of curvature of arc, e.g. (0,0)
#		theta1 = initial theta of arc (radians)
#		theta2 = ending theta of arc (radians)
#		rad = radius of arc
arc = function(x,y,theta1,theta2,rad,border,...)
{
	steps = (theta2-theta1)/100; noTips = steps > 0;
	steps = steps[noTips];
	theta1 = theta1[noTips];
	theta2 = theta2[noTips];
	rad = rad[noTips];
	border = border[noTips];
	for (i in 1:length(steps))
	{
		xv = x+rad[i]*cos(seq(theta1[i],theta2[i],steps[i]));
		yv = y+rad[i]*sin(seq(theta1[i],theta2[i],steps[i]));
		lines(xv,yv,lend=2,col=border[i],...);		
	}
}
