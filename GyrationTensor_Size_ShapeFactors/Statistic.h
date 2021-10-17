void moment(VecDoub_I &data, Doub &ave, Doub &sdev)
{
/*Given an array of data[0..n-1], this routine returns its mean ave, average deviation adev,
standard deviation sdev, variance var, skewness skew, and kurtosis curt.*/
	Doub adev;
	Doub var;
	Doub skew;
	Doub curt;
	Int j,n=data.size();
	Doub ep=0.0,s,p;
	if (n <= 1) throw("n must be at least 2 in moment");
	s=0.0; //First pass to get the mean.
	for (j=0;j<n;j++) s += data[j];
	ave=s/n;
	adev=var=skew=curt=0.0; //Second pass to get the first (absolute), second,
	//third, and fourth moments of thedeviation from the mean.
	for (j=0;j<n;j++) {
		adev += abs(s=data[j]-ave);
		ep += s;
		var += (p=s*s);
		skew += (p *= s);
		curt += (p *= s);
	}
	adev /= n;
	var=(var-ep*ep/n)/(n-1); //Corrected two-pass formula.
	sdev=sqrt(var); //Put the pieces together according to the conif
	/*if(var != 0.0) 
	{ //ventional definitions.
		skew /= (n*var*sdev);
		curt=curt/(n*var*var)-3.0;
	} else throw("No skew/kurtosis when variance = 0 (in moment)");*/
}