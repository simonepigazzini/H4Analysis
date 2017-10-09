double modelFunction(double* x, double* par)
{
    double f1=0, f2=0;
    
    f1 = par[0]*(1-1/(par[1]*log(x[0]/par[2])+1));

    if(f1 <= f2)
	return 0;
    else
	return f1;
}
	
