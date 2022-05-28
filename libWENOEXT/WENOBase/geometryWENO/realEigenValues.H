/*---------------------------------------------------------------------------*\
       ██╗    ██╗███████╗███╗   ██╗ ██████╗     ███████╗██╗  ██╗████████╗
       ██║    ██║██╔════╝████╗  ██║██╔═══██╗    ██╔════╝╚██╗██╔╝╚══██╔══╝
       ██║ █╗ ██║█████╗  ██╔██╗ ██║██║   ██║    █████╗   ╚███╔╝    ██║   
       ██║███╗██║██╔══╝  ██║╚██╗██║██║   ██║    ██╔══╝   ██╔██╗    ██║   
       ╚███╔███╔╝███████╗██║ ╚████║╚██████╔╝    ███████╗██╔╝ ██╗   ██║   
        ╚══╝╚══╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝     ╚══════╝╚═╝  ╚═╝   ╚═╝   
-------------------------------------------------------------------------------                                                                                                                                                      
License
    This file is part of WENO Ext.

    WENO Ext is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    WENO Ext is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with  WENO Ext.  If not, see <http://www.gnu.org/licenses/>.

Description
    Calculation of eigen values of a non symmetric real matrix A
    Based on Numerical Recipes 3rd edition 

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022

\*---------------------------------------------------------------------------*/


#ifndef mathFunctionsEigenValue_H 
#define mathFunctionsEigenValue_H

#include "geometryWENO.H"
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace mathFunctionsWENO
{
/*---------------------------------------------------------------------------*\
                        Class to Calculate Eigen Values
\*---------------------------------------------------------------------------*/

class realEigenValues 
{

    using complexVector = 
        blaze::DynamicVector<blaze::complex<double>,blaze::columnVector>;

    private:
    
    // Private member variables
    
        //- Number of rows
        const int n;
        
        geometryWENO::scalarSquareMatrix a;
        
        //- Eigen values
        complexVector sigma;

    // Private member functions
     
        //- Balance the matrix 
        //  Algorithm of Obsourne is used
        void balance();
        
        //- Reduce to Hessenberg form 
        void elmhes();
        
        //- QR decomposition
        void hqr();

        //- Sort the eigen values
        void sort();
        
        //- Returns the magnitude of a times the sign of b 
        double sign(const double& a, const double& b) const;

    public:
        //- Constructor with given matrix 
        realEigenValues(const geometryWENO::scalarSquareMatrix& aa) 
        :
            n(aa.rows()), a(aa), sigma(n)
        {
            // Balance matrix
            balance();
            
            // Reduce to Hessenberg form
            elmhes();
            
            // QR decompositon
            hqr();
            
            // Sort the eigen values
            sort();
        }
        
        //- Return the eigenvalues
        inline const complexVector& eig() {return sigma;}

        //- Return the condition of the matrix 
        inline scalar cond() 
        {
            if (min(abs(sigma)) < ROOTVSMALL)
                return GREAT;
            
            return max(abs(sigma))/min(abs(sigma));
        }
};

} // End of mathFunctionsWENO
} // End of Foam

void Foam::mathFunctionsWENO::realEigenValues::balance()
{
	const double RADIX = std::numeric_limits<double>::radix;
	bool done = false;
	double sqrdx = RADIX*RADIX;
	while (!done) 
    {
		done = true;
		for (int i=0; i < n; i++) 
        {
			double r = 0.0, c = 0.0;
			for (int j=0; j < n; j++)
				if (j != i) 
                {
					c += abs(a(j,i));
					r += abs(a(i,j));
				}
			if (c != 0.0 && r != 0.0) 
            {
				double g = r/RADIX;
				double f = 1.0;
				double s = c+r;
				while (c<g) 
                {
					f *= RADIX;
					c *= sqrdx;
				}
				g = r*RADIX;
				while (c>g) 
                {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s) 
                {
					done = false;
					g = 1.0/f;
                    // Not needed as we do not calculate eigen vectors
					// scale[i] *= f;
					for (int j=0;j<n;j++) a(i,j) *= g;
					for (int j=0;j<n;j++) a(j,i) *= f;
				}
			}
		}
	}
}


void Foam::mathFunctionsWENO::realEigenValues::elmhes()
{
    // Function to swap to values
    auto swap = [](double& a, double& b)
    {
        double temp = a;
        a = b;
        b = temp;
    };
    
    
    
	for (int m = 1; m < n-1; m++) 
    {
		double x = 0.0;
		int i=m;
		for (int j = m; j<n; j++) 
        {
			if (abs(a(j,m-1)) > abs(x)) 
            {
				x=a(j,m-1);
				i=j;
			}
		}
        // Not needed because we do not calculate eigen vectors
		//perm[m] = i;
		if (i != m) 
        {
			for (int j = m-1; j<n; j++) swap(a(i,j),a(m,j));
			for (int j = 0; j<n; j++) swap(a(j,i),a(j,m));
		}
		if (x != 0.0) 
        {
			for (i=m+1; i<n; i++) {
				double y=a(i,m-1);
				if (y != 0.0)
                {
					y /= x;
					a(i,m-1) = y;
					for (int j=m;j<n;j++) a(i,j) -= y*a(m,j);
					for (int j=0;j<n;j++) a(j,m) += y*a(j,i);
				}
			}
		}
	}
}


void Foam::mathFunctionsWENO::realEigenValues::hqr()
{
	int nn,m,l,k,j,its,i,mmin;
	double z,y,x,w,v,u,t,s,r,q,p,anorm=0.0;

	const double EPS = std::numeric_limits<double>::epsilon();
    
	for (i=0; i<n; i++)
		for (j=std::max(i-1,0); j<n; j++)
			anorm += std::abs(a(i,j));
	nn = n-1;
	t = 0.0;
	while (nn >= 0) 
    {
		its=0;
		do 
        {
			for (l=nn; l>0; l--) 
            {
				s=abs(a(l-1,l-1))+abs(a(l,l));
				if (s == 0.0) s=anorm;
                
				if (abs(a(l,l-1)) <= EPS*s) 
                {
					a(l,l-1) = 0.0;
					break;
				}
			}
			
            x=a(nn,nn);
			
            if (l == nn) 
            {
				sigma[nn--]=x+t;
			} 
            else 
            {
				y = a(nn-1,nn-1);
				w = a(nn,nn-1)*a(nn-1,nn);
				if (l == nn-1) 
                {
					p = 0.5*(y-x);
					q = p*p+w;
					z = std::sqrt(abs(q));
					x += t;
					if (q >= 0.0) 
                    {
						z= p + sign(z,p);
						sigma[nn-1] = sigma[nn] = x+z;
						if (z != 0.0) 
                            sigma[nn] = x-w/z;
					} 
                    else 
                    {
						sigma[nn] = blaze::complex<double>(x+p,-z);
						sigma[nn-1]=conj(sigma[nn]);
					}
					nn -= 2;
				} 
                else 
                {
					if (its == 30) 
                        FatalError 
                            << "Too many iterations in hqr" 
                            << exit(FatalError); 
					if (its == 10 || its == 20) 
                    {
						t += x;
						for (i=0; i<nn+1; i++) 
                            a(i,i) -= x;
						s = std::abs(a(nn,nn-1)) + std::abs(a(nn-1,nn-2));
						y = x= 0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=nn-2; m>=l; m--) 
                    {
						z = a(m,m);
						r = x-z;
						s = y-z;
						p = (r*s-w)/a(m+1,m)+a(m,m+1);
						q = a(m+1,m+1)-z-r-s;
						r = a(m+2,m+1);
						s = std::abs(p)+std::abs(q)+std::abs(r);
						p /= s;
						q /= s;
						r /= s;
                        
						if (m == l) 
                            break;
                            
						u = std::abs(a(m,m-1))*(std::abs(q)+std::abs(r));
						v = std::abs(p)*(std::abs(a(m-1,m-1))+std::abs(z)+std::abs(a(m+1,m+1)));
						
                        if (u <= EPS*v) 
                            break;
					}
					for (i=m; i<nn-1; i++) 
                    {
						a(i+2,i) = 0.0;
						if (i != m) 
                            a(i+2,i-1) = 0.0;
					}
					for (k=m;k<nn;k++) 
                    {
						if (k != m) 
                        {
							p = a(k,k-1);
							q = a(k+1,k-1);
							r = 0.0;
							if (k+1 != nn) 
                                r = a(k+2,k-1);
							if ((x=abs(p)+abs(q)+abs(r)) != 0.0) 
                            {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=sign(std::sqrt(p*p+q*q+r*r),p)) != 0.0) 
                        {
							if (k == m) 
                            {
								if (l != m)
                                    a(k,k-1) = -a(k,k-1);
							} 
                            else
                            {
								a(k,k-1) = -s*x;
                            }
                            
							p += s;
							x = p/s;
							y = q/s;
							z = r/s;
							q /= p;
							r /= p;
							for (j=k;j<nn+1;j++) 
                            {
								p=a(k,j)+q*a(k+1,j);
								if (k+1 != nn) 
                                {
									p += r*a(k+2,j);
									a(k+2,j) -= p*z;
								}
								a(k+1,j) -= p*y;
								a(k,j) -= p*x;
							}
                            
							mmin = nn < k+3 ? nn : k+3;
                            
							for (i=l; i<mmin+1; i++) 
                            {
								p=x*a(i,k)+y*a(i,k+1);
								if (k+1 != nn) 
                                {
									p += z*a(i,k+2);
									a(i,k+2) -= p*r;
								}
								a(i,k+1) -= p*q;
								a(i,k) -= p;
							}
						}
					}
				}
			}
		} while (l+1 < nn);
	}
}


//void Unsymmeig::hqr2()
//{
	//Int nn,m,l,k,j,its,i,mmin,na;
	//Doub z,y,x,w,v,u,t,s,r,q,p,anorm=0.0,ra,sa,vr,vi;

	//const Doub EPS=numeric_limits<Doub>::epsilon();
	//for (i=0;i<n;i++)
		//for (j=MAX(i-1,0);j<n;j++)
			//anorm += abs(a[i][j]);
	//nn=n-1;
	//t=0.0;
	//while (nn >= 0) {
		//its=0;
		//do {
			//for (l=nn;l>0;l--) {
				//s=abs(a[l-1][l-1])+abs(a[l][l]);
				//if (s == 0.0) s=anorm;
				//if (abs(a[l][l-1]) <= EPS*s) {
					//a[l][l-1] = 0.0;
					//break;
				//}
			//}
			//x=a[nn][nn];
			//if (l == nn) {
				//sigma[nn]=a[nn][nn]=x+t;
				//nn--;
			//} else {
				//y=a[nn-1][nn-1];
				//w=a[nn][nn-1]*a[nn-1][nn];
				//if (l == nn-1) {
					//p=0.5*(y-x);
					//q=p*p+w;
					//z=sqrt(abs(q));
					//x += t;
					//a[nn][nn]=x;
					//a[nn-1][nn-1]=y+t;
					//if (q >= 0.0) {
						//z=p+SIGN(z,p);
						//sigma[nn-1]=sigma[nn]=x+z;
						//if (z != 0.0) sigma[nn]=x-w/z;
						//x=a[nn][nn-1];
						//s=abs(x)+abs(z);
						//p=x/s;
						//q=z/s;
						//r=sqrt(p*p+q*q);
						//p /= r;
						//q /= r;
						//for (j=nn-1;j<n;j++) {
							//z=a[nn-1][j];
							//a[nn-1][j]=q*z+p*a[nn][j];
							//a[nn][j]=q*a[nn][j]-p*z;
						//}
						//for (i=0;i<=nn;i++) {
							//z=a[i][nn-1];
							//a[i][nn-1]=q*z+p*a[i][nn];
							//a[i][nn]=q*a[i][nn]-p*z;
						//}
						//for (i=0;i<n;i++) {
							//z=zz[i][nn-1];
							//zz[i][nn-1]=q*z+p*zz[i][nn];
							//zz[i][nn]=q*zz[i][nn]-p*z;
						//}
					//} else {
						//sigma[nn]=Complex(x+p,-z);
						//sigma[nn-1]=conj(sigma[nn]);
					//}
					//nn -= 2;
				//} else {
					//if (its == 30) throw("Too many iterations in hqr");
					//if (its == 10 || its == 20) {
						//t += x;
						//for (i=0;i<nn+1;i++) a[i][i] -= x;
						//s=abs(a[nn][nn-1])+abs(a[nn-1][nn-2]);
						//y=x=0.75*s;
						//w = -0.4375*s*s;
					//}
					//++its;
					//for (m=nn-2;m>=l;m--) {
						//z=a[m][m];
						//r=x-z;
						//s=y-z;
						//p=(r*s-w)/a[m+1][m]+a[m][m+1];
						//q=a[m+1][m+1]-z-r-s;
						//r=a[m+2][m+1];
						//s=abs(p)+abs(q)+abs(r);
						//p /= s;
						//q /= s;
						//r /= s;
						//if (m == l) break;
						//u=abs(a[m][m-1])*(abs(q)+abs(r));
						//v=abs(p)*(abs(a[m-1][m-1])+abs(z)+abs(a[m+1][m+1]));
						//if (u <= EPS*v) break;
					//}
					//for (i=m;i<nn-1;i++) {
						//a[i+2][i]=0.0;
						//if (i != m) a[i+2][i-1]=0.0;
					//}
					//for (k=m;k<nn;k++) {
						//if (k != m) {
							//p=a[k][k-1];
							//q=a[k+1][k-1];
							//r=0.0;
							//if (k+1 != nn) r=a[k+2][k-1];
							//if ((x=abs(p)+abs(q)+abs(r)) != 0.0) {
								//p /= x;
								//q /= x;
								//r /= x;
							//}
						//}
						//if ((s=SIGN(sqrt(p*p+q*q+r*r),p)) != 0.0) {
							//if (k == m) {
								//if (l != m)
								//a[k][k-1] = -a[k][k-1];
							//} else
								//a[k][k-1] = -s*x;
							//p += s;
							//x=p/s;
							//y=q/s;
							//z=r/s;
							//q /= p;
							//r /= p;
							//for (j=k;j<n;j++) {
								//p=a[k][j]+q*a[k+1][j];
								//if (k+1 != nn) {
									//p += r*a[k+2][j];
									//a[k+2][j] -= p*z;
								//}
								//a[k+1][j] -= p*y;
								//a[k][j] -= p*x;
							//}
							//mmin = nn < k+3 ? nn : k+3;
							//for (i=0;i<mmin+1;i++) {
								//p=x*a[i][k]+y*a[i][k+1];
								//if (k+1 != nn) {
									//p += z*a[i][k+2];
									//a[i][k+2] -= p*r;
								//}
								//a[i][k+1] -= p*q;
								//a[i][k] -= p;
							//}
							//for (i=0; i<n; i++) {
								//p=x*zz[i][k]+y*zz[i][k+1];
								//if (k+1 != nn) {
									//p += z*zz[i][k+2];
									//zz[i][k+2] -= p*r;
								//}
								//zz[i][k+1] -= p*q;
								//zz[i][k] -= p;
							//}
						//}
					//}
				//}
			//}
		//} while (l+1 < nn);
	//}
	//if (anorm != 0.0) {
		//for (nn=n-1;nn>=0;nn--) {
			//p=real(sigma[nn]);
			//q=imag(sigma[nn]);
			//na=nn-1;
			//if (q == 0.0) {
				//m=nn;
				//a[nn][nn]=1.0;
				//for (i=nn-1;i>=0;i--) {
					//w=a[i][i]-p;
					//r=0.0;
					//for (j=m;j<=nn;j++)
						//r += a[i][j]*a[j][nn];
					//if (imag(sigma[i]) < 0.0) {
						//z=w;
						//s=r;
					//} else {
						//m=i;
						
						//if (imag(sigma[i]) == 0.0) {
							//t=w;
							//if (t == 0.0)
								//t=EPS*anorm;
							//a[i][nn]=-r/t;
						//} else {
							//x=a[i][i+1];
							//y=a[i+1][i];
							//q=SQR(real(sigma[i])-p)+SQR(imag(sigma[i]));
							//t=(x*s-z*r)/q;
							//a[i][nn]=t;
							//if (abs(x) > abs(z))
								//a[i+1][nn]=(-r-w*t)/x;
							//else
								//a[i+1][nn]=(-s-y*t)/z;
						//}
						//t=abs(a[i][nn]);
						//if (EPS*t*t > 1)
							//for (j=i;j<=nn;j++)
								//a[j][nn] /= t;
					//}
				//}
			//} else if (q < 0.0) {
				//m=na;
				//if (abs(a[nn][na]) > abs(a[na][nn])) {
					//a[na][na]=q/a[nn][na];
					//a[na][nn]=-(a[nn][nn]-p)/a[nn][na];
				//} else {
					//Complex temp=Complex(0.0,-a[na][nn])/Complex(a[na][na]-p,q);
					//a[na][na]=real(temp);
					//a[na][nn]=imag(temp);
				//}
				//a[nn][na]=0.0;
				//a[nn][nn]=1.0;
				//for (i=nn-2;i>=0;i--) {
					//w=a[i][i]-p;
					//ra=sa=0.0;
					//for (j=m;j<=nn;j++) {
						//ra += a[i][j]*a[j][na];
						//sa += a[i][j]*a[j][nn];
					//}
					//if (imag(sigma[i]) < 0.0) {
						//z=w;
						//r=ra;
						//s=sa;
					//} else {
						//m=i;
						//if (imag(sigma[i]) == 0.0) {
							//Complex temp = Complex(-ra,-sa)/Complex(w,q);
							//a[i][na]=real(temp);
							//a[i][nn]=imag(temp);
						//} else {
							//x=a[i][i+1];
							//y=a[i+1][i];
							//vr=SQR(real(sigma[i])-p)+SQR(imag(sigma[i]))-q*q;
							//vi=2.0*q*(real(sigma[i])-p);
							//if (vr == 0.0 && vi == 0.0)
								//vr=EPS*anorm*(abs(w)+abs(q)+abs(x)+abs(y)+abs(z));
							//Complex temp=Complex(x*r-z*ra+q*sa,x*s-z*sa-q*ra)/
								//Complex(vr,vi);
							//a[i][na]=real(temp);
							//a[i][nn]=imag(temp);
							//if (abs(x) > abs(z)+abs(q)) {
								//a[i+1][na]=(-ra-w*a[i][na]+q*a[i][nn])/x;
								//a[i+1][nn]=(-sa-w*a[i][nn]-q*a[i][na])/x;
							//} else {
								//Complex temp=Complex(-r-y*a[i][na],-s-y*a[i][nn])/
									//Complex(z,q);
								//a[i+1][na]=real(temp);
								//a[i+1][nn]=imag(temp);
							//}
						//}
					//}
					//t=MAX(abs(a[i][na]),abs(a[i][nn]));
					//if (EPS*t*t > 1)
						//for (j=i;j<=nn;j++) {
							//a[j][na] /= t;
							//a[j][nn] /= t;
						//}
				//}
			//}
		//}
		//for (j=n-1;j>=0;j--)
			//for (i=0;i<n;i++) {
				//z=0.0;
				//for (k=0;k<=j;k++)
					//z += zz[i][k]*a[k][j];
				//zz[i][j]=z;
			//}
	//}
//}



void Foam::mathFunctionsWENO::realEigenValues::sort()
{
	int i;
	for (int j=1; j<n; j++) 
    {
		blaze::complex<double> x = sigma[j];
		for (i=j-1; i>=0; i--) 
        {
			if (real(sigma[i]) >= real(x)) 
                break;
			sigma[i+1]=sigma[i];
		}
		sigma[i+1]=x;
	}
}


double Foam::mathFunctionsWENO::realEigenValues::sign
(
    const double& a,
    const double& b
) const
{
    // Get magnitude of a 
    double c = std::abs(a);
    
    if (std::signbit(a) != std::signbit(b))
        c = -c;
    
    return c;
}

#endif
