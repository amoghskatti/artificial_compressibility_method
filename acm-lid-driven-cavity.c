#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define grid 128

int main()
{
    double u[grid][grid+1],un[grid][grid+1];//,u_n[grid][grid+1];
    double v[grid+1][grid],vn[grid+1][grid];//,v_n[grid+1][grid];
    double p[grid+1][grid+1],pn[grid+1][grid+1];//,p_n[grid+1][grid+1];
    double T[grid+1][grid+1],Tn[grid+1][grid+1];//,T_n[grid+1][grid+1];
    double m[grid+1][grid+1];
    double Re = 100.0;
    double Pr = 6.2;
    double delta = 0.5;
    double error = 1.0;
    double dx = 1.0/(grid-1), dy=1.0/(grid-1), dt=0.001;
    int i,j,step = 1;

    // initialization of u
    for(i=0;i<=grid-1;i++)
    {
        for(j=0;j<=grid-2;j++)
            u[i][j]=0.0;
        u[i][grid]=1.0;
        u[i][grid-1]=1.0;
    }

    // initialization of v
    for(i=0;i<=grid;i++)
    {
        for(j=0;j<=grid-1;j++)
            v[i][j]=0.0;
    }

    //initialization of p
    for(i=0;i<=grid;i++)
    {
        for(j=0;j<=grid;j++)
            p[i][j]=1.0;
    }
    //initialization of T
    for(i=0;i<=grid;i++)
    {
        for(j=0;j<=grid-2;j++)
            T[i][j] = 0.0;
        T[i][grid-1] = 1.0;
        T[i][grid] = 1.0;
    }

    while(error > 0.00000001)
    {

        //solving x-momentum corrector equation
        for(i=1;i<=grid-2;i++)
        {
            for(j=1;j<=grid-1;j++)
            {
                un[i][j] = u[i][j]-(dt)*( (((u[i+1][j]*u[i+1][j])-(u[i-1][j]*u[i-1][j]))/(2*dx))
                                       + (((u[i][j+1]+u[i][j])*(v[i][j]+v[i+1][j]) - (u[i][j]+u[i][j-1])*(v[i][j-1]+v[i+1][j-1]))/(4.0*dy))
                                       + ((p[i+1][j]-p[i][j])/dx)
                                       - (1.0/Re)*(((u[i+1][j]-2.0*u[i][j]+u[i-1][j])/(dx*dx))+((u[i][j+1]-2.0*u[i][j]+u[i][j-1])/(dy*dy))) );

            }
        }

        //x-momentum boundary conditions
        for(i=0;i<=grid-1;i++)
        {
            un[i][0]= -un[i][1];
            un[i][grid] = 2.0-un[i][grid-1];
        }
        for(j=1;j<=grid-1;j++)
        {
            un[0][j] = 0.0;
            un[grid-1][j] = 0.0;
        }


        //solving y-momentum equation
        for(i=1;i<=grid-1;i++)
        {
            for(j=1;j<=grid-2;j++)
            {
                vn[i][j] = v[i][j]-(dt)*( (((v[i][j+1]*v[i][j+1])-(v[i][j-1]*v[i][j-1]))/(2.0*dy))
                                       + (((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1]) - (v[i][j]+v[i-1][j])*(u[i-1][j]+u[i-1][j+1]))/(4.0*dx))
                                       + ((p[i][j+1]-p[i][j])/dy)
                                       - (1.0/Re)*(((v[i][j+1]-2.0*v[i][j]+v[i][j-1])/(dy*dy))+((v[i+1][j]-2.0*v[i][j]+v[i-1][j])/(dx*dx))) );
            }
        }

        //y-momentum boundary conditions
        for(j=0;j<=grid-1;j++)
        {
            vn[0][j]= -vn[1][j];
            vn[grid][j] = -vn[grid-1][j];
        }
        for(i=1;i<=grid-1;i++)
        {
            vn[i][0] = 0.0;
            vn[i][grid-1] = 0.0;
        }

        //continuity equation - solve for p
        for(i=1;i<=grid-1;i++)
        {
            for(j=1;j<=grid-1;j++)
            {
                pn[i][j] = p[i][j] - delta*(dt)*( ((un[i][j] - un[i-1][j])/dx) + ((vn[i][j] - vn[i][j-1])/dy) );
            }
        }

        //p boundary conditions
        for(i=1;i<=grid-1;i++)
        {
            pn[i][0] = pn[i][1];
            pn[i][grid] = pn[i][grid-1];
        }
        for(j=0;j<=grid;j++)
        {
            pn[0][j] = pn[1][j];
            pn[grid][j] = pn[grid-1][j];
        }

        for(i=1;i<=grid-1;i++)
        {
            for(j=1;j<=grid-1;j++)
            {
                Tn[i][j] = T[i][j] + dt*( (1.0/(Re*Pr))*( (T[i+1][j]-2.0*T[i][j]+T[i-1][j])/(dx*dx) + (T[i][j+1]-2.0*T[i][j]+T[i][j-1])/(dy*dy) )
                           - ((u[i-1][j]+u[i][j])/2.0)*((T[i+1][j]-T[i-1][j])/(2.0*dx))
                           - ((v[i][j-1]+v[i][j])/2.0)*((T[i][j+1]-T[i][j-1])/(2.0*dy)) );
            }
        }
        for(j=1;j<=grid-1;j++)
        {
            Tn[0][j] = Tn[1][j];
            Tn[grid][j] = Tn[grid-1][j];
        }
        for(i=0;i<=grid;i++)
        {
            Tn[i][0] = -Tn[i][1];
            Tn[i][grid] = 2.0-Tn[i][grid-1];
        }

        //error
        error = 0.0;
        for(i=1;i<=grid-1;i++)
        {
            for(j=1;j<=grid-1;j++)
            {
                m[i][j] = ( ((un[i][j]-u[i][j])/dt) + ((vn[i][j]-v[i][j])/dt) ) + (Tn[i][j] - T[i][j])/dt;
                error = error + fabs(m[i][j]);
            }
        }

        if(step%1000 == 1)
            printf("Error is %f for step %d\n",error,step);

        //next iteration
        for(i=0;i<=grid-1;i++)
        {
            for(j=0;j<=grid;j++)
                u[i][j] = un[i][j];
        }

        for(i=0;i<=grid;i++)
        {
            for(j=0;j<=grid-1;j++)
                v[i][j] = vn[i][j];
        }
        for(i=0;i<=grid;i++)
        {
            for(j=0;j<=grid;j++)
                p[i][j] = pn[i][j];
        }
        for(i=0;i<=grid;i++)
        {
            for(j=0;j<=grid;j++)
                T[i][j] = Tn[i][j];
        }
        step++;

    }

    double uc[grid][grid], vc[grid][grid], pc[grid][grid], Tc[grid][grid];

    for(i=0;i<=grid-1;i++)
    {
        for(j=0;j<=grid-1;j++)
        {
            uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
            vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
            pc[i][j] = 0.25*(p[i][j]+p[i][j+1]+p[i+1][j]+p[i+1][j+1]);
            Tc[i][j] = 0.25*(T[i][j]+T[i][j+1]+T[i+1][j]+T[i+1][j+1]);
        }
    }

    double Nu[grid], NuAv, sumNu = 0.0;
    for(i=0;i<=grid-1;i++)
    {
        Nu[i] = (Tc[i][grid-1]-Tc[i][grid-3])/(2.0*dy);
        sumNu+=Nu[i];
    }
    NuAv = sumNu/grid;
    printf("\n%5.8lf\n",NuAv);

    // OUTPUT DATA
	FILE *fout2, *fout3, *fout4, *fout5, *fout6;
	fout2 = fopen("UVPT100.plt","w+t");
	fout3 = fopen("Central_U100.plt","w+t");
	fout4 = fopen("Central_T100.plt","w+t");
	fout5 = fopen("Nu100_Plot.plt","w+t");
    fout6 = fopen("Average_Nu_Plot100.plt","w+t");

	if ( fout2 == NULL )
	{
        printf("\nERROR when opening file\n");
        fclose( fout2 );
	}

    else
	{
        fprintf( fout2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\",\"T\"\n");
        fprintf( fout2, "ZONE  F=POINT\n");
        fprintf( fout2, "I=%d, J=%d\n", grid, grid );

        for ( j = 0 ; j < (grid) ; j++ )
        {
            for ( i = 0 ; i < (grid) ; i++ )
            {
                double xpos, ypos;
                xpos = i*dx;
                ypos = j*dy;

                fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, uc[i][j], vc[i][j], pc[i][j], Tc[i][j] );
            }
        }
	}

	fclose( fout2 );

	// CENTRAL --U
    fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
    fprintf(fout3, "ZONE F=POINT\n");
    fprintf(fout3, "I=%d\n", grid );

    for ( j = 0 ; j < grid ; j++ )
    {
        double ypos;
        ypos = (double) j*dy;

        fprintf( fout3, "%5.8lf\t%5.8lf\n", (uc[grid/2][j] + uc[(grid/2)+1][j])/(2.), ypos );
    }

    fclose(fout3);

    fprintf(fout4, "VARIABLES=\"T\",\"Y\"\n");
    fprintf(fout4, "ZONE F=POINT\n");
    fprintf(fout4, "I=%d\n", grid );
    for ( j = 0 ; j < grid ; j++ )
    {
        double ypos;
        ypos = (double) j*dy;

        fprintf( fout4, "%5.8lf\t%5.8lf\n", (Tc[grid/2][j] + Tc[(grid/2)+1][j])/(2.), ypos );
    }
    fclose(fout4);

    fprintf(fout5, "VARIABLES=\"Nu\",\"X\"\n");
    fprintf(fout5, "ZONE F=POINT\n");
    fprintf(fout5, "I=%d\n", grid );
    for ( i = 0 ; i < grid ; i++ )
    {
        double xpos;
        xpos = (double) i*dy;

        fprintf( fout5, "%5.8lf\t%5.8lf\n", Nu[i], xpos );
    }
    fclose(fout5);

    fprintf(fout6, "VARIABLES=\"NuAv\"\n");
    fprintf(fout6, "ZONE F=POINT\n");
    fprintf(fout6, "I=%d\n", 1);
    fprintf(fout6, "%5.8lf\n", NuAv);
    fclose(fout6);
}
