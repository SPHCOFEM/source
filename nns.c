#include "header.h"

#undef __FUNC__
#define __FUNC__ "nns"

void nns(int n_s,int opt,int *nn_s,double *x_s,double *h_s,int **neig_s)
{
	int i,j,count_neig=0;
  	double aux1,aux2;

  	fprintf(stdout,
	      	//"\n"
			"nns: looking for nearest neighbours ..."
	      	"\n"
	      	"\n");

	/* setting nearest neighbour counts vector */
	memset(nn_s,0,n_s*sizeof(int));

	for (i=0;i<n_s;i++)
    {
		for (j=i+1;j<n_s;j++)
	  	{
	    	aux1=norm(x_s[i]-x_s[j],x_s[i+n_s]-x_s[j+n_s],x_s[i+2*n_s]-x_s[j+2*n_s]);
	    	aux2=(h_s[i]+h_s[j])/2.0;

	    	if (aux1<2.0*opt*aux2)
			{
		  		nn_s[i]++;
		  		count_neig++;
			}
	    }
    }

    /* nearest neighbours exist */
	if (count_neig>0)
	{
		*neig_s=createMemMore(int,count_neig);
	}
	else
	{
		fprintf(stdout,
			"ERROR: no nearest neighbours!"
			"\n"
			"\n");
			exit(0);	  
	}

    /* reset counter */
	count_neig=0;

    for (i=0;i<n_s;i++)
	{
		for (j=i+1;j<n_s;j++)
	    {
	    	aux1=norm(x_s[i]-x_s[j],x_s[i+n_s]-x_s[j+n_s],x_s[i+2*n_s]-x_s[j+2*n_s]);
	      	aux2=(h_s[i]+h_s[j])/2.0;

	      	if (aux1<2.0*opt*aux2)
			{
		  		(*neig_s)[count_neig]=j;
		  		count_neig++;
			}
	    }
	}

	/* print records 
    count_neig=0;
    for (i=0;i<n_s;i++)
	{
		fprintf(stdout,"Particle %d has %d neighbours : ",i,nn_s[i]);
		for (j=0;j<nn_s[i];j++)
		{
			fprintf(stdout,"%d ",(*neig_s)[count_neig+j]);
		}
		fprintf(stdout,"\n");
		count_neig=count_neig+nn_s[i];
	}
	exit(0);
	*/
}
