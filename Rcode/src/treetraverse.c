#include <R.h>

void postorder_tree_traverse(int * anc, int * desc, int * node, int * nnode, int * order);
void preorder_tree_traverse(int * anc, int * desc, int * node, int * nnode, int * order);
void setrecursivesequence(int * anc, int * desc, int * node, int * ne, int * downseq, int * lastvisit);
void rootward(int * anc, int * desc, int * node, int * nnode, int * ndorder);
void tipward(int * anc, int * desc, int * node, int * nnode, int * ndorder);
void recursivesequence(int * anc, int * desc, int * node, int * ne, int * downseq, int * lastvisit);
void setpolartreecoords(int * anc, int * desc, int * ndorder, int * ntip, int * rootnd, int * nnode, double * ths, double * theta, double * root);
void setphylotreecoords(int * anc, int * desc, int * ndorder, double * bl, int * ntip, int * rootnd, int * nnode, double * bar, double * xy, double * root);
void fetchmrca(int * anc, int * desc, int * root, int * ne, int * npair, int * t1, int * t2, int * ret);

static int zzz;

void postorder_tree_traverse(int * anc, int * desc, int * node, int * nnode, int * ndorder)
{
	zzz = 0;
	rootward(anc, desc, node, nnode, ndorder);
}

void preorder_tree_traverse(int * anc, int * desc, int * node, int * nnode, int * ndorder)
{
	zzz = 0;
	tipward(anc, desc, node, nnode, ndorder);
}

void setrecursivesequence(int * anc, int * desc, int * node, int * ne, int * downseq, int * lastvisit)
{
	zzz = 0;
	recursivesequence(anc,desc,node,ne,downseq,lastvisit);
}

void rootward(int * anc, int * desc, int * node, int * nnode, int * ndorder)
{
	int i, d = 0;
	int * children;
	children = Calloc(2, int);
	for (i = 0; i < (*nnode-1); i++)
	{
		if (anc[i] == *node)
		{
			children[d] = desc[i];
			d++;
		}
		if (d == 2) break;
	}
	
	if (children[0] != 0 && children[1] != 0)
	{
		int * child;
		child = Calloc(1,int); 
		for (i = 0; i < 2; i++)
		{
			*child = children[i];
			rootward(anc, desc, child, nnode, ndorder);
			//Free(child);
		}
		Free(child);
	}
	ndorder[zzz] = *node; zzz++;
	//Rprintf("%d\n", zzz);
	Free(children);
}

void tipward(int * anc, int * desc, int * node, int * nnode, int * ndorder)
{
	ndorder[zzz] = *node; zzz++;
	int i, d = 0;
	int * children;
	children = Calloc(2, int);
	for (i = 0; i < (*nnode-1); i++)
	{
		if (anc[i] == *node)
		{
			children[d] = desc[i];
			d++;
		}
		if (d == 2) break;
	}
	
	if (children[0] != 0 && children[1] != 0)
	{
		int * child;
		child = Calloc(1,int); 
		for (i = 0; i < 2; i++)
		{
			*child = children[i];
			tipward(anc, desc, child, nnode, ndorder);
			//Free(child);
		}
		Free(child);
	}
	//Rprintf("%d\n", zzz);
	Free(children);
}

void recursivesequence(int * anc, int * desc, int * node, int * ne, int * downseq, int * lastvisit)
{
	downseq[zzz] = *node; zzz++; 
	int i, d = 0;
	int * children;
	children = Calloc(2, int);
	for (i = 0; i < *ne; i++)
	{
		if (anc[i] == *node)
		{
			children[d] = desc[i];
			d++;
		}
		if (d == 2) break;
	}
	if (children[0] != 0 && children[1] != 0)
	{
		int * child;
		child = Calloc(1, int);
		for (i = 0; i < 2; i++)
		{
			*child = children[i];
			recursivesequence(anc,desc,child,ne,downseq,lastvisit);
		}
		Free(child);
	}
	for (i = 0; i < *ne+1; i++) if (downseq[i] == 0) break;
	lastvisit[*node-1] = downseq[i-1];  
	Free(children);	
}

void setpolartreecoords(int * anc, int * desc, int * ndorder, int * ntip, int * rootnd, int * nnode, double * ths, double * theta, double * root)
{
	int i, j, d, ne = *nnode - 1;
	int * ci;
	ci = Calloc(2, int);
	
	for (i = 0; i < *nnode; i++)
	{
		d = 0;
		if (ndorder[i] <= *ntip)
		{
			for (j = 0; j < ne; j++)
			{
				if (desc[j] == ndorder[i])
				{
					theta[j + ne * 0] = (ndorder[i] - 1) * (*ths);
					theta[j + ne * 1] = (ndorder[i] - 1) * (*ths);
					theta[j + ne * 2] = (ndorder[i] - 1) * (*ths); 
					break;
				}
			}
		}
		else
		{
			for (j = 0; j < ne; j++)
			{
				if (anc[j] == ndorder[i])
				{
					ci[d] = j;
					d++;
				}
				if (d == 2) break;
			}
			if (ndorder[i] == *rootnd)
			{
				root[0] = theta[ci[0] + ne * 0]/2. + theta[ci[1] + ne * 0]/2.;
				root[1] = theta[ci[0] + ne * 0];
				root[2] = theta[ci[1] + ne * 0];
				continue;
			}
			for (j = 0; j < ne; j++) if (desc[j] == ndorder[i]) break;
			theta[j + ne * 0] = theta[ci[0] + ne * 0]/2. + theta[ci[1] + ne * 0]/2.;
			theta[j + ne * 1] = theta[ci[0] + ne * 0];
			theta[j + ne * 2] = theta[ci[1] + ne * 0];
		}
	}
	Free(ci);
}

void setphylotreecoords(int * anc, int * desc, int * ndorder, double * bl, int * ntip, int * rootnd, int * nnode, double * bar, double * xy, double * root)
{
	int i, j, d, ne = *nnode - 1;
	double dy = 1./(*ntip);
	int * ci;
	ci = Calloc(2, int);
	
	for (i = 0; i < *nnode; i++)
	{
		d = 0;
		if (ndorder[i] <= *ntip)
		{
			for (j = 0; j < ne; j++)
			{
				if (desc[j] == ndorder[i])
				{
					xy[j + ne * 0] = 1. - bl[j];
					xy[j + ne * 1] = (ndorder[i] - 1) * dy;
					xy[j + ne * 2] = 1.;
					xy[j + ne * 3] = (ndorder[i] - 1) * dy;
					
					bar[j + ne * 0] = 0.;
					bar[j + ne * 1] = 0.;
					bar[j + ne * 2] = 0.;
					bar[j + ne * 3] = 0.;

					break;
				}
			}
		}
		else
		{
			for (j = 0; j < ne; j++)
			{
				if (anc[j] == ndorder[i])
				{
					ci[d] = j;
					d++;
				}
				if (d == 2) break;
			}
			if (ndorder[i] == *rootnd)
			{
				root[0] = xy[ci[0] + ne * 0];
				root[1] = xy[ci[0] + ne * 3];
				root[2] = xy[ci[0] + ne * 0];
				root[3] = xy[ci[1] + ne * 3];
				continue;
			}
			for (j = 0; j < ne; j++) if (desc[j] == ndorder[i]) break;
			
			xy[j + ne * 0] = xy[ci[0] + ne * 0] - bl[j];
			xy[j + ne * 1] = xy[ci[0] + ne * 3]/2. + xy[ci[1] + ne * 3]/2.;
			xy[j + ne * 2] = xy[ci[0] + ne * 0];
			xy[j + ne * 3] = xy[ci[0] + ne * 3]/2. + xy[ci[1] + ne * 3]/2.;
					
			bar[j + ne * 0] = xy[ci[0] + ne * 0];
			bar[j + ne * 1] = xy[ci[0] + ne * 3];
			bar[j + ne * 2] = xy[ci[1] + ne * 0];
			bar[j + ne * 3] = xy[ci[1] + ne * 3];
		}
	}
	Free(ci);
}

void fetchmrca(int * anc, int * desc, int * root, int * ne, int * npair, int * t1, int * t2, int * ret)
{	
	int i,j,k;

	int cnt, node, mrca;
	int * path;
	
	for (k=0; k < *npair; k++)
	{
		if (t2[k] == 0)
		{
			ret[k] = t1[k];
			continue;
		}
		path = Calloc(*ne, int);
		cnt = 0; mrca = 0;
		node = t1[k];
		while (node != *root)
		{
			for (i = 0; i < *ne; i++)
			{
				if (desc[i] == node)
				{
					node = anc[i];
					path[cnt] = node;
					cnt++;
					break;
				}
			}
		}
		
		node = t2[k];
		while (node != *root)
		{
			for (i = 0; i < *ne; i++)
			{
				if (desc[i] == node)
				{
					node = anc[i];
					for (j = 0; j < *ne; j++)
					{
						if (node == path[j])
						{
							mrca = 1;
							break;
						}
					}
				}
				if (mrca == 1) break;
			}
			if (mrca == 1) break;
		}
	
		if (mrca == 1)
		{
			ret[k] = node;
		}
		else
		{
			ret[k] = *root;
		}
		Free(path);
	}
}
