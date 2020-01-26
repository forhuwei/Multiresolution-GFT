/*************************************************************************
 * Hochbaum's Pseudo-flow (HPF) Algorithm Matlab implementation          *
 * ************************************************************          *
 * The HPF algorithm for finding Minimum-cut in a graph is described in: *                                        *
 * [1] D.S. Hochbaum, "The Pseudoflow algorithm: A new algorithm for the *
 * maximum flow problem", Operations Research, 58(4):992-1009,2008.      *
 *                                                                       *
 * The algorithm was found to be fast in theory (see the above paper)    *
 * and in practice (see:                                                 *
 * [2] D.S. Hochbaum and B. Chandran, "A Computational Study of the      *
 * Pseudoflow and Push-relabel Algorithms for the Maximum Flow Problem,  *
 * Operations Research, 57(2):358-376, 2009.                             *
 *
 * and                                                                   *
 *                                                                       *
 * [3] B. Fishbain, D.S. Hochbaum, S. Mueller, "Competitive Analysis of  *
 * Minimum-Cut Maximum Flow Algorithms in Vision Problems,               *
 * arXiv:1007.4531v2 [cs.CV]                                             *
 *                                                                       *
 * Usage: Within Matlab environment:                                     *
 * [value,cut] = hpf(sim_mat,source,sink);                               *
 *                                                                       *
 * INPUTS                                                                *
 * ******                                                                *
 * sim_mat - similarity matrix - a_{i,j} is the capacity of the arc (i,j)*
 *           a_{i,j} are non-negatives; the self-similarities (diagonal  *
 *           values) are zero.                                           *
 *           the sim_mat should be sparse (see Matlab's help)            *
 * source - The numeric label of the source node                         *
 * sink   - The numeric label of the sink node                           *
 *                                                                       *
 * OUTPUTS                                                               *
 * *******                                                               *
 * value - the capacity of the cut                                       *
 * cut   - the source set (see [1]), where x_i = 1, if i \in S ; 0 o/w   *
 *                                                                       *
 * Set-up                                                                *
 * ******                                                                *
 * Uncompress the MatlabHPF.zip file into the Matlab's working directory *
 * The zip file contains the following files:                            *
 * hpf.c - source code                                                   *
 * hpf.m - Matlab's help file                                            *
 * hpf.mexmaci - The compiled code for Mac OS 10.0.5 (Intel)/ Matlab     *
 *               7.6.0.324 (R2008a).                                     *
 * hpf.mexw32  - The compiled code for Windows 7 / Matlab 7.11.0.584     *
 *               (R2010b).                                               *
 * demo_general - Short Matlab code that generates small network and     *
 *                computes the minimum flow                              *
 * demo_vision - Short Matlab code that loads a Multiview reconstruction *
 *               vision problem (see: [3]) and computes its minimum cut. *
 * gargoyle-smal.mat - The vision problem.                               *
 *                                                                       *
 * When using this code, please cite:                                    *
 * References [1], [2] and [3] above and:                                *
 * B. Fishbain and D.S. Hochbaum, "Hochbaum's Pseudo-flow Matlab         *
 * implementation", http://riot.ieor.berkeley.edu/riot/Applications/     *
 * Pseudoflow/maxflow.html                                               *
 *************************************************************************/

#include "mex.h"
#include "matrix.h"
#include <stdio.h>
//#include <sys/time.h>
//#include <sys/resource.h>
#include <stdlib.h>
//#include <unistd.h>

/*************************************************************************
Definitions
*************************************************************************/
#define  MAX_LEVELS  300
#define VERSION 3.3

typedef unsigned int uint;
typedef long int lint;
typedef long long int llint;
typedef unsigned long long int ullint;

struct node;

typedef struct arc 
	{
		struct node *from;
		struct node *to;
		uint flow;
		uint capacity;
		uint direction;
	} Arc;

typedef struct node 
	{
		uint visited;
		uint numAdjacent;
		uint number;
		uint label;
		int excess;
		struct node *parent;
		struct node *childList;
		struct node *nextScan;
		uint numOutOfTree;
		Arc **outOfTree;
		uint nextArc;
		Arc *arcToParent;
		struct node *next;
	} Node;

typedef struct root 
	{
		Node *start;
		Node *end;
	} Root;

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

/* Input Arguments */
#define	AFFIN_MAT	prhs[0]
#define SOURCE		prhs[1]
#define SINK		prhs[2]

/* Output Arguments */
#define	CUT_VAL	plhs[0]
#define	SEGMENTS plhs[1]

/*************************************************************************
Global variables
*************************************************************************/
static uint numNodes = 0;
static uint numArcs = 0;
static uint source = 0;
static uint sink = 0;
static uint highestStrongLabel = 1;

static Node *nodesList = NULL;
static Root *strongRoots = NULL;
static uint *labelCount = NULL;
static Arc *arcList = NULL;
static uint lowestPositiveExcessNode = 0;

bool isNonNegative(const mxArray* affinityMat)
/*************************************************************************
isNonNegative
*************************************************************************/
{
  double  *pr, *pi;
  mwIndex  *ir, *jc;
  mwSize      col, total=0;
  mwIndex   starting_row_index, stopping_row_index, current_row_index;
  mwSize      n;
  
  /* Get the starting positions of all four data arrays. */ 
  pr = mxGetPr(affinityMat);
  pi = mxGetPi(affinityMat);
  ir = mxGetIr(affinityMat);
  jc = mxGetJc(affinityMat);
  
  /* Display the nonzero elements of the sparse array. */ 
  n = mxGetN(affinityMat);
  for (col=0; col<n; col++)  { 
    starting_row_index = jc[col]; 
    stopping_row_index = jc[col+1]; 
    if (starting_row_index == stopping_row_index)
      continue;
    else 
	{
      for (current_row_index = starting_row_index; current_row_index < stopping_row_index; current_row_index++)  
	  {
		  if (pr[total++] < 0) 
		  {
			  return(false);			
		  }
	  }
	}
  }
  return(true);
}

bool isDiagonalZero(const mxArray* affinityMat)
/*************************************************************************
isDiagonalZero
*************************************************************************/
{
double  *pr, *pi;
  mwIndex  *ir, *jc;
  mwSize      col, total=0;
  mwIndex   starting_row_index, stopping_row_index, current_row_index;
  mwSize      n;
  
  /* Get the starting positions of all four data arrays. */ 
  pr = mxGetPr(affinityMat);
  pi = mxGetPi(affinityMat);
  ir = mxGetIr(affinityMat);
  jc = mxGetJc(affinityMat);
  
  /* Display the nonzero elements of the sparse array. */ 
  n = mxGetN(affinityMat);
  for (col=0; col<n; col++)  { 
    starting_row_index = jc[col]; 
    stopping_row_index = jc[col+1]; 
    if (starting_row_index == stopping_row_index)
      continue;
    else {
      for (current_row_index = starting_row_index; current_row_index < stopping_row_index; current_row_index++)  
	  {
		  if ((ir[current_row_index] == col) && (pr[total++] != 0))
		  {
			  return(false);
		  }
	  }
	}
  }
  return(true);
}

static void createOutOfTree (Node *nd)
{
/*************************************************************************
createOutOfTree
*************************************************************************/
	if (nd->numAdjacent)
	{
		if ((nd->outOfTree = (Arc **) mxMalloc (nd->numAdjacent * sizeof (Arc *))) == NULL)
		{
			mexErrMsgTxt("Out of memory\n");
		}
	}
}

static void initializeArc (Arc *ac)
{
/*************************************************************************
initializeArc
*************************************************************************/
	ac->from = NULL;
	ac->to = NULL;
	ac->capacity = 0;
	ac->flow = 0;
	ac->direction = 1;
}

static void liftAll (Node *rootNode) 
{
/*************************************************************************
liftAll
*************************************************************************/
	Node *temp, *current=rootNode;
	
	current->nextScan = current->childList;
	
	-- labelCount[current->label];
	current->label = numNodes;	
	
	for ( ; (current); current = current->parent)
	{
		while (current->nextScan) 
		{
			temp = current->nextScan;
			current->nextScan = current->nextScan->next;
			current = temp;
			current->nextScan = current->childList;
			
			-- labelCount[current->label];
			current->label = numNodes;	
		}
	}
}
static void addOutOfTreeNode (Node *n, Arc *out) 
{
/*************************************************************************
addOutOfTreeNode
*************************************************************************/
	n->outOfTree[n->numOutOfTree] = out;
	++ n->numOutOfTree;
}

static void addToStrongBucket (Node *newRoot, Root *rootBucket) 
{
/*************************************************************************
addToStrongBucket
*************************************************************************/
	if (rootBucket->start)
	{
		rootBucket->end->next = newRoot;
		rootBucket->end = newRoot;
		newRoot->next = NULL;
	}
	else
	{
		rootBucket->start = newRoot;
		rootBucket->end = newRoot;
		newRoot->next = NULL;
	}
}

static __inline int addRelationship (Node *newParent, Node *child) 
{
/*************************************************************************
addRelationship
*************************************************************************/
	child->parent = newParent;
	child->next = newParent->childList;
	newParent->childList = child;
	
	return 0;
}

static __inline void breakRelationship (Node *oldParent, Node *child) 
{
/*************************************************************************
breakRelationship
*************************************************************************/
	Node *current;
	
	child->parent = NULL;
	
	if (oldParent->childList == child) 
	{
		oldParent->childList = child->next;
		child->next = NULL;
		return;
	}
	
	for (current = oldParent->childList; (current->next != child); current = current->next);
	
	current->next = child->next;
	child->next = NULL;
}

static void merge (Node *parent, Node *child, Arc *newArc) 
{
/*************************************************************************
merge
*************************************************************************/
	Arc *oldArc;
	Node *current = child, *oldParent, *newParent = parent;
	
	while (current->parent) 
	{
		oldArc = current->arcToParent;
		current->arcToParent = newArc;
		oldParent = current->parent;
		breakRelationship (oldParent, current);
		addRelationship (newParent, current);
		newParent = current;
		current = oldParent;
		newArc = oldArc;
		newArc->direction = 1 - newArc->direction;
	}
	
	current->arcToParent = newArc;
	addRelationship (newParent, current);
}


static __inline void pushUpward (Arc *currentArc, Node *child, Node *parent, const uint resCap) 
{
/*************************************************************************
pushUpward
*************************************************************************/
	if ((int)resCap >= child->excess) 
	{
		parent->excess += child->excess;
		currentArc->flow += child->excess;
		child->excess = 0;
		return;
	}
	
	currentArc->direction = 0;
	parent->excess += resCap;
	child->excess -= resCap;
	currentArc->flow = currentArc->capacity;
	parent->outOfTree[parent->numOutOfTree] = currentArc;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);
	
	addToStrongBucket (child, &strongRoots[child->label]);
}


static __inline void pushDownward (Arc *currentArc, Node *child, Node *parent, uint flow) 
{
/*************************************************************************
pushDownward
*************************************************************************/
	if ((int)flow >= child->excess) 
	{
		parent->excess += child->excess;
		currentArc->flow -= child->excess;
		child->excess = 0;
		return;
	}
	
	currentArc->direction = 1;
	child->excess -= flow;
	parent->excess += flow;
	currentArc->flow = 0;
	parent->outOfTree[parent->numOutOfTree] = currentArc;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);
	
	addToStrongBucket (child, &strongRoots[child->label]);
}

static void pushExcess (Node *strongRoot) 
{
/*************************************************************************
pushExcess
*************************************************************************/
	Node *current, *parent;
	Arc *arcToParent;
	int prevEx=1;
	
	for (current = strongRoot; (current->excess && current->parent); current = parent) 
	{
		parent = current->parent;
		prevEx = parent->excess;
		
		arcToParent = current->arcToParent;
		
		if (arcToParent->direction)
		{
			pushUpward (arcToParent, current, parent, (arcToParent->capacity - arcToParent->flow)); 
		}
		else
		{
			pushDownward (arcToParent, current, parent, arcToParent->flow); 
		}
	}
	
	if ((current->excess > 0) && (prevEx <= 0))
	{
		addToStrongBucket (current, &strongRoots[current->label]);
	}
}


static Arc * findWeakNode (Node *strongNode, Node **weakNode) 
{
/*************************************************************************
findWeakNode
*************************************************************************/
	uint i, size;
	Arc *out;
	
	size = strongNode->numOutOfTree;
	
	for (i=strongNode->nextArc; i<size; ++i) 
	{	
		if (strongNode->outOfTree[i]->to->label == (highestStrongLabel-1)) 
		{
			strongNode->nextArc = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->to;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		} else if (strongNode->outOfTree[i]->from->label == (highestStrongLabel-1)) {
			strongNode->nextArc = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->from;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		}
	}
	
	strongNode->nextArc = strongNode->numOutOfTree;
	
	return NULL;
}


static void checkChildren (Node *curNode) 
{
/*************************************************************************
checkChildren
*************************************************************************/
	for ( ; (curNode->nextScan); curNode->nextScan = curNode->nextScan->next)
	{
		if (curNode->nextScan->label == curNode->label)
		{
			return;
		}
		
	}	
	
	-- labelCount[curNode->label];
	++	curNode->label;
	++ labelCount[curNode->label];
	
	curNode->nextArc = 0;
}


static void simpleInitialization (void) 
{
/*************************************************************************
simpleInitialization
*************************************************************************/
	uint i, size;
	Arc *tempArc;
	
	size = nodesList[source-1].numOutOfTree;
	for (i=0; i<size; ++i) // Saturating source adjacent nodes
	{
		tempArc = nodesList[source-1].outOfTree[i];
		tempArc->flow = tempArc->capacity;
		tempArc->to->excess += tempArc->capacity;
	}
	
	size = nodesList[sink-1].numOutOfTree;
	for (i=0; i<size; ++i) // Pushing maximum flow on sink adjacent nodes
	{
		tempArc = nodesList[sink-1].outOfTree[i];
		tempArc->flow = tempArc->capacity;
		tempArc->from->excess -= tempArc->capacity;
	}
	
	nodesList[source-1].excess = 0; // zeroing source excess 
	nodesList[sink-1].excess = 0;	// zeroing sink excess
	
	for (i=0; i<numNodes; ++i) 
	{
		if (nodesList[i].excess > 0) 
		{
		    nodesList[i].label = 1;
			++ labelCount[1];
			
			addToStrongBucket (&nodesList[i], &strongRoots[1]);
		}
	}
	
	nodesList[source-1].label = numNodes;	// Set the source label to n
	nodesList[sink-1].label = 0;			// set the sink label to 0
	labelCount[0] = (numNodes - 2) - labelCount[1];
}


static Node* getHighestStrongRoot (void) 
{
/*************************************************************************
getHighestStrongRoot
*************************************************************************/
	uint i;
	Node *strongRoot;
	
	for (i=highestStrongLabel; i>0; --i) 
	{
		if (strongRoots[i].start)  
		{
			highestStrongLabel = i;
			if (labelCount[i-1]) 
			{
				strongRoot = strongRoots[i].start;
				strongRoots[i].start = strongRoot->next;
				strongRoot->next = NULL;
				return strongRoot;				
			}
			
			while (strongRoots[i].start) 
			{
				strongRoot = strongRoots[i].start;
				strongRoots[i].start = strongRoot->next;
				liftAll (strongRoot);
			}
		}
	}
	
	if (!strongRoots[0].start) 
	{
		return NULL;
	}
	
	while (strongRoots[0].start) 
	{
		strongRoot = strongRoots[0].start;
		strongRoots[0].start = strongRoot->next;
		strongRoot->label = 1;
		-- labelCount[0];
		++ labelCount[1];
		
		addToStrongBucket (strongRoot, &strongRoots[strongRoot->label]);		
	}	
	
	highestStrongLabel = 1;
	
	strongRoot = strongRoots[1].start;
	strongRoots[1].start = strongRoot->next;
	strongRoot->next = NULL;
	
	return strongRoot;	
}



static void initializeRoot (Root *rt) 
{
/*************************************************************************
initializeRoot
*************************************************************************/
	rt->start = NULL;
	rt->end = NULL;
}


static void initializeNode (Node *nd, const uint n)
{
/*************************************************************************
initializeNode
*************************************************************************/
	nd->label = 0;
	nd->excess = 0;
	nd->parent = NULL;
	nd->childList = NULL;
	nd->nextScan = NULL;
	nd->nextArc = 0;
	nd->numOutOfTree = 0;
	nd->arcToParent = NULL;
	nd->next = NULL;
	nd->visited = 0;
	nd->numAdjacent = 0;
	nd->number = n;
	nd->outOfTree = NULL;
}
static void freeRoot (Root *rt) 
{
/*************************************************************************
freeRoot
*************************************************************************/
	rt->start = NULL;
	rt->end = NULL;
}
static void freeMemory (void)
{
/*************************************************************************
freeMemory
*************************************************************************/
	uint i;
	
	for (i=0; i<numNodes; ++i)
	{
		freeRoot (&strongRoots[i]);
	}
	
	mxFree(strongRoots);
	
	/*for (i=0; i<numNodes; ++i)
	{
		if (nodesList[i].outOfTree)
		{
			mxFree(nodesList[i].outOfTree);
		}
	}*/
	
	mxFree(nodesList);
	mxFree(labelCount);
	mxFree(arcList);
}
static void processRoot (Node *strongRoot) 
{
/*************************************************************************
processRoot
*************************************************************************/
	Node *temp, *strongNode = strongRoot, *weakNode;
	Arc *out;
	
	strongRoot->nextScan = strongRoot->childList;
	
	if ((out = findWeakNode (strongRoot, &weakNode)))
	{
		merge (weakNode, strongNode, out);
		pushExcess (strongRoot);
		return;
	}
	
	checkChildren (strongRoot);
	
	while (strongNode)
	{
		while (strongNode->nextScan) 
		{
			temp = strongNode->nextScan;
			strongNode->nextScan = strongNode->nextScan->next;
			strongNode = temp;
			strongNode->nextScan = strongNode->childList;
			
			if ((out = findWeakNode (strongNode, &weakNode)))
			{
				merge (weakNode, strongNode, out);
				pushExcess (strongRoot);
				return;
			}
			
			checkChildren (strongNode);
		}
		
		if ((strongNode = strongNode->parent))
		{
			checkChildren (strongNode);
		}
	}
	
	addToStrongBucket (strongRoot, &strongRoots[strongRoot->label]);
	++ highestStrongLabel;
}


static __inline void quickSort (Arc **arr, const uint first, const uint last)
{
/*************************************************************************
quickSort
*************************************************************************/
	int i=0, j, L=first, R=last, beg[MAX_LEVELS], end[MAX_LEVELS], temp=0;
	Arc *swap;
	
	if ((R-L) <= 5)
	{// Bubble sort if 5 elements or less
		for (i=R; (i>L); --i)
		{
			swap = NULL;
			for (j=L; j<i; ++j)
			{
				if (arr[j]->flow < arr[j+1]->flow)
				{
					swap = arr[j];
					arr[j] = arr[j+1];
					arr[j+1] = swap;
				}
			}
			
			if (!swap)
			{
				return;
			}
		}
		
		return;
	}
	
	beg[0]=first; 
	end[0]=last+1;
		
	while (i>=0) 
	{
		L=beg[i]; 
		R=end[i]-1;
		
		if (L<R) 
		{
			swap=arr[L];
			while (L<R) 
			{
				while ((arr[R]->flow >= swap->flow) && (L<R)) 
					R--; 
				
				if (L<R) 
				{
					arr[L]=arr[R];
					L++;
				}
					
				while ((arr[L]->flow <= swap->flow) && (L<R)) 
					L++; 
					
				if (L<R) 
				{
					arr[R]=arr[L]; 
					R--;
				}
			}
				
			arr[L] = swap;
				
			beg[i+1] = L+1; 
			end[i+1] = end[i]; 
			end[i++] = L;
				
			if ((end[i]-beg[i]) > (end[i-1]-beg[i-1])) 
			{
				temp=beg[i]; 
				beg[i]=beg[i-1]; 
				beg[i-1]=temp;
				temp=end[i]; 
				end[i]=end[i-1]; 
				end[i-1]=temp;
			}
		}
		else 
		{
			i--; 
		}
	}
}

static __inline void sort (Node * current)
{
/*************************************************************************
sort
*************************************************************************/
	if (current->numOutOfTree > 1)
	{
		quickSort (current->outOfTree, 0, (current->numOutOfTree-1));
	}
}

static __inline void minisort (Node *current) 
{
/*************************************************************************
minisort
*************************************************************************/
	Arc *temp = current->outOfTree[current->nextArc];
	uint i, size = current->numOutOfTree, tempflow = temp->flow;

	for(i=current->nextArc+1; ((i<size) && (tempflow < current->outOfTree[i]->flow)); ++i)
	{
		current->outOfTree[i-1] = current->outOfTree[i];
	}
	current->outOfTree[i-1] = temp;
}


static ullint checkOptimality (const uint gap) 
{
/*************************************************************************
checkOptimality
*************************************************************************/
	uint i, check = 1;
	ullint mincut = 0;
	llint *excess = NULL; 
	
	Arc *tempArc;
	
	excess = (llint *) mxMalloc (numNodes * sizeof (llint));
	if (!excess)
		mexErrMsgTxt("Out of memory\n");
	
	// Pushing depicits from all sink adjacent nodes to the sink 
	for (i=0; i<nodesList[sink-1].numOutOfTree; ++i) 
	{
		tempArc = nodesList[sink-1].outOfTree[i];
		if (tempArc->from->excess < 0) 
		{
			if ((tempArc->from->excess + (int) tempArc->flow)  < 0)
			{
				// Excess is high enough to saturate the arc => Flow on residual arc is zeroed 
				tempArc->from->excess += (int) tempArc->flow;				
				tempArc->flow = 0;
			}
			else
				// Excess is NOT high enough to saturate the arc => Excess is zeroed
			{
				tempArc->flow = (uint) (tempArc->from->excess + (int) tempArc->flow);
				tempArc->from->excess = 0;
			}
		}	
	}
	
	for (i=0; i<numNodes; ++i)
	{
		excess[i] = 0;
	}
	
	for (i=0; i<numArcs; ++i) 
	{
		if ((arcList[i].from->label >= gap) && (arcList[i].to->label < gap))
		{
			mincut += arcList[i].capacity;
		}
		
		if ((arcList[i].flow > arcList[i].capacity) || (arcList[i].flow < 0)) 
		{
			check = 0;
			mexPrintf("Warning - Capacity constraint violated on arc (%d, %d). Flow = %d, capacity = %d\n", 
				   arcList[i].from->number,
				   arcList[i].to->number,
				   arcList[i].flow,
				   arcList[i].capacity);
		}
		excess[arcList[i].from->number - 1] -= arcList[i].flow;
		excess[arcList[i].to->number - 1] += arcList[i].flow;
	}
	
	for (i=0; i<numNodes; i++) 
	{
		if ((i != (source-1)) && (i != (sink-1))) 
		{
			if (excess[i]) 
			{
				check = 0;
				mexPrintf ("Warning - Flow balance constraint violated in node %d. Excess = %lld\n", 
						i+1,
						excess[i]);
			}
		}
	}

	check = 1;
	
	if (excess[sink-1] != mincut) 
	{
		check = 0;
		mexWarnMsgTxt("Flow is not optimal - max flow does not equal min cut!\n");
	}

// 	if (check) 
// 	{
// 		mexPrintf ("Solution checks as optimal. \t Max Flow: \t %lld\n", mincut);
// 	}
	
	mxFree (excess);
	excess = NULL;
	return mincut;
}

static void decompose (Node *excessNode, const uint source, uint *iteration) 
{
/*************************************************************************
decompose
*************************************************************************/
	Node *current = excessNode;
	Arc *tempArc;
	uint bottleneck = excessNode->excess;

	// Find the bottleneck along a path to the source or on a cycle
	for ( ;(current->number != source) && (current->visited < (*iteration)); 
		 current = tempArc->from)
	{
		current->visited = (*iteration);
		tempArc = current->outOfTree[current->nextArc];
		
		if (tempArc->flow < bottleneck)
		{
			bottleneck = tempArc->flow;
		}
	}
	
	if (current->number == source) // the DFS reached the source
	{
		excessNode->excess -= bottleneck;
		current = excessNode;
		
		// Push the excess all the way to the source
		while (current->number != source) 
		{
			tempArc = current->outOfTree[current->nextArc]; // Pick arc going out of node to push excess to 
			tempArc->flow -= bottleneck; // Push back bottleneck excess on this arc
			
			if (tempArc->flow) // If there is still flow on this arc do sort on this current node
			{
				minisort(current);
			} else {
				++ current->nextArc;
			}

			current = tempArc->from;
		}
		return;
	}
		
	++ (*iteration);
	
	bottleneck = current->outOfTree[current->nextArc]->flow;

	while (current->visited < (*iteration))
	{
		current->visited = (*iteration);
		tempArc = current->outOfTree[current->nextArc];
		
		if (tempArc->flow < bottleneck)
		{
			bottleneck = tempArc->flow;
		}
		current = tempArc->from;
	}	
	
	++ (*iteration);
	
	while (current->visited < (*iteration))
	{
		current->visited = (*iteration);
		
		tempArc = current->outOfTree[current->nextArc];
		tempArc->flow -= bottleneck;
		
		if (tempArc->flow) 
		{
			minisort(current);
			current = tempArc->from;
		} else {
			++ current->nextArc;
			current = tempArc->from;
		}
	}
}

static void recoverFlow (const uint gap)
{
/*************************************************************************
recoverFlow
*************************************************************************/
	uint iteration = 1;
	uint i, j;
	Arc *tempArc;
	Node *tempNode;

	Node **nodePtrArray;
	
	if ((nodePtrArray = (Node **) mxMalloc (numNodes * sizeof (Node *))) == NULL)
	{
		mexErrMsgTxt("Out of memory\n");
	}
	
	for (i=0; i < numNodes ; i++)
		nodePtrArray[i] = &nodesList[i];

	// Adding arcs FROM the Source to Source adjacent nodes. 
	for (i=0; i<nodesList[source-1].numOutOfTree; ++i) 
	{
		tempArc = nodesList[source-1].outOfTree[i];
		addOutOfTreeNode (tempArc->to, tempArc);
	}
	
	// Zeroing excess on source and sink nodes
	nodesList[source-1].excess = 0;
	nodesList[sink-1].excess = 0;
	
	for (i=0; i<numNodes; ++i) 
	{
		tempNode = &nodesList[i];
		
		if ((i == (source-1)) || (i == (sink-1)))
		{
			continue;
		}
		
		if (tempNode->label >= gap) //tempNode is in SINK set
		{
			tempNode->nextArc = 0;
			if ((tempNode->parent) && (tempNode->arcToParent->flow))
			{
				addOutOfTreeNode (tempNode->arcToParent->to, tempNode->arcToParent);
			}
			
			for (j=0; j<tempNode->numOutOfTree; ++j) 
			{ // go over all sink-set-node's arcs and look for arc with NO flow
				if (!tempNode->outOfTree[j]->flow) 
				{	// Remove arc with no flow
					-- tempNode->numOutOfTree;
					tempNode->outOfTree[j] = tempNode->outOfTree[tempNode->numOutOfTree];
					-- j;
				}
			}
			
			sort(tempNode);
		}
	}

	for (i=lowestPositiveExcessNode ; i < numNodes ; ++i) 
	{
		tempNode = nodePtrArray[i];
		while (tempNode->excess > 0) 
		{
			++ iteration;
			decompose(tempNode, source, &iteration);
		}
	}
}

static void readData (const mxArray* array_ptr)
/*************************************************************************
readData
*************************************************************************/
{
	double  *pr, *pi;
	mwIndex  *ir, *jc;
	mwIndex   starting_row_index, 
		stopping_row_index, 
		current_row_index;
	mwSize col, 
		n = 0, 
		total=0;
	uint from, 
		to, 
		first=0,
		capacity,
		i = 0;
	

	/* Get the starting positions of all four data arrays. */ 
	pr = mxGetPr(array_ptr);
	pi = mxGetPi(array_ptr);
	ir = mxGetIr(array_ptr);
	jc = mxGetJc(array_ptr);
  
	/* Number of nodes */ 
	numNodes = mxGetN(array_ptr);

	if ((nodesList = (Node *) mxMalloc (numNodes * sizeof (Node))) == NULL)
		mexErrMsgTxt("Could not allocate memory.\n");				
	if ((strongRoots = (Root *) mxMalloc (numNodes * sizeof (Root))) == NULL)
		mexErrMsgTxt("Could not allocate memory.\n");				
	if ((labelCount = (uint *) mxMalloc (numNodes * sizeof (uint))) == NULL)
		mexErrMsgTxt("Could not allocate memory.\n");				

	/* Number of Arcs */ 
	numArcs = *(mxGetJc(array_ptr) + numNodes);
	if ((arcList = (Arc *) mxMalloc (numArcs * sizeof (Arc))) == NULL)
		mexErrMsgTxt("Could not allocate memory.\n");

	/* Initialization */
	for (i=0; i<numNodes; ++i)
	{
		initializeRoot (&strongRoots[i]);
		initializeNode (&nodesList[i], (i+1));
		labelCount[i] = 0;
	}

	for (i=0; i<numArcs; ++i)
	{
		initializeArc (&arcList[i]);
	}
	
	/* Reading Arcs Info */
	for (col=0; (uint)col<numNodes; col++)  
	{
		starting_row_index = jc[col]; 
		stopping_row_index = jc[col+1]; 
		if (starting_row_index != stopping_row_index)
		{
			for (current_row_index = starting_row_index; current_row_index < stopping_row_index; current_row_index++)  
			{
				/* Reading arc parameters */
				from = ir[current_row_index];
				to = col;
				capacity = (uint)pr[total++];

				arcList[first].from = &nodesList[from];
				arcList[first].to = &nodesList[to];
				arcList[first].capacity = capacity;
				
				++ first;
				
				++ nodesList[from].numAdjacent;
				++ nodesList[to].numAdjacent;
			}
		}
	}

	/* Create memory structures */
	for (i=0; i<numNodes; ++i) 
	{
		createOutOfTree (&nodesList[i]);
	}
	
	for (i=0; i<numArcs; i++) 
	{
		to = arcList[i].to->number;
		from = arcList[i].from->number;
		capacity = arcList[i].capacity;
	
		if (!((source == to) || (sink == from) || (from == to))) 
		{
			if ((source == from) && (to == sink)) 
			{
				arcList[i].flow = capacity;
			} else if (to == sink) {
				addOutOfTreeNode (&nodesList[to-1], &arcList[i]);
			} else {
				addOutOfTreeNode (&nodesList[from-1], &arcList[i]);
			}
		}
	}
}


static void pseudoflowPhase1 (void) 
{
/*************************************************************************
pseudoflowPhase1
*************************************************************************/
	Node *strongRoot;
	while ((strongRoot = getHighestStrongRoot ()))  
	{ 
		processRoot (strongRoot);
	}
}

static void displayCut (double* segs, double* cut, const uint gap) 
{
/*************************************************************************
displayCut
*************************************************************************/
	uint i;
	
	for (i=0; i<numNodes; ++i) 
	{
		if (nodesList[i].label >= gap) 
		{
			segs[nodesList[i].number-1] = 1;
		}
	}
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
/*************************************************************************
mexFunction - Main function
*************************************************************************/
{ 
	double *segs = 0,
		*cut = 0;
	uint n = 0, // Number of nodes / vertices in the graph.
		m = 0; // Number of ARCS in the graph (edge is two arcs)
	uint gap;	
	ullint flow = 0;

    /* Check for proper arguments */ 
    if (nrhs != 3) { 
        mexErrMsgTxt("Usage [segmentations, cut] = hpf(Similarity Matrix, source node, sink node)"); 
    } else if (nlhs > 2) {
        mexErrMsgTxt("Too many output arguments."); 
    } else if (mxGetNumberOfDimensions(AFFIN_MAT) > 2) {
        mexErrMsgTxt("Affinity matrix must be a 2D array");
    } else if (mxGetM(AFFIN_MAT) != mxGetN(AFFIN_MAT)) {
        mexErrMsgTxt("Affinity matrix must a square matrix");
    } else if (mxIsEmpty(AFFIN_MAT)) {
        mexErrMsgTxt("Input argument is empty\n");
    } else if (mxIsComplex(AFFIN_MAT)) {
        mexErrMsgTxt("Affinity matrix must consists of real number (non-complex)");
    } else if (!mxIsNumeric(AFFIN_MAT)) {
        mexErrMsgTxt("Affinity matrix must consists of numeric values");
    } else if (!mxIsSparse(AFFIN_MAT)) {
        mexErrMsgTxt("Affinity matrix must be sparse (use sparse(Affinity matrix) function)"); 
    } else if (!isDiagonalZero(AFFIN_MAT)) {
        mexErrMsgTxt("The disimilarity of an elelment to itself Must be Zero");
	} else if (!isNonNegative(AFFIN_MAT)) {
        mexErrMsgTxt("Affinity matrix must consists of non-negative values");
	}

	n = mxGetN(AFFIN_MAT);

	source = (uint)mxGetScalar(SOURCE);	
	sink = (uint)mxGetScalar(SINK);
	
	if ((source > n) || (source < 0) || (sink > n) || (sink < 0))
		mexErrMsgTxt("Error assigning source or sink nodes");

	/* Create a matrix for the return argument */ 
    SEGMENTS = mxCreateDoubleMatrix(1, n, mxREAL); 
    CUT_VAL  = mxCreateDoubleMatrix(1, 1, mxREAL); 

    /* Assign pointers to the output parameters */ 
    segs = mxGetPr(SEGMENTS);
    cut  = mxGetPr(CUT_VAL);
    
	readData(AFFIN_MAT);
	simpleInitialization ();
	pseudoflowPhase1 ();
	
	gap = numNodes;
	recoverFlow(gap);
	flow = checkOptimality (gap);
	
	displayCut (segs, cut, gap);
	*cut = flow; 

	freeMemory ();
	
	return;

}