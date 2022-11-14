#include "page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double *solution, double damping, double convergence)
{

  // initialize vertex weights to uniform probability. Double
  // precision scores are used to avoid underflow for large graphs

  int numNodes = num_nodes(g);
  double equal_prob = 1.0 / numNodes;
  for (int i = 0; i < numNodes; ++i)
  {
    solution[i] = equal_prob;
  }

  double *score_new = (double *)malloc(sizeof(double) * numNodes);
  double *score_old = (double *)malloc(sizeof(double) * numNodes);
  double sum_no_outgoing;

  //init
  for(int i=0; i<numNodes; ++i)
  {
    score_old[i] = equal_prob;
  }


  // #pragma omp for
  double global_diff;

  do
  {
    #pragma omp parallel for
    for(int i=0; i<numNodes; ++i)
    {
      score_new[i] = 0.0;
      for(int j=0; j<incoming_size(g, i); ++j)
      {
        int intcoming_point = g->incoming_edges[g->incoming_starts[i] + j];
        score_new[i] += score_old[intcoming_point]/outgoing_size(g, intcoming_point);
      }
      score_new[i] = (damping * score_new[i]) + (1.0 - damping) / numNodes;
    }

    sum_no_outgoing = 0.0;
    // #pragma omp parallel for
    for(int i=0; i<numNodes; ++i)
    {
      if(outgoing_size(g, i) == 0)
        sum_no_outgoing = sum_no_outgoing + (damping * score_old[i] / numNodes);
    }
    #pragma omp parallel for
    for(int i=0; i<numNodes; ++i)
    {
      score_new[i] += sum_no_outgoing;
    }

    global_diff = 0.0;
    for(int i=0; i<numNodes; ++i)
    {
      global_diff += abs(score_new[i] - score_old[i]);
      score_old[i] = score_new[i];
    }

  }while(!(global_diff < convergence));
  #pragma omp parallel for
  for(int i=0; i<numNodes; ++i)
  {
    solution[i] = score_new[i];
  }

  free(score_new);
  free(score_old);
  /*
     For PP students: Implement the page rank algorithm here.  You
     are expected to parallelize the algorithm using openMP.  Your
     solution may need to allocate (and free) temporary arrays.

     Basic page rank pseudocode is provided below to get you started:

     // initialization: see example code above
     score_old[vi] = 1/numNodes;

     while (!converged) {

       // compute score_new[vi] for all nodes vi:
       score_new[vi] = sum over all nodes vj reachable from incoming edges
                          { score_old[vj] / number of edges leaving vj  }
       score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

       score_new[vi] += sum over all nodes v in graph with no outgoing edges
                          { damping * score_old[v] / numNodes }

       // compute how much per-node scores have changed
       // quit once algorithm has converged

       global_diff = sum over all nodes vi { abs(score_new[vi] - score_old[vi]) };
       converged = (global_diff < convergence)
     }

   */
}
