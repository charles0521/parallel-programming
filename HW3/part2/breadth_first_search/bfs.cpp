#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1
#define ALPHA 12
#define BETA 24

void vertex_set_clear(vertex_set *list)
{
    list->count = 0;
}

void vertex_set_init(vertex_set *list, int count)
{
    list->max_vertices = count;
    list->vertices = (int *)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.

void top_down_step(Graph g, vertex_set *frontier, int *distances, int front_id)
{
    int next_count = 0;

#pragma omp parallel for reduction(+:next_count)
    for (int node = 0; node < g->num_nodes; node++)
    {
        if (frontier->vertices[node] == front_id)
        {

            int start_edge = g->outgoing_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                               ? g->num_edges
                               : g->outgoing_starts[node + 1];

            // attempt to add all neighbors to the new frontier
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
            {
                int outgoing = g->outgoing_edges[neighbor];
                if (distances[outgoing] == NOT_VISITED_MARKER)
                {
                    distances[outgoing] = distances[node] + 1;
                    next_count = next_count + 1;
                    frontier->vertices[outgoing] = front_id + 1;
                }
            }
        }
    }

    frontier->count = next_count;
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution *sol)
{

    vertex_set list1;
    vertex_set_init(&list1, graph->num_nodes);

    int front_id = 1;
    vertex_set *frontier = &list1;

    memset(frontier->vertices, NOT_VISITED_MARKER, sizeof(int) * graph->num_nodes);
    frontier->vertices[ROOT_NODE_ID] = front_id;
    frontier->count++;

    // initialize all nodes to NOT_VISITED
    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0)
    {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(frontier);

        top_down_step(graph, frontier, sol->distances, front_id);
        front_id++;
#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif
    }
}

// void bfs_bottom_step(Graph g, vertex_set *frontier, int *distances, int front_id)
// {
//     int next_count = 0;

// #pragma omp parallel for reduction(+ : next_count)
//     for (int node = 0; node < g->num_nodes; node++)
//     {
//         if (frontier->vertices[node] == NOT_VISITED_MARKER)
//         {
//             int start_edge = g->incoming_starts[node];
//             int end_edge = (node == g->num_nodes - 1) ? g->num_edges : g->incoming_starts[node + 1];

//             for (int neighbor = start_edge; neighbor < end_edge; neighbor++)
//             {
//                 int neighbor_node = g->incoming_edges[neighbor];
//                 if (frontier->vertices[neighbor_node] == front_id)
//                 {
//                     distances[node] = distances[neighbor_node] + 1;
//                     frontier->vertices[node] = front_id + 1;
//                     next_count = next_count + 1;
//                     break;
//                 }
//             }
//         }
//     }
//     frontier->count = next_count;
// }

void bottom_up_step(
    Graph g,
    vertex_set* frontier,
    int *distances,
    int frontier_id)
{
    int next_frontier_cnt = 0;

    #pragma omp parallel
    {
	#pragma omp for reduction(+:next_frontier_cnt)
	for (int i=0; i < g->num_nodes; i++) {
    	    if (frontier->vertices[i] == NOT_VISITED_MARKER){
    	        int start_edge = g->incoming_starts[i];
    	        int end_edge = (i == g->num_nodes-1) ? g->num_edges : g->incoming_starts[i+1];

    	        for (int neighbor = start_edge; neighbor < end_edge; neighbor++) {
		    int neighbor_node = g->incoming_edges[neighbor];

		    if (frontier->vertices[neighbor_node] == frontier_id){
			distances[i] = distances[neighbor_node] + 1;
			frontier->vertices[i] = frontier_id + 1;
    	    	      	next_frontier_cnt++;
    	    	      	break;
		    }
    	        }
    	    }
    	}
    }
    frontier->count += next_frontier_cnt;
}

void bfs_bottom_up(Graph graph, solution *sol)
{
    vertex_set list;
    vertex_set_init(&list, graph->num_nodes);

    vertex_set *frontier = &list;
    int front_id = 1;

    memset(frontier->vertices, NOT_VISITED_MARKER, sizeof(int) * graph->num_nodes);
    frontier->vertices[ROOT_NODE_ID] = front_id;
    frontier->count++;

    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0)
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(frontier);

        // bfs_bottom_step(graph, frontier, sol->distances, front_id);
        bottom_up_step(graph, frontier, sol->distances, front_id);
        front_id++;
#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif
    }
}

void parallel_chooseDirection(bool &currentDirection, int sizeGraph, int sizeFrontier, int sizeNext, int alpha, int beta)
{

    int edgesToCheck;
    double branching_factor = (double)(sizeNext - sizeFrontier) / sizeFrontier;

    // this is the case we the graph is growing
    if (currentDirection && branching_factor > 0)
    {
        edgesToCheck = sizeNext * branching_factor;
        currentDirection = (edgesToCheck < (sizeGraph * branching_factor / alpha));
        // Here the graph is shrinking
    }
    else if (!currentDirection && branching_factor < 0)
    {
        edgesToCheck = sizeFrontier;
        currentDirection = (sizeFrontier < (sizeGraph / beta));
    }
}

void bfs_hybrid(Graph graph, solution *sol)
{
    vertex_set list;
    vertex_set_init(&list, graph->num_nodes);

    vertex_set *frontier = &list;
    int front_id = 1;
    int current_graph_size = graph->num_nodes - 1; // -1 => root
    memset(frontier->vertices, NOT_VISITED_MARKER, sizeof(int) * graph->num_nodes);
    frontier->vertices[ROOT_NODE_ID] = front_id;
    frontier->count++;

    bool directionDown = true; // choose topdown or bottomup

    for (int i = 0; i < graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    sol->distances[ROOT_NODE_ID] = 0;

    int now_frontier_count = frontier->count;

    while (frontier->count != 0)
    {
#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(frontier);
        if (directionDown)
        {
            top_down_step(graph, frontier, sol->distances, front_id);
        }
        else
        {
            // bfs_bottom_step(graph, frontier, sol->distances, front_id);
            bottom_up_step(graph, frontier, sol->distances, front_id);
        }

        front_id++;
        current_graph_size -= frontier->count;
        parallel_chooseDirection(directionDown, current_graph_size, now_frontier_count, frontier->count, ALPHA, BETA);
        now_frontier_count = frontier->count;

#ifdef VERBOSE
        double end_time = CycleTimer::currentSeconds();
        printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif
    }
}
