#ifndef update_h
#define update_h
#include <gsl/gsl_rng.h>

#include "dtype.h"

double ninfection_value();

double nrecover_value();

void remove_vertices(world_line* w);
/** 
 * This function removes vertices from the world-line of the simulation that do not contribute to state changes.
 *
 * Parameters:
 *   w (world_line*): Pointer to the world_line structure that represents the current state of the system.
 *
 * Behavior:
 *   - The function iterates over all vertices in the active sequence (either sequenceA or sequenceB, depending on the flag).
 *   - It checks each vertex to determine if it causes a change in the state of any spins. A vertex is retained if it changes the state
 *     of at least one spin; otherwise, it is removed.
 *   - Counters for the number of infection and recovery events (ninfection and nrecover) are updated based on the types of vertices
 *     that are retained.
 *   - The function toggles the active sequence flag at the end, swapping the roles of sequenceA and sequenceB for the next operation.
 *
 * Outputs:
 *   - The function modifies the world_line structure in-place, reducing the number of vertices and potentially altering which sequence
 *     is active. It also updates global counters for the number of infections and recoveries observed.
 */


void remove_only_fixed_vertices(world_line* w);
/**
 * This function removes fixed (unchanging) vertices from the world-line of a Monte Carlo simulation.
 *
 * Parameters:
 *   w (world_line*): Pointer to the world_line structure representing the current state and configuration of the simulation.
 *
 * Behavior:
 *   - The function iterates over all vertices in the currently active sequence (either sequenceA or sequenceB, depending on the flag).
 *   - It examines each vertex to determine if any changes occur across its associated legs or if any of its legs belong to a
 *     cluster with zero weight, indicating potential for state change.
 *   - Vertices that exhibit any change in state or are part of a dynamic cluster (non-zero weight) are copied to the other sequence for retention.
 *   - The count of infection-related and recovery-related vertices is updated based on the type of interaction they represent.
 *   - This process reduces the number of vertices in the sequence, potentially enhancing performance by focusing computational efforts on dynamic parts of the system.
 *
 * Outputs:
 *   - Modifies the world_line structure in-place, reducing the number of vertices and toggling the active sequence flag to switch between sequences.
 *   - Updates global counters for the number of infections and recoveries observed during the process.
 */


void swapping_graphs(world_line* w, model* m, gsl_rng* rng);
/**
 * This function performs random swaps of vertex types within a world-line according to specific rules based on the type of bond.
 * It's used in the simulation to introduce randomness and to explore different configurations in the phase space of the model.
 *
 * Parameters:
 *   w (world_line*): Pointer to the world_line structure representing the simulation's current state.
 *   m (model*): Pointer to the model structure containing information about the simulation's sites and bonds.
 *   rng (gsl_rng*): Pointer to a GSL random number generator, used to introduce randomness into the graph swapping process.
 *
 * Behavior:
 *   - The function iterates through each vertex in the active sequence (sequenceA or sequenceB, depending on the flag).
 *   - For vertices associated with bonds of types 1, 3, and 5 or types 2, 4, and 6, it randomly assigns a new bond of the same general type
 *     (odd or even) based on a uniformly distributed random value.
 *   - For vertices associated with bond types 7 or 8, it randomly adjusts the bond by adding or subtracting the number of nodes,
 *     effectively swapping between these two types.
 *   - These swaps are designed to maintain the overall connectivity and type balance of the graph while exploring new configurations.
 *
 * Outputs:
 *   - The function modifies the bonds of vertices in the active world-line sequence directly, altering the graph structure used in
 *     subsequent simulation steps.
 */


void insert_vertices(world_line* w, model* m, gsl_rng* rng);
/**
 * This function inserts new vertices into the world-line of the simulation based on a sampling of a uniform sequence,
 * which is influenced by the model's site weight and the simulation's inverse temperature.
 *
 * Parameters:
 *   w (world_line*): Pointer to the world_line structure representing the simulation's current state.
 *   m (model*): Pointer to the model structure containing information about the system's sites and bonds.
 *   rng (gsl_rng*): Pointer to a GSL random number generator, used for generating the random positions of new vertices.
 *
 * Behavior:
 *   - The function first adjusts the length of the world-line sequence to accommodate the insertion of new vertices.
 *   - It copies the current state of the system into a temporary state array.
 *   - The function iterates over the sampled sequence of insertion times and inserts new vertices at these times
 *     if they satisfy certain conditions defined by the model's rules for insertion.
 *   - Each new vertex inserted is initialized with the appropriate bond and state information, and placed in the sequence
 *     in chronological order.
 *   - The original vertices are also copied into the new sequence, maintaining their original order.
 *   - At the end of the function, the sequence with the newly inserted vertices becomes the active sequence.
 *
 * Outputs:
 *   - This function modifies the world_line structure in-place by reallocating its vertex sequence and updating
 *     the count and arrangement of vertices. It also toggles the flag that determines which of two possible
 *     vertex sequences is active.
 */


void clustering(world_line* w, model* m);
/**
 * This function implements the clustering algorithm for the world-line Monte Carlo simulation,
 * linking vertices based on their interactions to form clusters.
 *
 * Parameters:
 *   w (world_line*): Pointer to the world_line structure representing the current state of the simulation.
 *   m (model*): Pointer to the model structure containing details about the bonds and rules for linking vertices.
 *
 * Behavior:
 *   - The function initializes the first and last indices for each site to track the start and end of clusters.
 *   - It iterates through each vertex in the active sequence (sequenceA or sequenceB, depending on the flag).
 *   - For each vertex, it applies the linking rules defined in the model based on the type of bond associated with the vertex.
 *   - These rules determine how vertices are connected within the cluster framework, setting up the foundation for collective
 *     updates during the simulation.
 *   - The function updates the `cluster` and `weight` arrays in the world_line structure to reflect the connections and weights
 *     between vertices as dictated by the rules.
 *   - Optionally (as noted by commented code), it can handle open boundary conditions by linking the first and last vertices
 *     in each cluster, though this is disabled by default in the provided code snippet.
 *
 * Outputs:
 *   - The function modifies the world_line structure in-place by setting up links between vertices based on the model's rules.
 *     These links are used in later steps of the simulation to perform updates across connected vertices simultaneously.
 */

void cluster_statistic(world_line* w, model* m);
/**
 * This function computes statistics related to clusters within the world-line Monte Carlo simulation, tracking the distribution
 * and dynamics of clusters over time. It evaluates the ratio of free clusters to total clusters and their respective sizes,
 * as well as the temporal dynamics of infections.
 *
 * Parameters:
 *   w (world_line*): Pointer to the world_line structure representing the simulation's current state and setup.
 *   m (model*): Pointer to the model structure providing details on bonds and their indices.
 *
 * Behavior:
 *   - Initializes or updates storage arrays for cluster and infection statistics as needed, checking for memory allocation success.
 *   - Resets counters and accumulators for statistical metrics at the start of computation.
 *   - Iterates over all vertices in the active sequence to evaluate and update cluster and infection statistics based on the state
 *     of each vertex and its relationship to others via bonds.
 *   - Tracks changes in cluster and infection states over time, accumulating data on cluster sizes and the duration of infection states.
 *   - Computes ratios of free to total clusters and their sizes to assess the dynamism and spread within the simulation.
 *   - Logs computed statistics to the console and appends detailed records to a file for further analysis.
 *
 * Outputs:
 *   - Modifies global and static variables to store and update statistical data.
 *   - Outputs to the console for immediate observation of the simulation's state and progression.
 *   - Writes detailed cluster and infection statistics to a file named 'cluster_statistic.txt' for persistence and later analysis.
 */


void flip_cluster(world_line* w, gsl_rng* rng);
/**
 * This function performs the flip operation on clusters within a world-line Monte Carlo simulation.
 * It determines whether each cluster will flip its state based on random choices and the cluster's
 * associated weight.
 *
 * Parameters:
 *   w (world_line*): Pointer to the world_line structure representing the current state of the simulation.
 *   rng (gsl_rng*): Pointer to a GSL random number generator used to introduce randomness in the flip decision.
 *
 * Behavior:
 *   - The function iterates through all vertices in the active sequence (sequenceA or sequenceB, depending on the flag).
 *   - For each vertex, it processes each leg, determining the root of its cluster and deciding if the cluster's state will flip.
 *   - The decision to flip is based on the cluster's weight and a random value generated for each cluster.
 *   - If a cluster is determined to flip, all states in the cluster are inverted.
 *   - After processing the vertices, it updates the initial and final states of each site in the simulation based on the active sequence
 *     or random values if no active vertex influences the site.
 *
 * Outputs:
 *   - The function modifies the state arrays within the world-line structure directly, affecting the simulation's subsequent behavior.
 *   - It also updates the initial and projected state arrays (`istate` and `pstate`) for each site, ensuring that the simulation
 *     reflects the changes made during this operation.
 */


//int check_periodic(world_line* w, model* m);

void snapshot_show(world_line* w, model* m, FILE* file);

#endif
