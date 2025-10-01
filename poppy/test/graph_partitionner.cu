#include <cstddef> /* NULL */
#include <metis.h>
#include <iostream>
#include <vector>
#include <metis.h>
#include "../src/utils/include/utils.cuh"

int main(){

    Graph g = read_graph("../Graphs/students.csv", ",", false);

std::multimap<int, int> map_out;
    for (auto arc : g.arcs) {
        map_out.insert(arc);
        map_out.insert({arc.second, arc.first});
    }
    // Count edges
    int count = 0;
    for (auto node : g.nodes) {
        count += map_out.count(node);
    }
    idx_t nVertices = g.nodes.size();
    idx_t nEdges    = count / 2;
    idx_t nWeights  = 1;
    idx_t nParts    = 2;

    idx_t objval;
    std::vector<idx_t> part(nVertices, 0);


    // Indexes of starting points in adjacent array
    std::vector<idx_t> xadj;
    // Adjacent vertices in consecutive index order
    std::vector<idx_t> adjncy;
    // Weights of vertices
    // if all weights are equal then can be set to NULL
    std::vector<idx_t> vwgt(nVertices * nWeights, 0);

    // Create the xadj and adjncy arrays
    
    for (auto node : g.nodes) {
        xadj.push_back(adjncy.size());
        for (auto arc : map_out) {
            if (arc.first == node) {
                adjncy.push_back(arc.second);
            }
        }
    }
    xadj.push_back(adjncy.size());
    std::cout << "xadj :" << xadj << std::endl;
    std::cout << "adjncy :" << adjncy << std::endl;

    assert(xadj.size() == nVertices + 1);
    assert(adjncy.size() == 2*nEdges);
    
    int ret = METIS_PartGraphKway(&nVertices,& nWeights, xadj.data(), adjncy.data(),
				       NULL, NULL, NULL, &nParts, NULL,
       				  NULL, NULL, &objval, part.data());

    std::cout << ret << std::endl;
    
    std::cout << "[";
    for(unsigned part_i = 0; part_i < part.size(); part_i++){
	std::cout << part[part_i] << ", ";
    }
    std::cout << "]" << std::endl;
    
    return 0;
}
