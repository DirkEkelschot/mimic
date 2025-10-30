#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>


// Data structure for a node
struct Node {
    int id;
    double x, y, z;
};



// Data structure for an element (simplified)
struct Element {
    int id;
    int type;
    std::vector<int> node_ids;
};



bool read_gmsh(const std::string& filename,
               std::vector<Node>& nodes,
               std::vector<Element>& elements) {
    std::ifstream infile(filename);
    if (!infile.is_open()) return false;

    std::string line;
    enum Section { NONE, NODES, ELEMENTS };
    Section current = NONE;
    int nodes_to_read = 0, elements_to_read = 0;

    while (std::getline(infile, line)) {
        if (line == "$Nodes") {
            current = NODES;
            std::getline(infile, line);
            nodes_to_read = std::stoi(line);
            continue;
        }
        if (line == "$EndNodes") {
            current = NONE;
            continue;
        }
        if (line == "$Elements") {
            current = ELEMENTS;
            std::getline(infile, line);
            elements_to_read = std::stoi(line);
            continue;
        }
        if (line == "$EndElements") {
            current = NONE;
            continue;
        }

        if (current == NODES && nodes_to_read > 0) {
            std::istringstream iss(line);
            Node node;
            iss >> node.id >> node.x >> node.y >> node.z;
            nodes.push_back(node);
            if (--nodes_to_read == 0) current = NONE;
        }

        if (current == ELEMENTS && elements_to_read > 0) {
            std::istringstream iss(line);
            Element elem;
            int num_tags;
            iss >> elem.id >> elem.type >> num_tags;
            // Skip tags
            for (int i=0; i<num_tags; ++i) { int tmp; iss >> tmp; }
            int node_id;
            while (iss >> node_id) {
                elem.node_ids.push_back(node_id);
            }
            elements.push_back(elem);
            if (--elements_to_read == 0) current = NONE;
        }
    }
    infile.close();
    return true;
}





int main(int argc, char** argv)
{

    const char* fm = "mesh/hemisphere.msh";

    std::vector<Element> elems;
    std::vector<Node> nodes;
    read_gmsh(fm,nodes,elems);



}