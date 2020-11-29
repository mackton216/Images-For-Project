# EE599 Final Project Report
## Travelling Trojan Problem
### Ashwin Sivakumar, Mackton Vishal Mendonca

## Overview of Design

**insert image here**

- The travelling salesman problem is a well known [NP-hard problem](https://en.wikipedia.org/wiki/NP-hardness) that is used in computer science for combinatorial optimizations.
- Our design employs the use of structures and helper functions to try and solve this problem.
- We have made the use of the Haversine formula to find the distance in miles between various points differing in longitude and latitudes.
- We use the Dijkstra algorithm along with a priority queue to find the Shortest Path between 2 points on the map.
- We use Brute Force and several helper functions to reach our goal in helping our Travelling Trojan reach all the nodes in a given map and return in the shortest path possible.
- We have also used two heuristic methods to find “good enough outputs” for the Travelling Salesman Problem namely:
    - 2 - Opt
    - Genetic Algorithm
## The Node Class
```cpp
class Node {
  public:
    Node(){};
    Node(const Node &n){id = n.id; lat = n.lat; lon = n.lon; name = n.name; neighbors = n.neighbors;};
    std::string id;    // A unique id assign to each point
    double lat;        // Latitude
    double lon;        // Longitude
    std::string name;  // Name of the location. E.g. "Bank of America".
    std::vector<std::string>neighbors;  // List of the ids of all neighbor points.
};
```
- Contains a parameterized constructor with datatype Node
- The constructor defines objects for Node IDs, latitudes, longitudes, names of places represented by the nodes and neighbouring nodes

## Primary Code Functions:
```cpp
// Returns a vector of names given a partial name.
std::vector<std::string> Autocomplete(std::string name);
// Given the name of two locations, it should return the **ids** of the nodes on the shortest path using Dijkstra Priority Queue
std::vector<std::string> CalculateShortestPath(std::string location1_name, std::string location2_name);
// Travelling Salesman Using Brute Force
std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan(std::vector<std::string> &location_ids);
// Travelling Salesman Using 2-Opt
std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan_2opt(std::vector<std::string> &location_ids);
// Travelling Salesman Using Genetic Algorithm
std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan_GeneticAlgo(std::vector<std::string> &location_ids);
```

## Generating the weight matrix
**insert image**
- Objective was to create an adjacency matrix while mapping every location id to an index.
- This way the core functions needed minimal alteration.
- Made implementation more modular.

## Shortest Path - Dijkstra Priority Queue
**insert image**
- Idea is to traverse through graph and find next minimum distanced node and update distances of all children.
- Store parent nodes for each node in every iteration.
- Leverages Priority queue functionality to efficiently obtain next minimum unvisited node.
- Parent nodes are accessed and recursed to obtain the path in addition to the shortest distance.

## Travelling Salesman – Brute Force
**insert image**
- Aimed at solving the TSP using a greedy technique. 
- Method always returns an optimal result but computationally expensive.
- Intuitive technique inspired from DFS search + Permutations.

## Travelling Salesman – 2 opt
**insert image**
- Aimed at solving the TSP using an approximate technique.
- Utilizes 2-edge swapping
- Method does not always return optimal distance / path , but it is expected to give a “very good” approximate of the ideal output.
- Computationally less expensive in comparison with Brute Force.

## Travelling Salesman – Genetic Algorithm
**insert image**


---

## Every function explained:
### GetLat, GetLon, GetName:
```cpp
double TrojanMap::GetLat(std::string id) 
{ 
 return data[id].lat;
}

double TrojanMap::GetLon(std::string id) 
{  
   return data[id].lon;
}

std::string TrojanMap::GetName(std::string id) 
{ 
  return data[id].name;
}
```
- Using the attributes of the variable ‘data’ defined in the cpp.h file.
- Data is a map consisting of strings and Nodes defined in the .csv file.
- String “id” is user input and is mapped to the latitude, longitude and name (if any) of the respective Node id.
- Each of these have a time complexity of O(1).

### GetNeighborIDs
```cpp
std::vector<std::string> TrojanMap::GetNeighborIDs(std::string id) 
{   

    std::vector<std::string> result;

    for(int i = 0; i < data[id].neighbors.size(); i++)
    {
      result.push_back(data[id].neighbors[i]);
    }
  return result;
} 
```
- Uses the neighbors attribute of the map ‘data’ in the cpp.h file to retrieve the neighboring IDs of each Node ID as defined the .CSV file.
- Returns a vector of strings.
- O(N) where N is the number of nodes

### CalculateDistance
```cpp
  // TODO: Use Haversine Formula:
  // dlon = lon2 - lon1;
  // dlat = lat2 - lat1;
  // a = (sin(dlat / 2)) ^ 2 + cos(lat1) * cos(lat2) * (sin(dlon / 2)) ^ 2;
  // c = 2 * arcsin(min(1, sqrt(a)));
  // distances = 3961 * c;

  // where 3961 is the approximate radius of the earth at the latitude of
  // Washington, D.C., in miles
double TrojanMap::CalculateDistance(const Node &a, const Node &b) 
{

  double dlon, dlat, a1, c, distance;

  dlon = b.lon - a.lon;
  dlon = dlon * (M_PI/180);
  dlat = b.lat - a.lat;
  dlat = dlat * (M_PI/180);

  a1 = (pow((std::sin(dlat / 2)),2)) + std::cos(a.lat* (M_PI/180)) * std::cos(b.lat* (M_PI/180)) * (pow((std::sin(dlon / 2)),2));
  c = 2 * std::asin(MIN(1, sqrt(a1)));

  distance = (3961 * c);

  return distance;
}
```
- Using the Haversine formula to calculate distance between 2 Nodes which contain respective latitudes and longitudes. 
- Returns a double.
- O(1)

### CalculatePathLength
```cpp
double TrojanMap::CalculatePathLength(const std::vector<std::string> &path) 
{  

    int i = 0, j = 1;
    double sum = 0;
    while(j<path.size())
    {
      sum += CalculateDistance(data[path[i]], data[path[j]]);
      i++;
      j++; 
    }
   
  return sum;
}
```
- Uses CalculateDistance function to find distance between the latitudes and longitudes of neighboring nodes that are in the path.
- Sums it up to find total length of a path taken.
- Returns a double.
- O(path.size() - 1).

### GetNode and GetNodeFromName
```cpp
Node TrojanMap::GetNode(std::string id)
{
  std::map<std::string, Node>::iterator it;

  for(it = data.begin(); it!=data.end();it++)
     {
       if(it->first.compare(id) == 0)
       {
         return it->second;
         break;         
       }
     }
}


Node TrojanMap::GetNodeFromName(std::string locname)
{
  std::map<std::string, Node>::iterator it;
  int flag = 0;
  for(it = data.begin(); it!=data.end();it++)
     {
       if(it->second.name.compare(locname) == 0)
       {
         flag = 1;
         return it->second;
         break;         
       }
     }
}
```
- The id of the node is compared to the ids present in the .CSV file and the subsequent node is returned.
- O()
- The location name is searched for in the .CSV file and then corresponding name is returned. 
- O()








