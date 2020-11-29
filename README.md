# EE599 Final Project Report
## Travelling Trojan Problem
### Ashwin Sivakumar, Mackton Vishal Mendonca

## Overview of Design

<p align="left"><img src="https://user-images.githubusercontent.com/47607653/100549270-d6462780-3226-11eb-9fbb-97508216af2c.png" alt="TSP" width="500"/></p>

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
- **Autocomplete** - To take a part of a location name from the user and return all places that contain that part of the name.
- **GetPosition** - To find the latitude and longitude of a given location
- **CalculateShortestPath** - Calculate the shortest path between 2 locations
- **TravellingTrojan** -  Given a set of locations, find the shortest route between after visiting every location once and returning using Brute Force.
- **Heuristical solutions to TravellingTrojan using 2_opt and GeneticAlgo.**

## Updated Menu
<p align="center"><img src="https://user-images.githubusercontent.com/47607653/100550142-b285e000-322c-11eb-9e02-90bee1b77302.PNG" alt="NewMenu" width="500"/></p>

- Our updated menu allowing 7 options.

## Generating the weight matrix
<p align="center"><img src="https://user-images.githubusercontent.com/47607653/100549604-c92a3800-3228-11eb-8bb8-d6264d32434b.PNG" alt="GenWeightMatrix" width="500"/></p>

- Objective was to create an adjacency matrix while mapping every location id to an index.
- This way the core functions needed minimal alteration.
- Made implementation more modular.

## Shortest Path - Dijkstra Priority Queue
<p align="center"><img src="https://user-images.githubusercontent.com/47607653/100549661-15757800-3229-11eb-9825-2d8727353a41.PNG" alt="SP_Djikstra" width="500"/></p>

- Idea is to traverse through graph and find next minimum distanced node and update distances of all children.
- Store parent nodes for each node in every iteration.
- Leverages Priority queue functionality to efficiently obtain next minimum unvisited node.
- Parent nodes are accessed and recursed to obtain the path in addition to the shortest distance.

## Travelling Salesman – Brute Force
<p align="center"><img src="https://user-images.githubusercontent.com/47607653/100549672-2aeaa200-3229-11eb-80af-e367c020347d.PNG" alt="TSP_brute" width="500"/></p>

- Aimed at solving the TSP using a greedy technique. 
- Method always returns an optimal result but computationally expensive.
- Intuitive technique inspired from DFS search + Permutations.

## Travelling Salesman – 2 opt
<p align="center"><img src="https://user-images.githubusercontent.com/47607653/100549679-40f86280-3229-11eb-94ab-d3c72d6d8ab5.PNG" alt="TSP_2Opt" width="500"/></p>

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

### Autocomplete
```cpp
std::vector<std::string> TrojanMap::Autocomplete(std::string name) {

  int n = name.size();
  std::string temp;
  std::vector<std::string> results = {};
  
  std::transform(name.begin(), name.end(), name.begin(),
    [](unsigned char c){ return std::tolower(c); });

  std::map<std::string, Node>::iterator it;

  for(it = data.begin(); it!=data.end();it++)
  { 

      if(it->second.name.empty() == 0)
      {
           
        temp.append(it->second.name.begin(),it->second.name.begin() + n);

        std::transform(temp.begin(), temp.end(), temp.begin(),
        [](unsigned char c1){ return std::tolower(c1); });

        if(name.compare(temp)==0)
          {
            results.push_back(it->second.name);
          }
      }
      temp.clear();
  }
  return results;
}
```
- An input from the user is stored as a string and converted into lower case.
- A map is created with the string as the key and the Node as the data.
- A temporary string is created to hold a part of the name. This temp string is used to scan the .CSV file and return strings which match the temp string.
- If they match the temp string, the string is pushed into the results vector and returned to the user.
- O()

### GetPosition
```cpp
std::pair<double, double> TrojanMap::GetPosition(std::string name) {
  std::pair<double, double> results(-1, -1);
  std::map<std::string, Node>::iterator it;

  for(it = data.begin(); it!=data.end();it++)
     {
       if(it->second.name.compare(name) == 0)
       {
         results.first = it->second.lat;
         results.second = it->second.lon;
         return results;
         break;         
       }
     }
}
```
- An iterator is used to traverse the .CSV file holding the nodes.
- This iterator compares the name string that the user inputs with the name string in the file.
- If there is a match in the name, the corresponding latitudes and longitudes of the name are retrieved from the .CSV file, stored in a pair and displayed to the user.
- O(data.size -1)

### DikjstraPriorityQueue
```cpp
std::vector<double> TrojanMap::DijkstraPriorityQueue(int source, std::vector<std::vector<double>> weight_, std::map<int, int> &prev)
{
  std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>>q;
  std::vector<double> d(weight_.size(), DBL_MAX), d_temp(weight_.size(), DBL_MAX);

  d[source] = 0;
  q.push(std::make_pair(0, source));

  while(!q.empty())
  { 

    int u = q.top().second;
    q.pop();
    for(int j = 0; j < weight_.size(); j++)
    {
      if(d[j] > d[u] + weight_[u][j])
      {
        d[j] = d[u] + weight_[u][j];
        q.push(std::make_pair(d[j], j));
        prev[j] = u;
      }
    }
  }
  return d;
}
```
- Uses the Dijkstra algorithm to find the shortest path between all the nodes in the map.
- The priority queue is used to make it more efficient to find the shortest path.
- The function starts at the source node.
- It travels to next node which is the closer than its other neighbouring nodes and updates the distance between the source node and next node.
- It then moves on to the neighbouring node which is the shorter distance from the previous node and updates that distance.
- Ultimately returns the shortest path taken to cover all the nodes.
- We use this function to find the shortest path taken between 2 nodes defined by the user. 
- O(n^2)

### CalculateShortestPath

```cpp
std::vector<std::string> TrojanMap::CalculateShortestPath(
  
  std::string location1_name, std::string location2_name) {
  std::vector<std::string> path = {};
  Node a, b; 
  double dist;
  int m = 0;
  std::map<std::string, Node>::iterator it;

  int flag = 0;
  for(it = data.begin(); it!=data.end();it++)
  {
    if(it->second.name.compare(location2_name) == 0)
    {
      flag = 1;
      break;         
    }
  }

  if(flag == 0)
  {
    path = {};
    return path;
  }


  a = GetNodeFromName(location1_name);
  b = GetNodeFromName(location2_name); 
  
  std::vector<std::vector<double>> weight_(2237, std::vector<double>(2237, DBL_MAX));

  std::map<std::string, Node>::iterator itr1, itr2;

  std::map<std::string, int>index;
  std::map<std::string, int>::iterator ie;
  std::map<int, std::string> rev_in;

  

  std::map<int, int> prev;

  for(itr1 = data.begin(); itr1!=data.end(); itr1++) // current node
  {
    index[itr1->first] = m;
    m++;
  }

   for(ie = index.begin(); ie != index.end(); ie++)
  {
      rev_in[ie->second] = ie->first;
  }

  int source = index[a.id];
  int destination = index[b.id];

  if(source == destination)
  {
    path.push_back(rev_in[source]);
    return path;
  }

 
  int found = 0, mx = 0;

  for(itr1 = data.begin(); itr1!=data.end(); itr1++) // current node
  { 

    for(itr2= data.begin(); itr2!=data.end(); itr2++) // iterate through all the nodes
    {
      for(int k = 0; k < itr1->second.neighbors.size(); k++) // check neighbors of current node
      {
        if(found == 1)
        {
          continue;
        }
        if(itr1->second.neighbors[k].compare(itr2->first) == 0) // check if neighborID and node ID match, if so find distance and update AM
        {
          found = 1;
          dist = CalculateDistance(itr1->second, itr2->second);
          weight_[index[itr1->first]][index[itr2->first]] = dist;
        }
      }
      if(found == 0)
      {
        weight_[index[itr1->first]][index[itr2->first]] = DBL_MAX;
      }
      found = 0;
    }

  }
  
   for(itr1 = data.begin(); itr1!=data.end(); itr1++) // diagonal elements set to zero
  {
    weight_[index[itr1->first]][index[itr1->first]] = 0;
  }

  for(int i = 0; i < 2237; i++)
  {
      if(weight_[source][i]!= DBL_MAX)
      {
        prev[i] = source;
      }
      else
      {
        prev[i] = DBL_MAX;
      }
  }
 
  std::vector<double> x = DijkstraPriorityQueue(source, weight_, prev);  

  int mi = destination;

  path.push_back(rev_in[destination]);
  while(mi!=source)
  { 
    path.push_back(rev_in[prev[mi]]);
    mi = prev[mi];
  }

  std::reverse(path.begin(),path.end());
  return path;
  
}
```
- This function takes in 2 strings from the user. These strings represent the names of the places on the map.
- The GetNodeFromName function then retrieves the nodes with their respective latitudes and longitudes from the .CSV file and they are assigned to 2 integer values.
- Starting from the source node, the distances of the neighbouring nodes are calculated. 
- The neighbouring node that is the closest to the source is found and the distance between them is stored. This is repeated for all the nodes connecting the source to the destination node. 
- A weight matrix is created from evaluation of distances of neighbouring nodes by pushing back these distances into a vector of double vectors.
- The reverse of the path generated from the source to destination is evaluated to verify the path.
- These values are passed as parameters for the DijkstraPriorityQueue function.
- The path is returned.

### TSPHELPER: For Brute Force Method

```cpp
double TrojanMap::TSPHELPER(int start, int cur_node, double cur_cost, std::vector<std::string> path, std::vector<std::vector<double>> weight_, std::map<int, std::string>& rev_index, std::pair<double, std::vector<std::vector<std::string>>>&progress)
{
 
  double result = INT16_MAX;

    path.push_back(rev_index[cur_node]);

    if(path.size() == weight_.size())
    {
      path.push_back(rev_index[start]);
      progress.second.push_back(path);
      return cur_cost + weight_[cur_node][start];
    }

    for(int i = 0; i < weight_.size(); i++)
    {
      if(i!=cur_node && std::find(path.begin(), path.end(), rev_index[i]) == path.end())
      { 
        result = std::min(result, TSPHELPER(start, i, cur_cost + weight_[cur_node][i], path, weight_, rev_index, progress));
      }
    }

  progress.first = result;
  return result;
}
```
- The TSPHELPER function is a function written to evaluate the various paths that the Trojan can take to reach a destination and come back to the source node.
- The parameters passed here are the start int, the current node, the “cost” of travelling through the nodes, the path taken, the a vector of vectors of weights, a map of integers and strings containing the indices in reverse and a pair of double and vector of vector of strings containing the total progress.
- Each node that the Trojan visits, the distance between the nodes are evaluated and stored.
- Each path that the Trojan visits is stored in a vector.
- This happens recursively for each of the neighbouring nodes from the source to the destination.
- Once all the paths have been traversed, the path with the least distance is chosen.

### TravellingTrojan

```cpp
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan(std::vector<std::string> &location_ids) 
{
 
  std::pair<double, std::vector<std::vector<std::string>>> progress;
  std::map<std::string, int>index;
  std::map<int, std::string>rev_index;

  
  double dist;
  std::vector<std::vector<double>> weight_(location_ids.size(), std::vector<double> (location_ids.size()));
  // create indices

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    index[location_ids[i]] = i;
  }

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    rev_index[i] = location_ids[i];
  }
  
  // generate weight matrices
  for(int i = 0; i < location_ids.size(); i++)
  { 
    for(int j = 0; j < location_ids.size(); j++)
    {
      if(i == j)
      {
        weight_[i][j] = 0;
      }

      else
      {      
      dist = CalculateDistance(data[location_ids[i]], data[location_ids[j]]);      
      weight_[i][j] = dist;
      }
    }
  }
  
  std::vector<std::string> path;

  progress.first = 0;
  int start = 0;
  double d;
  double rs = TSPHELPER(start, start, 0, path, weight_, rev_index, progress);
  progress.first = rs;

  for(int i = 0; i < progress.second.size(); i++)
  {
    d = CalculatePathLength(progress.second[i]);
    if(d == progress.first)
    {
    progress.second.push_back(progress.second[i]);
    progress.second.erase(progress.second.begin()+i);
    }
  }

  
  return progress;
}
```
- The parameters passed are a vector of strings containing the location ids.
- The function is used to create a weight matrix. This weight matrix contains the various distances covered by the Trojan while traversing through the map
- These sizes are passed to the CalculateDistance function to find the total distance the Trojan travels after passing through each node between source and destination once. 
- The progress i.e. node ids in every path chosen are stored. 
- The weight matrix along with the first node, final node, paths taken, the reverse path and the progress of each path is passed to the TSPHELPER function to pick the shortest path.
- The CalculatePathLength function is called to find the total distance of the shortest path in miles.
- The various nodes covered in path determined as the shortest are returned.
- O(n!)

### NN
```cpp
std::vector<int> TrojanMap::NN(std::vector<std::vector<double>> weight_, double &dist)
{
  std::vector<int> visited;
  std::vector<int> path;
  int cur = 0, j, start = 0;

  path.push_back(cur);
  visited.push_back(cur);
  
  
  int min_id;
  double min;  

 
  while(visited.size() < weight_.size())
  {
    min = DBL_MAX;
    for(j = 0; j < weight_.size(); j++)
    {
      if(std::find(visited.begin(), visited.end(), j) == visited.end())
      {
        if(weight_[cur][j] < min)
        {
          min = weight_[cur][j];
          min_id = j;
        }
      }
    }
    dist += weight_[cur][min_id];
    cur = min_id;
    path.push_back(min_id);
    visited.push_back(min_id);
  }
  dist+=weight_[cur][start];

  path.push_back(start);

  return path;
}
```


### Distance
```cpp
double TrojanMap::distance(std::vector<int> pt, int i, std::vector<std::vector<double>> weight_)
{
  double dist = 0;
  int start = i;

  for(int j = 0; j < weight_.size() - 1;j++)
  {
    dist += weight_[pt[i]][pt[i+1]];
    i++;
  }

  dist+=weight_[pt[i]][pt[start]];

 
  return dist;
}
```

### TWO_OPT

```cpp
std::pair<double, std::vector<std::vector<int>>> TrojanMap::TWO_OPT(std::vector<std::vector<double>> weight_, std::pair<double, std::vector<std::vector<int>>> result)
{
  
  std::vector<int> shortest = result.second[0];
  std::vector<int> s_temp;
  std::vector<int> mp, m_temp;
  std::vector<int> temp_store;
  int j = 0, p;
  double dist = result.first;

  int flag = 1,x=0;
  while(flag == 1)
  {
    mp = result.second[result.second.size()-1];
    mp.pop_back();
    mp.insert(std::end(mp), std::begin(mp), std::end(mp));
    s_temp = shortest;
    flag = 0;
    for(j = 0; j < weight_.size(); j++)
    { 
      m_temp = mp;
      std::reverse(m_temp.begin() + j + 1, m_temp.begin() + j + 4);

      // calculate new distance
      dist = distance(m_temp, j, weight_);
      if(dist < result.first)
      {
        flag = 1;
        
        result.first = dist;
        p = j;
        for(int k = 0; k < weight_.size(); k++)
        {
          temp_store.push_back(m_temp[p]);
          p++;
        }
        temp_store.push_back(m_temp[j]);
        result.second.push_back(temp_store);
        shortest = temp_store;        
        temp_store.clear();
        break;
      }
    }
  }

return result;
}
```

### TravellingTrojan_2opt
```cpp
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_2opt(std::vector<std::string> &location_ids)
{
  std::pair<double, std::vector<std::vector<std::string>>> progress;
  std::map<std::string, int>index;
  std::map<int, std::string>rev_index;
  

  std::vector<std::string> loc_temp = location_ids;
  loc_temp.push_back(location_ids[0]);

  // double d_check = CalculatePathLength(loc_temp);
  // std::cout<<d_check<<"\n";
  

  double d = 0;  

  std::vector<std::vector<double>> weight_(location_ids.size(), std::vector<double> (location_ids.size()));

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    index[location_ids[i]] = i;
  }


  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    rev_index[i] = location_ids[i];
  }
  
  // generate weight matrices
  for(int i = 0; i < location_ids.size(); i++)
  { 
    for(int j = 0; j < location_ids.size(); j++)
    {
      if(i == j)
      {
        weight_[i][j] = 0;
      }

      else
      {      
      d = CalculateDistance(data[location_ids[i]], data[location_ids[j]]);      
      weight_[i][j] = d;
      }
    }
  }

  // Start 2-opt


  double path_dist = 0;
  
  //Obtain nearest neighbor heuristic before starting 2 opt

  auto nn = NN(weight_, path_dist);
  std::cout<<"\nDistance before 2 opt"<<" "<<path_dist<<"miles\n";
  std::pair<double, std::vector<std::vector<int>>>result;
 

  std::vector<int> path_temp;

  result.first = path_dist;
  result.second.push_back(nn);
  

  // Call 2opt 

  result = TWO_OPT(weight_, result);


  int z = result.second.size();
  int o = weight_.size()+1;
  std::vector<std::vector<std::string>> path_history;

  std::vector<std::string> temp_path;
  std::vector<int> temp_ind;

  for(int u = 0; u < result.second.size(); u++)
  {
    temp_ind = result.second[u];
    
    for(int l = 0; l < temp_ind.size(); l++)
    {
      temp_path.push_back(rev_index[temp_ind[l]]);
    }
    path_history.push_back(temp_path);
    temp_path.clear();
    temp_ind.clear();
  }


  progress.first = result.first;

  std::cout<<"Distance after 2 opt"<<" "<<progress.first<<" miles\n";
  std::cout<<"Improvement of"<<" "<<path_dist - progress.first<<" miles\n\n";
  progress.second = path_history;

  return progress;
}
```

### TSPUtil
```cpp
void TrojanMap::TSPUtil(std::vector<std::vector<double>> map) 
{ 
    // Generation Number 
    int gen = 1; 
    // Number of Gene Iterations 
    int gen_thres = 5; 
    // Population Size
    int POP_SIZE = 10;
  
    std::vector<struct individual> population; 
    struct individual temp; 
  
    // Populating the GNOME pool. 
    for (int i = 0; i < POP_SIZE; i++) { 
        temp.gnome = create_gnome(map); 
        temp.fitness = cal_fitness(temp.gnome, map); 
        population.push_back(temp); 
    } 
  
    std::cout << "\nInitial population: " << std::endl 
         << "GNOME     FITNESS VALUE\n"; 
    for (int i = 0; i < POP_SIZE; i++) 
        std::cout << population[i].gnome << "   "
             << population[i].fitness << std::endl; 
    std::cout << "\n"; 
  
    bool found = false; 
    int temperature = 10000; 
  
    // Iteration to perform 
    // population crossing and gene mutation. 
    while (temperature > 1000 && gen <= gen_thres) { 
        std ::sort(population.begin() , population.end() ,[](individual ptr_l , individual ptr_r) { return ptr_l < ptr_r;} );
        std::cout << "\nCurrent temp: " << temperature << "\n"; 
        std::vector<struct individual> new_population; 
  
        for (int i = 0; i < POP_SIZE; i++) { 
            struct individual p1 = population[i]; 
          
            while (true) { 
               
                std::string new_g = mutatedGene(p1.gnome, map); 
                struct individual new_gnome; 
                new_gnome.gnome = new_g; 
                new_gnome.fitness = cal_fitness(new_gnome.gnome, map); 
                
                if (new_gnome.fitness <= population[i].fitness) { 
                    new_population.push_back(new_gnome); 
                    break; 
                } 
                else { 
                      // Accepting the rejected children at 
                    // a possible probablity above threshold. 
                    float prob = pow(2.7, 
                                     -1 * ((float)(new_gnome.fitness 
                                                   - population[i].fitness) 
                                           / temperature)); 
                    if (prob > 0.5) { 
                        new_population.push_back(new_gnome); 
                        break; 
                    } 
                } 
            } 
        } 
  
        temperature = cooldown(temperature); 
        population = new_population; 
        std::cout << "Generation " << gen << " \n"; 
        std::cout << "GNOME     FITNESS VALUE\n"; 
  
        for (int i = 0; i < POP_SIZE; i++) 
            std::cout << population[i].gnome << "   "
                 << population[i].fitness << std::endl; 
        gen++; 
    } 
}
```

### TravellingTrojan_GeneticAlgo
```cpp
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_GeneticAlgo(std::vector<std::string> &location_ids) 

{
 
  std::pair<double, std::vector<std::vector<std::string>>> progress;
  std::map<std::string, int>index;
  std::map<int, std::string>rev_index;

  
  double dist;
  std::vector<std::vector<double>> weight_(location_ids.size(), std::vector<double> (location_ids.size()));
  // create indices

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    index[location_ids[i]] = i;
  }

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    rev_index[i] = location_ids[i];
  }
  
  // generate weight matrices
  for(int i = 0; i < location_ids.size(); i++)
  { 
    for(int j = 0; j < location_ids.size(); j++)
    {
      if(i == j)
      {
        weight_[i][j] = 0;
      }

      else
      {      
      dist = CalculateDistance(data[location_ids[i]], data[location_ids[j]]);      
      weight_[i][j] = dist;
      }
    }
  }

  TSPUtil(weight_);
  return progress;
}
```
---
### Results
#### Shortest Path Highlighted
##### Target to Popeyes Louisiana Kitchen - A distance of 1.53852 miles
<p align="center"><img src="https://user-images.githubusercontent.com/47607653/100550161-dcd79d80-322c-11eb-9bfa-84012dbac1d6.PNG" alt="TargetToPopeyes" width="500"/></p>

##### Ralphs to ChickFil-A - A distance of 0.74044 miles
<p align="center"><img src="https://user-images.githubusercontent.com/47607653/100550163-e234e800-322c-11eb-8b47-0180cb299dee.PNG" alt="RalphsToChickFilA" width="500"/></p>

---
### Testing the Algorithms
- The algorithms were tested using the tests provided by the TAs for validation of the solution.
- Designed tests ourselves to verify our code’s error handling. This was done in 2 ways.
    - Modifying the test cases generated by the TAs.
    - Writing our own test cases to better understand our code.
- Designed test cases to handle corner cases so that outlying inputs/conditions would not affect our final result.
    - For example, tested for Node ID not available in the .CSV file to alert the user.
-

---

### Learnings - _Ashwin Sivakumar_
- **Intuition behind shortest path algorithms.** Knowledge gained will help take an efficient and informed approach before tackling a problem.
- **Power of 2-opt.** A good and quicker approximate might come in handy. Not making perfect the enemy of the good!
- **Obtained an understanding of the Genetic Algorithm.** In addition to understanding how it could be used to solve TSP, I am fascinated by its application in data driven control engineering.
- **Handling, organizing and accessing data** efficiently and effectively by leveraging the power of STLs.
- **Fine tuned my thought process** while addressing a problem. Has become more structured and systematic.

### Learnings - _Mackton Vishal Mendonca_

- Research on this problem revealed some non computer science related problems that are solved using TSP. (**Physarum polycephalum**, an amoeboid that adapts its morphology to create an efficient path between its food sources).
- Substituting cities and travelling distances for chip components and connecting circuits, the potential of TSP in chip design to integrate more and more transistors on a chip.
- Deeper understanding of object oriented programming for data extraction and usage of structures in C++.
- Understanding the importance of testing an application to not just provide a solution but also catch “not solutions”. To not only catch bugs but also to improve the quality of code.









