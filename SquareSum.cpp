/*A program that is intended to use graph traversal techniques to
  find solutions to the square sum problem working up from 1 to
  whatever desired integer. It will output whether the square sum
  problem is solvable with the numbers available from 1 to that
  integer and print out the list that is a solution for that set of
  numbers.

  The Square Sum Problem: We need to take all integers from 1 to some
  integer and arrange in such an order that every integer's neighbors
  are numbers that add with it to make perfect squares, for example
  3, 1, 8 since 1 and 3 make 4, and 1 and 8 make 9. We solve this
  by creating a graph of the integers and connecting each number
  vertex to any other number vertex that adds to it to make a perfect
  square. */


#include <math.h>
#include <iostream>
#include <sstream>
#include <stack>
#include <vector>

using namespace std;

#define seqLen 75

void populateNeighbors(vector < vector < int > > &, const int);
void findSequence(vector < vector < int > > &, stack < int > &);

int main(void){

vector <vector < int > > SquareSums;	
stack <int> Traversal;

/*iterates through finding sequences all the way to seqLen which
  is the defined stop point*/
for(int Seq = 1; Seq <= seqLen ; Seq++){	
	SquareSums.push_back( vector < int >() );
	populateNeighbors(SquareSums, Seq);
	findSequence(SquareSums, Traversal);
	//The stack is dumped if no valid path exists
	if(Traversal.empty()){
		cerr << "No Solution for integers from 1 to " << Seq << endl;
	}
	//print out the sequence in the stack if a solution is found
	else{
		while(!Traversal.empty()){
			cerr << Traversal.top() << "  ";
			Traversal.pop();
		}
		cout << endl;
	}
}

//prints out the graph structure at the end into stderr
for(int j = 0; j < seqLen; j++){
	
	cerr << j+1 << ": " << endl;
	for(int i = 0; i < SquareSums[j].size(); i++){
	
		cerr << " " << SquareSums[j][i];

	}
	cerr << endl;
}
return 0;
}
/*Goes through the vector testing each index which represents the current node
  and tests it with each of the other existing nodes, if a match is found they add
  each other as neighbors*/
void populateNeighbors(vector < vector < int > > &squareSums, const int newInt){
	float intRoot = 0;
	for(int i = 1 ; i < newInt; i++){
		intRoot = sqrt(newInt + (i));
		if(intRoot - floor(intRoot) == 0){
			squareSums[i-1].push_back(newInt);
			squareSums[newInt - 1].push_back(i);
		}
	}
	return;
}

/*Not properly iterating through neighbors when pushing and popping things on and off the stack.
  Numbers that are not neighbors are being pushed together when they shouldn't be*/
bool sequenceRec(vector <vector < int > > &squareSums, stack < int > &Traversal, bool *visited, int curInt){
	if(visited[curInt-1]){
		return false;
	}
	//push the current node on the stack and mark it as visited
	Traversal.push(curInt);
	visited[curInt-1] = true;
	bool allVisited = true;
	for(int k = 0; k < squareSums.size(); k++){
		allVisited = allVisited & visited[k];	
	}
	//base case, all nodes are marked as visited so we've found a path through the graphh
	if(allVisited){
		return true;
	}
	//check all the neighbors of a node, which are its square sum pairs and move to the first unvisited node
	bool solution = false;
	for(int j = 0; j < squareSums[curInt-1].size(); j++){
		int newInt = squareSums[curInt-1][j];
		solution = sequenceRec(squareSums, Traversal, visited, newInt);
		if(solution){
				return true;
		}
	}
	//if we fall out of the loop all neighbors are visited so we have to go back and try a new path
	visited[Traversal.top()-1] = false;
	Traversal.pop();
	return false;
}

/*Iterate through the first node of a potential path until we find the first valid path through the graph
  when a solution is fouund the recursive sequence funtion will return true*/
void findSequence(vector <vector < int > > &squareSums, stack < int > &Traversal){
	int n = squareSums.size();
	bool visited[n];
	int curInt = 1;
	cout << "Progress " << n << " out of " << seqLen << "\r" << endl;
	for(int i = 0; i < n; i++, curInt++){
		for(int j = 0; j < n; j++){
			visited[j] = false;
		}
		bool solution = sequenceRec(squareSums, Traversal, visited, curInt);
		if(solution){
			cerr << "There is a solution from 1 to " << n << "\r" << endl;
			return;
		}
	}
	return;
}
