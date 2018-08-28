// ====================================================================
//
//	A network optimization algorithm using Fibonacci Heap
//
//	Written by: Max Winkler
//
// ====================================================================

#ifndef _FIBONACCI_HEAP_H_
#define _FIBONACCI_HEAP_H_

#include <stdio.h>
#include <vector> 

namespace FibonacciHeap
{
	enum State
	{
		LABELED, UNLABELED, SCANNED, UNSCANNED
	};

	class Node;

	class Edge
	{
	private:
	public:
		Node* tail;
		Node* head;
		double length;
		double delta;
	public:
		~Edge();
		Edge();
		Edge(Node* tail, Node* head, double length);
		Edge(const Edge& rhs);
		Edge& operator=(const Edge& rhs);
	};


	class Node
	{
	private:
	public:
		Node * parent;
		Node * leftSibling, * rightSibling;
		Node * children; 
		Node * pred;


		int data;
		double key;
		int rank;

		std::vector<Edge*> incomingEdges;
		std::vector<Edge*> outgoingEdges;

		int color;
		int index;

		State state;

	
	public:

		Node(int data, double key);
		Node();
		Node(const Node& rhs);
		Node& operator=(const Node& rhs);
		~Node();

		bool addChild(Node * node);
		bool addSibling(Node * node);
		bool remove();
		
		Node* leftMostSibling();
		Node* rightMostSibling();

		void addIncomingEdge(Edge * edge);
		void addOutgoingEdge(Edge * edge);

	};

	class FibonacciHeap
	{
	private:
		Node* rootListByRank[100];

		bool link(Node* root);	
		Node * minRoot;

	public:

		FibonacciHeap();
		FibonacciHeap(Node *root);

		~FibonacciHeap();

		bool isEmpty();
		bool insertVertex(Node * node);	
		void decreaseKey(double delta, Node* vertex);

		Node* findMin();
		Node* deleteMin();
	};
};

#endif
