/**
 * A tile puzzle game in which the user must sort numerical tiles in a given order.
 * 
 * This program is a demonstration of some basic search methods used in classical AI.
 * When the solution space is enumerable, finite, and the combinatorial expansion can be kept
 * under control, search algorithms are a very effective method to find solutions to hard
 * problems. Games such as this are a perfect example for some of these techniques.
 * 
 * This very simple game is also a demonstration of the OneLoneCoder (javidx9, David Barr) Pixel
 * Game Engine, a simple and versatile engine to build algorithm visualizations, UI, and games.
 * 
 * Copyright (c) 2020 Christian Bruccoleri
 * LICENSE: MIT
 * 
 * The Pixel Game Engine is licensed under its own terms (see OLC-3 in the header file).
 * 
 * Version: 0.1
 * Date:    2020-08-01
 */
#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include <random>
#include <string>
#include <array>
#include <deque>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>
#include <chrono>


// Generic unsigned integer
typedef unsigned int uint_t;

// Basic size of the grid is twice the size of a character, which in PGE is 8 pixels
constexpr uint_t uFixedFontSize = 8;
constexpr uint_t uBlockSize = 2*uFixedFontSize;
constexpr uint_t uPuzzleSideSize = 3;
constexpr uint_t uNumTiles = uPuzzleSideSize * uPuzzleSideSize;
constexpr uint_t uScreenWidth = uBlockSize*uPuzzleSideSize;
constexpr uint_t uScreenHeight = uBlockSize*uPuzzleSideSize;


/**
 * A tile puzzle, like "8" or "15".
 * 
 * The user can move the blank tile in order to sort the tiles in the desired order.
 * User actions are added to a queue and executed whenever the status of the UI is updated.
 * With this approach, it is possible for either the user or another algorithm to enqueue actions
 * to be performed.
 */
class BlockPuzzle: public olc::PixelGameEngine
{
public:
	enum Actions {
		Nothing,
		MoveUp,
		MoveRight,
		MoveDown,
		MoveLeft
	};
	// Number of available actions
	static constexpr int iNumActions = 5;
	// Text to describe the actions
	static const std::string sActionText[iNumActions];

public:
	BlockPuzzle()
	{
		sAppName = "Block Puzzle";
		m_bShuffle = false;
	}


public:
	bool OnUserCreate() override
	{
		// Called once at the start, so create things here
		if (m_bShuffle) {
			std::random_shuffle(m_sTiles.begin(), m_sTiles.end());
		}
		m_bReplayMode = false;
		m_bGameOver   = IsEndState();
		m_posMovable  = FindMovable();
		return true;
	}

	void DrawBlocks()
	{
		constexpr int iBlockSize = static_cast<int>(uBlockSize);
		constexpr int iPuzzleSize = static_cast<int>(uPuzzleSideSize);
		constexpr int iTxtOff = 4;
		const olc::vi2d blockDim = {iBlockSize - 1, iBlockSize - 1};
		for (int i = 0; i < iPuzzleSize; ++i)
			for (int j = 0; j < iPuzzleSize; ++j) {
				char buf[] = " ";
				buf[0] = m_sTiles[ i*iPuzzleSize + j ];
				const olc::vi2d blockOrig = {j*iBlockSize, i*iBlockSize};
				if (buf[0] == ' ')
					FillRect(blockOrig, blockDim, olc::DARK_YELLOW);
				else {
					olc::Pixel tileTextColor;
					if (m_bGameOver)
						tileTextColor = olc::DARK_RED;
					else
						tileTextColor = olc::WHITE;
					FillRect(blockOrig, blockDim, olc::DARK_GREY);
					DrawRect(blockOrig, blockDim, olc::GREY);
					DrawString({blockOrig.x + iTxtOff, blockOrig.y + iTxtOff}, std::string(buf), tileTextColor);
				}
			}
	}
	
	bool CheckUserInput()
	{
		using namespace olc;
		if (GetKey(Key::ESCAPE).bPressed) {
			// Exit the application by returning false
			return false;
		}
		if (GetKey(Key::SPACE).bPressed) {
			m_bShuffle = true;
			OnUserCreate();
		}
		if (GetKey(Key::F2).bPressed) {
			auto solVector = DepthFirstSearch();
			m_actionListSol = std::move(solVector);
			//test_hashInsertion();
		}
		if (GetKey(Key::F3).bPressed) {
			// add the solution actions to the queue, in reverse order
			std::cout << "Adding solution to queue: " << m_actionListSol.size() << " nodes.\n";
			for (auto itActNode = m_actionListSol.rbegin(); 
				 itActNode != m_actionListSol.rend(); 
				 ++itActNode) {
				m_qActions.push_back((*itActNode)->action);
			}
			std::cout << "Queue length: " << m_qActions.size() << "\n";
			m_bReplayMode = true; // disable key strokes while replaying
		}
		if (! (m_bGameOver || m_bReplayMode)) {
			// Move the empty tile using the Arrow Keys
			if (GetKey(Key::UP).bPressed && IsValidAction(Actions::MoveUp)) {
				m_qActions.push_back(Actions::MoveUp);			
			}
			if (GetKey(Key::RIGHT).bPressed && IsValidAction(Actions::MoveRight)) {
				m_qActions.push_back(Actions::MoveRight);			
			}
			if (GetKey(Key::DOWN).bPressed && IsValidAction(Actions::MoveDown)) {
				m_qActions.push_back(Actions::MoveDown);			
			}
			if (GetKey(Key::LEFT).bPressed && IsValidAction(Actions::MoveLeft)) {
				m_qActions.push_back(Actions::MoveLeft);			
			}
		}
		return true;
	}

	/**
	 * Perform the actions in the actions queue.
	 */
	bool UpdateState()
	{
		if (! m_qActions.empty()) {
			DoAction(m_qActions.front());
			m_qActions.pop_front();
			if (m_bReplayMode) { // we are replaying actions
				using namespace std::chrono_literals;
				std::this_thread::sleep_for(500ms);
			}
			m_bGameOver = IsEndState();
		}
		return true;
	}

	/**
	 * Update the state of the game UI.
	 */
	bool OnUserUpdate(float fElapsedTime) override
	{
		bool bStatus = CheckUserInput();
		// allow exit on ESC pressed
		if (! bStatus)
			return bStatus;
		// perform actions
		UpdateState();

		// called once per frame
		DrawBlocks();	
		return true;
	}

private:
	// State of the game

	// Actions are queued as integers (see also Actions enumeration)
	typedef std::deque<Actions> 	QueueActions;

	// Forward declaration, see below.
	struct ActionNode;
	/// A (more readable) pointer to a node.
	typedef std::shared_ptr<ActionNode> ActionNodePtr;

	/// A vector of pointers to nodes
	typedef std::vector<ActionNodePtr> 	ActionPtrList;

	/// A queue of pointers to nodes
	typedef std::deque<ActionNodePtr>  	ActionPtrDeque;

	// defines the end-state of the game;
	static const std::string sEndState;
	// State of the tiles, by row.
	//std::string		m_sTiles = "12345678 ";
	std::string		m_sTiles = "1234567 8";
	// Position of the movable empty square.
	olc::vi2d 		m_posMovable = {-1, -1}; // set to invalid position initially
	// If true, the tiles are shuffled at startup.
	bool 			m_bShuffle;
	// If true, no more moves are processed and the game is over.
	bool 			m_bGameOver;
	// If true, enables replay mode.
	bool 			m_bReplayMode;
	// Queue containing the actions to be performed.
	QueueActions	m_qActions;
	// Solution list
	ActionPtrList 	m_actionListSol;

private:
	// Methods that manipulate the state


	[[nodiscard]] bool IsEndState() const
	{
		return IsEndState(m_sTiles);
	}

	
	/**
	 * Return true if the current state is the end state, false otherwise.
	 */
	[[nodiscard]] static bool IsEndState(const std::string& sBoard)
	{
		bool bEndState = true;
		size_t i = 0;
		while (i < sBoard.length() && (bEndState = (sBoard[i] == sEndState[i]))) i++;
		return bEndState;
	}

	[[nodiscard]] olc::vi2d FindMovable() const 
	{
		return FindMovable(m_sTiles);
	}

	/**
	 * Find the position of the ' ' cell within the tile string.
	 */
	[[nodiscard]] static olc::vi2d FindMovable(const std::string& sBoard) 
	{
		olc::vi2d pos = {-1, -1};
		for (size_t i = 0; i < sBoard.size(); i++)
		{
			if (sBoard[i] == ' ')
				pos = { int(i % uPuzzleSideSize), int(i / uPuzzleSideSize) };
		}
		return pos;		
	}


	[[nodiscard]] bool IsValidAction(Actions action) const
	{
		return IsValidAction(action, m_posMovable);
	}

	/**
	 * Check whether the action passed as argument results in a valid state of the board.
	 */
	[[nodiscard]] static bool IsValidAction(Actions action, const olc::vi2d& oldPos)
	{
		olc::vi2d newPos{oldPos};
		if (newPos.x == -1 || newPos.y == -1)
			return false;
		else {
			ApplyActionToPosition(action, newPos);
			const int maxSize = static_cast<int>(uPuzzleSideSize);
			return  0 <= newPos.x && newPos.x < maxSize && 
					0 <= newPos.y && newPos.y < maxSize;
		}
	}


	/**
	 * Move the blank tile in the new position.
	 */
	void DoAction(Actions action)
	{
		const olc::vi2d oldPos = m_posMovable;
		size_t oldOffset = CalcOffset(oldPos);
		ApplyActionToPosition(action, m_posMovable);
		size_t newOffset = CalcOffset(m_posMovable);
		std::swap(m_sTiles[oldOffset], m_sTiles[newOffset]);
	}

	/**
	 * Compute the offset of a board position into the representation string.
	 */
	[[nodiscard]] static inline int CalcOffset(const olc::vi2d& pos)
	{
		return pos.y * uPuzzleSideSize + pos.x;
	}

private:
	// Implementation of search algorithms


	/**
	 * Auxiliary node-structure used for searching.
	 */
	struct ActionNode {
		// Action performed
		Actions action;
		// State of the board after the action
		std::string sBoard;
		// Cost of reaching this configuration
		int cost;
		// Heuristic score of the board
		int score;
		// Pointer to the previous node
		std::shared_ptr<ActionNode> prev{nullptr};

		ActionNode() = default;

		explicit ActionNode(std::string _sBoard ): 
			sBoard{_sBoard}, 
			cost{0}, 
			score{0}, 
			action{Actions::Nothing} 
			{}

		ActionNode(Actions _act, const std::string& _board, int _c, int _score, ActionNodePtr _prev):
			action{_act},
			sBoard{_board},
			cost{_c},
			score{_score},
			prev{_prev}
			{}

		
		ActionNode(const ActionNode&) = default;

		ActionNode& operator=(const ActionNode&) = default;

		/**
		 * Apply the action of this node to its board representation.
		 * 
		 * This function is used in creating child nodes from parent nodes.
		 * Pre:
		 * 		- A valid board configuration.
		 * Post:
		 * 		- A valid board configuration after the application of the node's action.
		 * 		- `score` is updated to reflect the change.
		 */
		void ApplyAction()
		{
			olc::vi2d posBlank = FindMovable(sBoard);
			int iOldOffset = CalcOffset(posBlank);
			ApplyActionToPosition(action, posBlank);
			int iNewOffset = CalcOffset(posBlank);
			std::swap(sBoard[iNewOffset], sBoard[iOldOffset]);
			score = BlockPuzzle::ScoreBoard(sBoard);
		}

		/**
		 * Return true IFF the action argument is valid in the current state of this node.
		 * This function is used when creating child nodes of a parent node.
		 */
		[[nodiscard]] bool IsValidAction(Actions action)
		{
			if (sBoard.length() != uNumTiles)
				return false;
			else {
				olc::vi2d posEmptyTile = FindMovable(this->sBoard);
				ApplyActionToPosition(action, posEmptyTile);
				return (0 <= posEmptyTile.x && posEmptyTile.x < uPuzzleSideSize) &&
					   (0 <= posEmptyTile.y && posEmptyTile.y < uPuzzleSideSize);
			}
		}
	};


	/**
	 * Hash functor to be used with  std::unordered_set<>
	 */
	struct NodeHash
	{
		std::size_t operator()(ActionNodePtr node) const {
			if (node)
				return std::hash<std::string>{}(node->sBoard);
			else
				return std::hash<std::string>{}(std::string{});
		}
	};

	/**
	 * Equality functor to be used with std::unordered_set<>
	 */
	struct NodeEquality
	{
		bool operator()(const ActionNodePtr& node1, const ActionNodePtr& node2) const {
			return node1->sBoard == node2->sBoard;
		}
	};


	struct NodeScoreCompare
	{
		bool operator()(const ActionNodePtr& node1, const ActionNodePtr& node2) const {
			return node1->score < node2->score;
		}
	};

	/// A hash-map containing pointers to nodes
	typedef std::unordered_set<ActionNodePtr, NodeHash, NodeEquality>	ActionPtrHash;


	// Assign cost to each action: for now we just count each action as costing unity.
	[[nodiscard]] inline static int ActionCost(Actions action)
	{
		// Make sure the cost is never negative when changing any of these.
		switch (action) {
			case Actions::MoveUp:
			case Actions::MoveRight:
			case Actions::MoveDown:
			case Actions::MoveLeft:
				return 1;
			case Actions::Nothing:
				return 0;
			default:
				return 0;
		}
	}


	/**
	 * Calls ScoreBoard(const std::string& strBoard) with the current board status.
	 */
	[[nodiscard]] inline int ScoreBoard() const
	{
		return ScoreBoard(m_sTiles);
	}

	/**
	 * Assign a score to the board by counting the number of tiles in the right
	 * location. This is an example of Heuristic function.
	 */
	[[nodiscard]] static int ScoreBoard(const std::string& strBoard)
	{
		int sum = 0;
		for (int i = 0; i < strBoard.length(); i++)
			if (strBoard[i] == sEndState[i])
				sum++;
		return sum;
	}

	/**
	 * Create a vector of all the valid moves possible from the current state
	 */
	[[nodiscard]] ActionPtrList GetValidMoves(const ActionNodePtr& pCurNode)
	{
		ActionPtrList list;
		// iterate over all actions skipping the empty action
		for (size_t act = 1; act < iNumActions; act++) {
			auto action = Actions(act);
			if (pCurNode->IsValidAction(action))
				list.push_back(MakeChildNode(action, pCurNode));
		}
		return list;
	}


	/**
	 * Make a child node from a given action.
	 * 
	 * Assumes that the action is valid.
	 */
	[[nodiscard]] ActionNodePtr MakeChildNode(Actions action, const ActionNodePtr& pParentNode)
	{
		auto pNewNode = std::make_shared<ActionNode>(pParentNode->sBoard);
		nExpansions ++;
		pNewNode->action = action;
		//pNewNode->sBoard = pParentNode->sBoard;
		// apply the action, if valid, to the node's board representation; score also updated
		pNewNode->ApplyAction();
		// update the cost (score updated above)
		pNewNode->cost = pParentNode->cost + ActionCost(action);
		// link the node to its parent so that we can reconstruct solution paths
		pNewNode->prev = pParentNode;
		return pNewNode;
	}

	ActionPtrList DepthFirstSearch()
	{
		ActionPtrList listSolution;
		//ActionPtrHash hashVisited;
		// start recursion call
		ActionNode nodeGoal{Actions::Nothing, sEndState, 0, 0, nullptr};
		ActionNodePtr pRootNode = std::make_shared<ActionNode>();
		pRootNode->sBoard = m_sTiles;
		//DepthFirstSearchRec(nodeGoal, pRootNode, hashVisited, listSolution);
		for (int iDepth = 1; iDepth <= 105 && listSolution.size() == 0; iDepth++ ) {
			DepthFirstSearchIter(nodeGoal, pRootNode, listSolution, iDepth);
			std::cout << "D: " << iDepth << '\n';
		}
		std::cout << "Sol Length: " << listSolution.size() << '\n';
		std::cout << " Hash hits: " << nHashHits << '\n';
		std::cout << "  Num Exp.: " << nExpansions << std::endl;
		WriteNodeList(listSolution);
		return listSolution;
	}


	void test_hashInsertion()
	{
		ActionPtrHash hashTable;
		ActionNodePtr pNode1 = std::make_shared<ActionNode>();		
		pNode1->sBoard = "12345678 ";
		ActionNodePtr pNode2 = std::make_shared<ActionNode>();		
		pNode2->sBoard = "12345678 ";
		auto pair1 = hashTable.insert(pNode1);
		auto pair2 = hashTable.insert(pNode2);
		std::cout << "insert 1: " << pair1.second << '\n';
		std::cout << "insert 2: " << pair2.second << '\n';
	}

	/**
	 * Recursive Depth-First Search with Visited set algorithm.
	 */
	bool DepthFirstSearchRec(const ActionNode& nodeGoal, ActionNodePtr pCurNode, 
							 ActionPtrHash& hashVisited, ActionPtrList& listSolution)
	{
		if (pCurNode->sBoard == nodeGoal.sBoard) {
			// goal reached: build the solution backward
			listSolution.push_back(pCurNode);
			return true;
		}
		else {
			auto iterFoundNode = hashVisited.find(pCurNode);
			if (iterFoundNode == hashVisited.end()) {
				auto pair = hashVisited.insert(pCurNode);
				if (! pair.second) { // insertion failed
					// TODO: use a non-console logging mechanism
					std::cerr << "!! Warning: insertion of " << *pair.first << " failed.\n";
				}
				// expand the actions available in the current node, get a list of node-pointers
				// to the child action-nodes.
				auto validChildrenList = GetValidMoves(pCurNode);
				// sort them in descending order by Score
				std::sort(validChildrenList.begin(), validChildrenList.end(), NodeScoreCompare{});
				for (auto& childNode: validChildrenList) {
					// for each valid action, Call DepthFirst Recursively
					bool bSolFound = DepthFirstSearchRec(nodeGoal, childNode, hashVisited, listSolution);
					// if any branch reached the goal, add to the solution the current node and 
					// return success
					if (bSolFound) {
						listSolution.push_back(pCurNode);
						return true;
					}
					// otherwise, it is a dead-end: try next.
				}
			}
			else {
				nHashHits++;
			}
			// else either the node has already been visited or
			// no child branches reach a solution, skip it.
			return false;
		}			
	}

	/**
	 * Recursive Depth-First Search with Visited set algorithm.
	 * If the last argument, `iMaxDepth` is used, it will only search down to the specified depth. If this argument
	 * is negative, it will not apply a maximum depth limit. This is useful to call this function with the Iterative
	 * Depth-First Search paradigm.
	 */
	bool DepthFirstSearchIter(const ActionNode& nodeGoal, ActionNodePtr pCurNode, ActionPtrList& listSolution, int iMaxDepth=-1)
	{
		ActionPtrHash 	hashVisited;
		ActionPtrDeque	dequeExpanded;
		dequeExpanded.push_front(pCurNode);
		int iDepth = 1;
		while (dequeExpanded.size() > 0 && iDepth <= iMaxDepth)
		{
			pCurNode = dequeExpanded.front();
			dequeExpanded.pop_front();
			auto pair = hashVisited.insert(pCurNode);
			if (pair.second) { // the current node is not in the visited table
				if (pCurNode->sBoard == nodeGoal.sBoard) {
					// goal reached: build the solution backward
					listSolution.push_back(pCurNode);
					while (pCurNode->prev) {
						pCurNode = pCurNode->prev;
						listSolution.push_back(pCurNode);
					}
					return true;
				}
				else {
					// expand the actions available in the current node, get a list of node-pointers
					// to the child action-nodes.
					auto validChildrenList = GetValidMoves(pCurNode);
					// if iMaxDepth was -1, make sure you continue searching indefinitely
					if (validChildrenList.size() > 0) {
						// at every expansion, we descend one level down
						iDepth = iMaxDepth < 0 ? iDepth: iDepth + 1;
						// OPTIONAL: sort them in ascending order by Score
						// recall: in depth first you need FIFO order, so the best is last
						std::sort(validChildrenList.begin(), validChildrenList.end(), NodeScoreCompare{});
					}
					// enqueue at front if the node has not already been expanded before
					for (auto& childNode: validChildrenList) {
						auto iterFoundNode = hashVisited.find(childNode);
						if (iterFoundNode == hashVisited.end()) // not already expanded, add it
							dequeExpanded.push_front(childNode);
						else // skip it and keep tab
							nHashHits++;
					}
				}
			}
		}
		// no branch reaches a solution, clean up and exit
		for (auto& shrd_ptr: hashVisited)
			shrd_ptr->prev=nullptr; // make sure there are no "loops" of pointers
		return false;
	}

	// auxiliary variables to report algorithm performance
	size_t nHashHits = 0;
	size_t nExpansions = 0;

	/**
	 * Apply the specified action to the second argument.
	 */
	static void ApplyActionToPosition(Actions action, olc::vi2d& position) noexcept
	{
		switch (action) {
			case Actions::MoveUp:
				position.y -= 1;
				break;
			case Actions::MoveRight:
				position.x += 1;
				break;
			case Actions::MoveDown:
				position.y += 1;
				break;
			case Actions::MoveLeft:
				position.x -= 1;
				break;
			default:
				// do nothing
				break;
		}
	}


	/**
	 * Helper function to write a node to a stream.
	 * Default is output to stdout.
	 */
	static void WriteActionNode(const ActionNode& node, std::ostream& os = std::cout)
	{
		using namespace std;
		// action board cost score prev
		os  << "<addr: 0x"  << std::hex << std::uppercase << ulong(&node) << "> "
			<< std::dec << "[ action: " << sActionText[node.action] << ", board: '" 
			<< node.sBoard << "', cost: "  << node.cost << ", score: " << node.score 
			<< ", prev: 0x" << std::hex << std::uppercase << ulong(node.prev.get())
			<< std::dec;
	}

	/**
	 * Helper function to output a list of nodes into a stream.
	 * Default is output to stdout.
	 */
	static void WriteNodeList(const ActionPtrList& pNodeList, std::ostream& os = std::cout)
	{
		os << "--- Node ---\n";
		for (auto& pNode: pNodeList) {
			WriteActionNode(*pNode);
			os << '\n';
		}
		os << "\n";
	}
}; // End BlockPuzzle

// Static initializers
const std::string BlockPuzzle::sEndState = "12345678 ";
const std::string BlockPuzzle::sActionText[BlockPuzzle::iNumActions] = 	{
	"Nothing", 
	"Up", 
	"Right", 
	"Down", 
	"Left"
}; 


// Entry Point: MAIN
int main()
{
	BlockPuzzle demo;
	if (demo.Construct(uScreenWidth, uScreenHeight, 8, 8))
		demo.Start();

	return 0;
}
